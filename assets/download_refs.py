#!/usr/bin/env python3
"""Download a matched Ensembl reference set: genome FASTA, annotation GTF and
transcriptome (cDNA) FASTA.

Usage: download_refs.py <species> <source> <outdir> [release]

  species : Ensembl species name, lower-case underscore-separated (homo_sapiens)
  source  : 'ensembl' (only supported source)
  outdir  : directory to write the reference files into
  release : an Ensembl release number (e.g. 102) to pin for a reproducible
            download, or 'current' (default) for the latest release.
            A leading 'v' is accepted and stripped (v102 -> 102).

The three files must describe the same assembly or every downstream result is
silently wrong. The genome FASTA names the assembly; the GTF and the cDNA FASTA
are then selected by requiring that same assembly, and the script fails loudly
rather than guessing whenever a pattern does not match exactly one file.

Ensembl publishes several GTF variants per species directory, e.g.

    Homo_sapiens.GRCh38.116.abinitio.gtf.gz
    Homo_sapiens.GRCh38.116.chr.gtf.gz
    Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz
    Homo_sapiens.GRCh38.116.gtf.gz              <- the one we want

Picking 'chr_patch_hapl_scaff' pairs a patch/haplotype annotation with a
primary-assembly genome. Measured on release-116 human, that GTF describes 528
contigs against the FASTA's 70 -- 458 of them have no sequence at all -- and
carries 86,411 gene records instead of 78,941. STAR cannot place the records on
the missing contigs, so the index it builds does not match the annotation it was
given; the surplus genes can never receive a read; and duplicated gene symbols
jump from 484 to 3,384 (HLA-A alone appears 8 times instead of once), which
breaks every symbol-keyed downstream step. Gene *ids* stay unique -- ALT copies
get their own ENSG accessions -- so gene-level counts are not themselves
corrupted, but the reference is simply not the one the user asked for.

Two traps make the naive filename rules wrong, both verified against the live
Ensembl FTP:

  * The trailing number is not always the Ensembl release. Release-116 ships
    Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz, where 63 is an annotation
    version.
  * A release directory can hold GTFs for more than one assembly. Release-110
    carries both Drosophila_melanogaster.BDGP6.32.110.gtf.gz and
    ...BDGP6.46.110.gtf.gz while the genome FASTA is BDGP6.46.

Matching on the assembly handles both. Note also that pub/current_gtf/ no
longer exists (it 404s) even though pub/current_fasta/ does, so 'current' is
resolved to a concrete release number up front and every file is fetched from
the same release-<N>/ directory.
"""
import sys
import os
import re
import gzip
import time
import datetime

import requests

ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
ENSEMBL_REST = "https://rest.ensembl.org"

RETRIES = 4
BACKOFF = 5  # seconds, multiplied by attempt number


def log(log_file, msg):
    """Print to stdout and append to the run's download log."""
    print(msg, flush=True)
    with open(log_file, "a") as fh:
        fh.write(msg + "\n")


def die(log_file, msg):
    log(log_file, f"ERROR: {msg}")
    sys.exit(1)


def _get(url, **kwargs):
    """GET with retries and exponential-ish backoff."""
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            r = requests.get(url, timeout=60, **kwargs)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa: BLE001 - network errors are varied
            last = e
            if attempt < RETRIES:
                time.sleep(BACKOFF * attempt)
    raise RuntimeError(f"{url}: {last}")


def listing(url):
    """Fetch an FTP-over-HTTP directory listing as text."""
    return _get(url).text


def pick_one(pattern, text, what, url, log_file):
    """Return the single filename matching `pattern`, or fail loudly.

    Ambiguity here is never benign: it means Ensembl's layout changed and we no
    longer know which file we are downloading.
    """
    hits = sorted(set(re.findall(pattern, text)))
    if not hits:
        die(log_file, f"no {what} matching /{pattern}/ at {url}")
    if len(hits) > 1:
        die(log_file,
            f"{what} is ambiguous at {url}: {len(hits)} files match "
            f"/{pattern}/ -> {hits}. Refusing to guess.")
    return hits[0]


def verify_gzip(path, log_file):
    """Stream the whole gzip member to catch truncated / corrupt downloads."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 22):
                pass
    except Exception as e:  # noqa: BLE001
        die(log_file, f"{os.path.basename(path)} is not a valid gzip file: {e}")


def download_file(url, outfile, log_file):
    log(log_file, f"Downloading: {url}")
    for attempt in range(1, RETRIES + 1):
        try:
            with _get(url, stream=True) as r:
                expect = r.headers.get("Content-Length")
                with open(outfile, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            got = os.path.getsize(outfile)
            if expect is not None and int(expect) != got:
                raise IOError(f"truncated: expected {expect} bytes, got {got}")
            verify_gzip(outfile, log_file)
            log(log_file, f"  saved {outfile} ({got} bytes, gzip OK)")
            return
        except SystemExit:
            raise
        except Exception as e:  # noqa: BLE001
            log(log_file, f"  attempt {attempt}/{RETRIES} failed: {e}")
            if os.path.exists(outfile):
                os.remove(outfile)
            if attempt < RETRIES:
                time.sleep(BACKOFF * attempt)
    die(log_file, f"could not download {url} after {RETRIES} attempts")


def resolve_release(release, log_file):
    """Turn 'current' into a concrete release number.

    Everything is then fetched from an explicit release-<N>/ directory. This
    matters for two reasons:

      * Ensembl no longer serves a pub/current_gtf/ alias (it 404s), while
        pub/current_fasta/ still exists -- so a 'current' run cannot mix the two.
      * Pinning the number means the genome, GTF and cDNA all come from the same
        release directory, and the log records exactly which one, so a 'current'
        download stays interpretable after the next Ensembl release.
    """
    if release not in ("", "current"):
        return release

    try:
        r = _get(f"{ENSEMBL_REST}/info/data/?",
                 headers={"Content-Type": "application/json"})
        releases = r.json().get("releases") or []
        if releases:
            n = str(max(int(x) for x in releases))
            log(log_file, f"Resolved 'current' to Ensembl release-{n} (via REST).")
            return n
    except Exception as e:  # noqa: BLE001
        log(log_file, f"(Ensembl REST release lookup failed: {e}; falling back to FTP)")

    try:
        text = listing(f"{ENSEMBL_FTP}/")
        nums = [int(x) for x in re.findall(r'href="release-(\d+)/"', text)]
        if nums:
            n = str(max(nums))
            log(log_file, f"Resolved 'current' to Ensembl release-{n} (via FTP listing).")
            return n
    except Exception as e:  # noqa: BLE001
        log(log_file, f"(FTP release listing failed: {e})")

    die(log_file, "could not determine the current Ensembl release; "
                  "pass --download_release <N> to pin one explicitly")


def download_ensembl(species, outdir, log_file, release):
    release = resolve_release(release, log_file)
    log(log_file, f"Using Ensembl release-{release}.")

    base = f"{ENSEMBL_FTP}/release-{release}"
    dna_url = f"{base}/fasta/{species}/dna/"
    cdna_url = f"{base}/fasta/{species}/cdna/"
    gtf_url = f"{base}/gtf/{species}/"

    # Informational only; never blocks the download.
    try:
        r = _get(f"{ENSEMBL_REST}/info/assembly/{species}?",
                 headers={"Content-Type": "application/json"})
        log(log_file, f"Assembly for {species}: "
                      f"{r.json().get('default_coord_system_version')}")
    except Exception as e:  # noqa: BLE001
        log(log_file, f"(Could not query Ensembl REST for assembly info: {e})")

    sp = re.escape(species.capitalize())

    # ---- genome FASTA -------------------------------------------------------
    # Prefer the primary assembly (no patches, no haplotypes). Fall back to
    # toplevel for species that publish no whole-genome primary_assembly file
    # (Drosophila, for instance, only ships per-chromosome primary_assembly
    # files: '...dna.primary_assembly.2L.fa.gz' -- which this pattern correctly
    # does not match, because '.fa.gz' must follow 'primary_assembly' directly).
    #
    # Anchoring on the literal '.dna.' means the soft-masked (dna_sm) and
    # repeat-masked (dna_rm) variants can never be selected.
    log(log_file, "Finding genome FASTA...")
    dna_text = listing(dna_url)
    genome_pat = rf'href="({sp}\.[^"]*\.dna\.primary_assembly\.fa\.gz)"'
    if not re.search(genome_pat, dna_text):
        log(log_file, "  no whole-genome primary_assembly file; using toplevel")
        genome_pat = rf'href="({sp}\.[^"]*\.dna\.toplevel\.fa\.gz)"'
    genome_fn = pick_one(genome_pat, dna_text, "genome FASTA", dna_url, log_file)

    # The genome file names the assembly, and the assembly is what the GTF and
    # the transcriptome must agree with. Everything below is selected by it.
    m = re.match(rf'^{sp}\.(?P<asm>.+)\.dna\.(?:primary_assembly|toplevel)\.fa\.gz$',
                 genome_fn)
    if not m:
        die(log_file, f"cannot parse the assembly name out of '{genome_fn}'")
    assembly = m.group("asm")
    asm = re.escape(assembly)
    log(log_file, f"  assembly: {assembly}")
    download_file(dna_url + genome_fn, os.path.join(outdir, genome_fn), log_file)

    # ---- annotation GTF -----------------------------------------------------
    # Ensembl ships several GTFs per species directory, e.g. GRCh38 release-116:
    #   Homo_sapiens.GRCh38.116.abinitio.gtf.gz
    #   Homo_sapiens.GRCh38.116.chr.gtf.gz
    #   Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz
    #   Homo_sapiens.GRCh38.116.gtf.gz                    <- the canonical one
    #
    # Selecting on the assembly (not on the release number) is what keeps the
    # GTF consistent with the genome, and it is not optional:
    #
    #   * The trailing number is NOT always the Ensembl release. Release-116
    #     ships Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz -- 63 is an
    #     annotation version.
    #   * A release directory may hold GTFs for more than one assembly.
    #     Release-110 has both Drosophila_melanogaster.BDGP6.32.110.gtf.gz and
    #     ...BDGP6.46.110.gtf.gz, while the genome FASTA is BDGP6.46. Picking
    #     the wrong one pairs an annotation with a genome it does not describe.
    #
    # Requiring '<Species>.<assembly>.<digits>.gtf.gz' collapses every observed
    # case to exactly one file. The exclusion list is belt-and-braces.
    log(log_file, "Finding annotation GTF...")
    gtf_text = listing(gtf_url)
    hits = sorted(set(re.findall(rf'href="({sp}\.{asm}\.\d+\.gtf\.gz)"', gtf_text)))
    BAD = ("abinitio", "chr.gtf", "chr_patch_hapl_scaff")
    hits = [h for h in hits if not any(b in h for b in BAD)]
    if len(hits) > 1:
        # Same assembly, several annotation versions: prefer this release's.
        preferred = [h for h in hits if h.endswith(f".{release}.gtf.gz")]
        if len(preferred) == 1:
            log(log_file, f"  {len(hits)} annotation versions for {assembly}; "
                          f"taking the release-{release} one")
            hits = preferred
    if len(hits) != 1:
        die(log_file,
            f"expected exactly one canonical GTF for assembly '{assembly}' at "
            f"{gtf_url}, found {hits}. Refusing to guess -- an annotation that "
            "does not describe the genome FASTA (a patch/haplotype GTF, or a "
            "GTF for a different assembly) makes the index, the counts and the "
            "gene symbols disagree about what the reference is.")
    gtf_fn = hits[0]
    download_file(gtf_url + gtf_fn, os.path.join(outdir, gtf_fn), log_file)

    # ---- transcriptome (cDNA) FASTA ----------------------------------------
    # Needed by --build_indices for Salmon and Kallisto, and by
    # IsoformSwitchAnalyzeR. '.cdna.all.fa.gz' is the full transcriptome;
    # '.cdna.abinitio.fa.gz' is an ab-initio prediction set and must not be used.
    # Matched on the same assembly, for the same reason as the GTF.
    #
    # NOTE: this is cDNA only. A total-RNA / rRNA-depleted library should also
    # index the ncRNA FASTA (fasta/<species>/ncrna/), or lncRNAs go unquantified
    # and their reads get misassigned.
    log(log_file, "Finding transcriptome (cDNA) FASTA...")
    cdna_text = listing(cdna_url)
    cdna_pat = rf'href="({sp}\.{asm}\.cdna\.all\.fa\.gz)"'
    cdna_fn = pick_one(cdna_pat, cdna_text, f"cDNA FASTA for assembly '{assembly}'",
                       cdna_url, log_file)
    download_file(cdna_url + cdna_fn, os.path.join(outdir, cdna_fn), log_file)

    # ---- provenance ---------------------------------------------------------
    # Record the resolved release so a 'current' download stays interpretable
    # after Ensembl moves on.
    log(log_file, "")
    log(log_file, "Reference set (all three describe assembly "
                  f"'{assembly}' from release-{release}):")
    log(log_file, f"  ensembl release : {release}")
    log(log_file, f"  assembly        : {assembly}")
    log(log_file, f"  genome          : {genome_fn}")
    log(log_file, f"  annotation      : {gtf_fn}")
    log(log_file, f"  transcriptome   : {cdna_fn}")
    log(log_file, "")
    log(log_file, "Reuse this exact set with:  --download_release " + release)


def main():
    if len(sys.argv) < 4:
        print("Usage: download_refs.py <species> <source> <outdir> [release]")
        sys.exit(1)

    species = sys.argv[1].lower().replace(" ", "_")
    source = sys.argv[2].lower()
    outdir = sys.argv[3]
    release = sys.argv[4] if len(sys.argv) > 4 else "current"
    release = re.sub(r"^[vV]", "", release.strip()) or "current"

    if release not in ("", "current") and not release.isdigit():
        sys.exit(f"Invalid release '{release}': expected a number (e.g. 102) or 'current'")

    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, "download_log.txt")

    log(log_file, f"Download started at: {datetime.datetime.now()}")
    log(log_file, f"Species: {species}")
    log(log_file, f"Source: {source}")
    log(log_file, f"Ensembl release: {release}")

    if source == "ensembl":
        download_ensembl(species, outdir, log_file, release)
    elif source == "ncbi":
        die(log_file,
            "source 'ncbi' is not implemented. Reliable NCBI downloads need the "
            "ncbi-datasets CLI. Use --download_source ensembl, or supply "
            "--genome_fasta / --gtf / --transcript_fasta directly.")
    else:
        die(log_file, f"unknown source '{source}'. Supported: ensembl")

    log(log_file, f"Download complete at: {datetime.datetime.now()}")


if __name__ == "__main__":
    main()
