#!/usr/bin/env python3
"""
======================================================================
  Variant Calling Pipeline — Nili Ravi (or any WGS/WES dataset)
  Automates: QC → Trimming → Alignment → Variant Calling → Filtering
  Supports single-sample and multi-sample (sample sheet) modes.
======================================================================

Author  : Prerna
Version : 3.0.0

SAMPLE SHEET FORMAT (TSV or CSV)
---------------------------------
No header row required. Each line = one sample:

  Column 1 : path to forward reads (R1)
  Column 2 : path to reverse reads (R2)

Example (samples.tsv):
  /data/nili_ravi_1.fastq.gz  /data/nili_ravi_2.fastq.gz
  /data/murrah_1.fastq.gz     /data/murrah_2.fastq.gz

The sample name (used for output folders and file prefixes) is
auto-derived from the R1 filename by stripping the common paired-end
suffixes (_1.fastq.gz, _R1.fastq.gz, _1.fastq, _R1.fastq).

Lines starting with '#' are treated as comments and ignored.
Delimiter is auto-detected (tab or comma or whitespace).
"""

import argparse
import logging
import os
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────
# DEFAULTS
# ─────────────────────────────────────────────────────────────────────
DEFAULTS = {
    "threads":         8,
    "java_mem":        "4G",
    "trimmomatic_jar": "/home/apps/illumina/scratch/compile/Trimmomatic-0.39/trimmomatic-0.39.jar",
    "adapters":        "/home/apps/illumina/scratch/compile/Trimmomatic-0.39/TruSeq3-PE.fa",
    "picard_jar":      "/home/apps/pacbio/scratch/compile/picard/build/libs/picard.jar",
    "bwa":             "/home/apps/illumina/scratch/compile/bwa/bwa",
    "gatk":            "/home/apps/compilers/anaconda3/envs/exome/bin/gatk",
    "rgid":            "SAMPLE",
    "rglb":            "lib1",
    "rgpl":            "ILLUMINA",
    "rgpu":            "unit1",
    "rgsm":            "SAMPLE",
    "sliding_window":  "4:20",
    "min_len":         36,
    "phred":           "33",
    "qd_snp":          2.0,
    "fs_snp":          60.0,
    "sor_snp":         3.0,
    "mq_snp":          40.0,
}

BANNER = r"""
╔══════════════════════════════════════════════════════════════════╗
║         W G S  V A R I A N T  C A L L I N G  P I P E L I N E   ║
║         BWA-MEM  ▸  Picard  ▸  GATK HaplotypeCaller            ║
╚══════════════════════════════════════════════════════════════════╝
"""

# RG columns that can be overridden per-sample in the sheet
PER_SAMPLE_RG_COLS = ("rgid", "rglb", "rgpl", "rgpu", "rgsm")


# ─────────────────────────────────────────────────────────────────────
# SAMPLE SHEET PARSER
# ─────────────────────────────────────────────────────────────────────
def derive_sample_id(r1_path: str) -> str:
    """Auto-derive a sample name from an R1 FASTQ filename."""
    name = Path(r1_path).name
    for suffix in ("_1.fastq.gz", "_1.fastq", "_R1.fastq.gz", "_R1.fastq",
                   "_1.fq.gz", "_1.fq", "_R1.fq.gz", "_R1.fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    # Fallback: strip all extensions
    return Path(name.split(".")[0]).name


def parse_sample_sheet(path: str) -> list[dict]:
    """
    Read a headerless TSV, CSV, or whitespace-delimited sample sheet.

    Each non-comment line must have exactly 2 columns:
      column 1 → R1 path (forward reads)
      column 2 → R2 path (reverse reads)

    The sample_id is auto-derived from the R1 filename.
    """
    p = Path(path)
    if not p.exists():
        sys.exit(f"[ERROR] Sample sheet not found: {path}")

    lines = [l for l in p.read_text().splitlines()
             if l.strip() and not l.startswith("#")]
    if not lines:
        sys.exit(f"[ERROR] Sample sheet is empty: {path}")

    # Auto-detect delimiter: tab → comma → whitespace
    first = lines[0]
    if "\t" in first:
        delimiter = "\t"
    elif "," in first:
        delimiter = ","
    else:
        delimiter = None           # split() handles any whitespace

    samples = []
    for i, line in enumerate(lines, 1):
        cols = line.split(delimiter) if delimiter else line.split()
        cols = [c.strip() for c in cols if c.strip()]

        if len(cols) < 2:
            sys.exit(
                f"[ERROR] Sample sheet line {i} has fewer than 2 columns:\n"
                f"        {line!r}\n"
                f"        Expected: <read1_path>  <read2_path>"
            )

        r1, r2 = cols[0], cols[1]
        sid = derive_sample_id(r1)
        samples.append({"sample_id": sid, "read1": r1, "read2": r2})

    # Warn about duplicate sample IDs
    seen: dict[str, int] = {}
    for s in samples:
        sid = s["sample_id"]
        seen[sid] = seen.get(sid, 0) + 1
    dupes = [sid for sid, cnt in seen.items() if cnt > 1]
    if dupes:
        print(f"[WARNING] Duplicate sample IDs detected: {dupes}. "
              f"Output folders will overwrite each other!", file=sys.stderr)

    return samples


# ─────────────────────────────────────────────────────────────────────
# ARGUMENT PARSER
# ─────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant_pipeline.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input modes (choose one)
------------------------
  MODE A — Single sample via direct arguments:
    python variant_pipeline.py \\
        -1 nili_ravi_1.fastq.gz -2 nili_ravi_2.fastq.gz \\
        -g ref.fna -o /results/

  MODE B — Multiple samples via a sample sheet (NO header, NO sample_id):
    python variant_pipeline.py \\
        --sample-sheet samples.tsv \\
        -g ref.fna -o /results/

    samples.tsv (tab/space/comma separated, col1=R1, col2=R2):
        /data/nili_ravi_1.fastq.gz  /data/nili_ravi_2.fastq.gz
        /data/murrah_1.fastq.gz     /data/murrah_2.fastq.gz

    Sample names are auto-derived from R1 filenames
    (nili_ravi_1.fastq.gz → nili_ravi, etc.)

  MODE C — Multiple samples directly on the command line (R1:R2 pairs):
    python variant_pipeline.py \\
        --samples /data/nili_ravi_1.fastq.gz:/data/nili_ravi_2.fastq.gz \\
                  /data/murrah_1.fastq.gz:/data/murrah_2.fastq.gz \\
        -g ref.fna -o /results/

Other examples
--------------
  # 16 threads, skip QC:
    python variant_pipeline.py --sample-sheet samples.tsv -g ref.fna -o /results/ \\
        --threads 16 --skip-qc

  # Dry run across all samples:
    python variant_pipeline.py --sample-sheet samples.tsv -g ref.fna -o /results/ --dry-run
""",
    )

    # ── Input mode ───────────────────────────────────────────────────
    inp = parser.add_argument_group(
        "Input mode  (provide ONE of: -1/-2, --sample-sheet, or --samples)"
    )
    inp.add_argument("-1", "--input1", metavar="FASTQ_R1",
                     help="Forward reads for a SINGLE sample")
    inp.add_argument("-2", "--input2", metavar="FASTQ_R2",
                     help="Reverse reads for a SINGLE sample")
    inp.add_argument("--sample-sheet", metavar="FILE",
                     help=(
                         "Headerless TSV/CSV/whitespace file — "
                         "column 1: R1 path, column 2: R2 path. "
                         "Sample name auto-derived from R1 filename."
                     ))
    inp.add_argument(
        "--samples", nargs="+", metavar="R1:R2",
        help=(
            "One or more R1:R2 path pairs separated by colon  "
            "(e.g. /data/nr_1.fq.gz:/data/nr_2.fq.gz). "
            "Sample name auto-derived from R1 filename."
        ),
    )

    # ── Required (always needed) ─────────────────────────────────────
    req = parser.add_argument_group("Required")
    req.add_argument("-g", "--genome", required=True, metavar="FASTA",
                     help="Reference genome FASTA (.fna / .fa / .fasta)")
    req.add_argument("-o", "--outdir", required=True, metavar="DIR",
                     help="Root output directory (created if absent)")

    # ── Tool paths ───────────────────────────────────────────────────
    tools = parser.add_argument_group("Tool paths (override defaults)")
    tools.add_argument("--bwa",         default=DEFAULTS["bwa"],
                       help=f"BWA binary  [default: {DEFAULTS['bwa']}]")
    tools.add_argument("--trimmomatic", default=DEFAULTS["trimmomatic_jar"], metavar="JAR",
                       help=f"Trimmomatic JAR  [default: {DEFAULTS['trimmomatic_jar']}]")
    tools.add_argument("--adapters",    default=DEFAULTS["adapters"],
                       help=f"Adapter FASTA  [default: {DEFAULTS['adapters']}]")
    tools.add_argument("--picard",      default=DEFAULTS["picard_jar"], metavar="JAR",
                       help=f"Picard JAR  [default: {DEFAULTS['picard_jar']}]")
    tools.add_argument("--gatk",        default=DEFAULTS["gatk"],
                       help=f"GATK executable  [default: {DEFAULTS['gatk']}]")

    # ── Performance ──────────────────────────────────────────────────
    perf = parser.add_argument_group("Performance")
    perf.add_argument("-t", "--threads", type=int, default=DEFAULTS["threads"],
                      help=(
                          f"Threads used by ALL parallel steps "
                          f"(FastQC, Trimmomatic, BWA-MEM, HaplotypeCaller)  "
                          f"[default: {DEFAULTS['threads']}]"
                      ))
    perf.add_argument("--java-mem", default=DEFAULTS["java_mem"],
                      help=f"Java heap for Picard/GATK  [default: {DEFAULTS['java_mem']}]")

    # ── Default read group metadata ──────────────────────────────────
    rg = parser.add_argument_group(
        "Default read group metadata  "
        "(per-sample values in the sheet override these)"
    )
    rg.add_argument("--rgid", default=DEFAULTS["rgid"],
                    help="Read group ID    [default: same as sample_id]")
    rg.add_argument("--rglb", default=DEFAULTS["rglb"],
                    help=f"Library          [default: {DEFAULTS['rglb']}]")
    rg.add_argument("--rgpl", default=DEFAULTS["rgpl"],
                    help=f"Platform         [default: {DEFAULTS['rgpl']}]")
    rg.add_argument("--rgpu", default=DEFAULTS["rgpu"],
                    help=f"Platform unit    [default: {DEFAULTS['rgpu']}]")
    rg.add_argument("--rgsm", default=DEFAULTS["rgsm"],
                    help="Sample name      [default: same as sample_id]")

    # ── Trimmomatic parameters ───────────────────────────────────────
    trim = parser.add_argument_group("Trimmomatic parameters")
    trim.add_argument("--sliding-window", default=DEFAULTS["sliding_window"],
                      help=f"SLIDINGWINDOW  [default: {DEFAULTS['sliding_window']}]")
    trim.add_argument("--min-len", type=int, default=DEFAULTS["min_len"],
                      help=f"MINLEN  [default: {DEFAULTS['min_len']}]")
    trim.add_argument("--phred", default=DEFAULTS["phred"], choices=["33", "64"],
                      help=f"Phred encoding  [default: {DEFAULTS['phred']}]")

    # ── GATK SNP hard-filter thresholds ─────────────────────────────
    filt = parser.add_argument_group("GATK SNP hard-filter thresholds")
    filt.add_argument("--qd-snp",  type=float, default=DEFAULTS["qd_snp"],
                      help=f"Filter SNPs with QD <  [default: {DEFAULTS['qd_snp']}]")
    filt.add_argument("--fs-snp",  type=float, default=DEFAULTS["fs_snp"],
                      help=f"Filter SNPs with FS >  [default: {DEFAULTS['fs_snp']}]")
    filt.add_argument("--sor-snp", type=float, default=DEFAULTS["sor_snp"],
                      help=f"Filter SNPs with SOR > [default: {DEFAULTS['sor_snp']}]")
    filt.add_argument("--mq-snp",  type=float, default=DEFAULTS["mq_snp"],
                      help=f"Filter SNPs with MQ <  [default: {DEFAULTS['mq_snp']}]")

    # ── Pipeline control ─────────────────────────────────────────────
    ctrl = parser.add_argument_group("Pipeline control")
    ctrl.add_argument("--skip-qc",            action="store_true",
                      help="Skip FastQC and seqkit QC")
    ctrl.add_argument("--skip-trim",          action="store_true",
                      help="Skip Trimmomatic (use raw reads for alignment)")
    ctrl.add_argument("--skip-alignment",     action="store_true",
                      help="Skip BWA alignment (assumes BAM already in outdir)")
    ctrl.add_argument("--skip-vc",            action="store_true",
                      help="Skip variant calling")
    ctrl.add_argument("--dry-run",            action="store_true",
                      help="Log commands without executing them")
    ctrl.add_argument("--keep-intermediates", action="store_true",
                      help="Keep SAM, unsorted BAM, RG BAM intermediates")
    ctrl.add_argument("--stop-on-error",      action="store_true",
                      help="Abort the whole run if any sample fails (default: continue)")

    # ── Misc ─────────────────────────────────────────────────────────
    misc = parser.add_argument_group("Miscellaneous")
    misc.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                      help="Verbosity  [default: INFO]")
    misc.add_argument("-v", "--version", action="version", version="%(prog)s 3.0.0")

    return parser


# ─────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────
def setup_master_logging(log_dir: Path, level: str) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s  [%(levelname)s]  %(message)s"
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    # Console
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(getattr(logging, level))
    sh.setFormatter(logging.Formatter(fmt))
    root.addHandler(sh)
    # Master file
    fh = logging.FileHandler(log_dir / "pipeline_master.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(fmt))
    root.addHandler(fh)
    return logging.getLogger("pipeline")


def step_logger(log_path: Path) -> logging.Logger:
    """Child logger that writes full DEBUG output to its own step log file."""
    name   = log_path.stem
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = True
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s  [%(levelname)s]  %(message)s"))
    logger.addHandler(fh)
    return logger


# ─────────────────────────────────────────────────────────────────────
# COMMAND RUNNER
# ─────────────────────────────────────────────────────────────────────
def run(cmd: str, logger: logging.Logger, dry_run: bool = False) -> None:
    logger.info("CMD: %s", cmd)
    if dry_run:
        return
    result = subprocess.run(
        cmd, shell=True, executable="/bin/bash",
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    for line in result.stdout.splitlines():
        logger.debug("  %s", line)
    if result.returncode != 0:
        logger.error("Command failed (exit %d)", result.returncode)
        raise RuntimeError(f"Step failed — see log: {log_path_from_logger(logger)}")


def log_path_from_logger(logger: logging.Logger) -> str:
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            return h.baseFilename
    return "<unknown>"


# ─────────────────────────────────────────────────────────────────────
# BUILD PER-SAMPLE ARGS NAMESPACE
# ─────────────────────────────────────────────────────────────────────
def sample_args(global_args, sample: dict):
    """
    Return a copy of args with sample-specific RG fields applied.
    Per-sample columns (rgid, rglb, rgpl, rgpu, rgsm) in the sheet
    override the global CLI defaults.
    """
    sa = deepcopy(global_args)
    sid = sample["sample_id"]
    # Apply per-sample overrides from sheet; fall back to sample_id for rgid/rgsm
    sa.rgid = sample.get("rgid") or sid
    sa.rglb = sample.get("rglb") or global_args.rglb
    sa.rgpl = sample.get("rgpl") or global_args.rgpl
    sa.rgpu = sample.get("rgpu") or global_args.rgpu
    sa.rgsm = sample.get("rgsm") or sid
    return sa


# ─────────────────────────────────────────────────────────────────────
# PIPELINE STEPS
# ─────────────────────────────────────────────────────────────────────

def step_qc_raw(args, dirs, r1, r2, master):
    master.info("  ── STEP 01 : QC on raw reads")
    log = step_logger(dirs["logs"] / "step01_qc_raw.log")
    run(f"seqkit stats {r1} {r2} > {dirs['qc_raw']}/seqkit_stats.txt", log, args.dry_run)
    run(f"fastqc -t {args.threads} {r1} {r2} -o {dirs['qc_raw']}", log, args.dry_run)
    master.info("     log → %s", dirs['logs'] / 'step01_qc_raw.log')


def step_trim(args, dirs, base, r1, r2, master):
    master.info("  ── STEP 02 : Trimmomatic adapter trimming")
    log = step_logger(dirs["logs"] / "step02_trimmomatic.log")
    p1 = f"{dirs['trim']}/{base}_1_paired.fastq.gz"
    u1 = f"{dirs['trim']}/{base}_1_unpaired.fastq.gz"
    p2 = f"{dirs['trim']}/{base}_2_paired.fastq.gz"
    u2 = f"{dirs['trim']}/{base}_2_unpaired.fastq.gz"
    run(
        f"java -jar {args.trimmomatic} PE "
        f"-threads {args.threads} -phred{args.phred} "
        f"{r1} {r2} {p1} {u1} {p2} {u2} "
        f"ILLUMINACLIP:{args.adapters}:2:30:10 "
        f"SLIDINGWINDOW:{args.sliding_window} MINLEN:{args.min_len}",
        log, args.dry_run,
    )
    master.info("     log → %s", dirs['logs'] / 'step02_trimmomatic.log')
    return p1, p2


def step_qc_trimmed(args, dirs, p1, p2, master):
    master.info("  ── STEP 03 : QC on trimmed reads")
    log = step_logger(dirs["logs"] / "step03_qc_trimmed.log")
    run(f"seqkit stats {p1} {p2} > {dirs['qc_trim']}/seqkit_stats.txt", log, args.dry_run)
    run(f"fastqc -t {args.threads} {p1} {p2} -o {dirs['qc_trim']}", log, args.dry_run)
    master.info("     log → %s", dirs['logs'] / 'step03_qc_trimmed.log')


def step_align(args, dirs, base, r1, r2, master):
    master.info("  ── STEP 04 : BWA-MEM alignment")
    log = step_logger(dirs["logs"] / "step04_bwa_mem.log")
    sam = f"{dirs['aln']}/{base}.sam"
    run(f"{args.bwa} mem -t {args.threads} {args.genome} {r1} {r2} -o {sam}",
        log, args.dry_run)
    master.info("     log → %s", dirs['logs'] / 'step04_bwa_mem.log')
    return sam


def step_sort(args, dirs, base, sam, master):
    master.info("  ── STEP 05 : Picard SortSam")
    log = step_logger(dirs["logs"] / "step05_picard_sort.log")
    sorted_bam = f"{dirs['aln']}/{base}.sorted.bam"
    tmp = f"{dirs['aln']}/tmp"; os.makedirs(tmp, exist_ok=True)
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} SortSam "
        f"I={sam} SORT_ORDER=coordinate O={sorted_bam} TMP_DIR={tmp}",
        log, args.dry_run,
    )
    if not args.keep_intermediates and not args.dry_run:
        Path(sam).unlink(missing_ok=True)
    master.info("     log → %s", dirs['logs'] / 'step05_picard_sort.log')
    return sorted_bam


def step_add_rg(args, dirs, base, sorted_bam, master):
    master.info("  ── STEP 06 : Picard AddOrReplaceReadGroups")
    log = step_logger(dirs["logs"] / "step06_picard_rg.log")
    rg_bam = f"{dirs['aln']}/{base}.rg.bam"
    tmp = f"{dirs['aln']}/tmp"
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} AddOrReplaceReadGroups "
        f"I={sorted_bam} O={rg_bam} "
        f"RGID={args.rgid} RGLB={args.rglb} RGPL={args.rgpl} "
        f"RGPU={args.rgpu} RGSM={args.rgsm} TMP_DIR={tmp}",
        log, args.dry_run,
    )
    if not args.keep_intermediates and not args.dry_run:
        Path(sorted_bam).unlink(missing_ok=True)
    master.info("     log → %s", dirs['logs'] / 'step06_picard_rg.log')
    return rg_bam


def step_mark_dups(args, dirs, base, rg_bam, master):
    master.info("  ── STEP 07 : Picard MarkDuplicates + samtools index")
    log = step_logger(dirs["logs"] / "step07_picard_markdup.log")
    dedup_bam = f"{dirs['aln']}/{base}.bam"
    metrics   = f"{dirs['aln']}/{base}.dup_metrics.txt"
    tmp = f"{dirs['aln']}/tmp"
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} MarkDuplicates "
        f"I={rg_bam} O={dedup_bam} M={metrics} "
        f"REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR={tmp}",
        log, args.dry_run,
    )
    run(f"samtools index {dedup_bam}", log, args.dry_run)
    if not args.keep_intermediates and not args.dry_run:
        Path(rg_bam).unlink(missing_ok=True)
    master.info("     log → %s", dirs['logs'] / 'step07_picard_markdup.log')
    return dedup_bam


def step_haplotype_caller(args, dirs, base, bam, master):
    master.info("  ── STEP 08 : GATK HaplotypeCaller")
    log = step_logger(dirs["logs"] / "step08_haplotypecaller.log")
    gvcf = f"{dirs['vc']}/{base}.g.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" HaplotypeCaller "
        f"-R {args.genome} -I {bam} -O {gvcf} "
        f"-ERC GVCF --native-pair-hmm-threads {args.threads}",
        log, args.dry_run,
    )
    master.info("     log → %s", dirs['logs'] / 'step08_haplotypecaller.log')
    return gvcf


def step_genotype_gvcfs(args, dirs, base, gvcf, master):
    master.info("  ── STEP 09 : GATK GenotypeGVCFs")
    log = step_logger(dirs["logs"] / "step09_genotypegvcfs.log")
    vcf = f"{dirs['vc']}/{base}.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" GenotypeGVCFs "
        f"-R {args.genome} -V {gvcf} -O {vcf}",
        log, args.dry_run,
    )
    master.info("     log → %s", dirs['logs'] / 'step09_genotypegvcfs.log')
    return vcf


def step_select_snps(args, dirs, base, vcf, master):
    """STEP 10 : Extract only SNPs from the genotyped VCF using SelectVariants."""
    master.info("  ── STEP 10 : GATK SelectVariants — extract SNPs")
    log = step_logger(dirs["logs"] / "step10_select_snps.log")
    snps_vcf = f"{dirs['vc']}/{base}_SNPs_raw.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" SelectVariants "
        f"-R {args.genome} "
        f"-V {vcf} "
        f"--select-type-to-include SNP "
        f"-O {snps_vcf}",
        log, args.dry_run,
    )
    master.info("     output → %s", snps_vcf)
    master.info("     log    → %s", dirs['logs'] / 'step10_select_snps.log')
    return snps_vcf


def step_filter_snps(args, dirs, base, snps_vcf, master):
    """STEP 11 : Apply hard filters to the raw SNP VCF using VariantFiltration."""
    master.info("  ── STEP 11 : GATK VariantFiltration — hard-filter SNPs")
    log = step_logger(dirs["logs"] / "step11_variantfiltration_snps.log")
    filt_vcf = f"{dirs['vc']}/{base}_SNPs_filtered.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" VariantFiltration "
        f"-R {args.genome} "
        f"-V {snps_vcf} "
        f"-O {filt_vcf} "
        f"--filter-name \"QD_lt_{args.qd_snp}\"   --filter-expression \"QD < {args.qd_snp}\" "
        f"--filter-name \"FS_gt_{args.fs_snp}\"   --filter-expression \"FS > {args.fs_snp}\" "
        f"--filter-name \"SOR_gt_{args.sor_snp}\" --filter-expression \"SOR > {args.sor_snp}\" "
        f"--filter-name \"MQ_lt_{args.mq_snp}\"   --filter-expression \"MQ < {args.mq_snp}\"",
        log, args.dry_run,
    )
    master.info("     output → %s", filt_vcf)
    master.info("     log    → %s", dirs['logs'] / 'step11_variantfiltration_snps.log')
    return filt_vcf


def step_pass_snps(args, dirs, base, filt_vcf, master):
    """STEP 12 : Keep only PASS variants from the filtered SNP VCF using SelectVariants."""
    master.info("  ── STEP 12 : GATK SelectVariants — extract PASS SNPs")
    log = step_logger(dirs["logs"] / "step12_pass_snps.log")
    pass_vcf = f"{dirs['vc']}/{base}_SNPs_final_PASS.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" SelectVariants "
        f"-V {filt_vcf} "
        f"--exclude-filtered "
        f"-O {pass_vcf}",
        log, args.dry_run,
    )
    master.info("     output → %s", pass_vcf)
    master.info("     log    → %s", dirs['logs'] / 'step12_pass_snps.log')
    return pass_vcf


# ─────────────────────────────────────────────────────────────────────
# SINGLE SAMPLE RUNNER
# ─────────────────────────────────────────────────────────────────────
def run_sample(args, sample: dict, master: logging.Logger) -> bool:
    """
    Run the full pipeline for one sample.
    Returns True on success, False on failure.
    args has already been patched with per-sample RG fields via sample_args().
    """
    base = sample["sample_id"]
    r1   = sample["read1"]
    r2   = sample["read2"]

    master.info("")
    master.info("╔══ SAMPLE : %s ══", base)
    master.info("  Read1 : %s", r1)
    master.info("  Read2 : %s", r2)

    # Validate files
    for label, path in [("R1", r1), ("R2", r2)]:
        if not Path(path).exists() and not args.dry_run:
            master.error("  %s file not found: %s  — skipping sample %s", label, path, base)
            return False

    # Create per-sample output directory tree
    root_out = Path(args.outdir) / base
    dirs = {
        "root":    root_out,
        "logs":    root_out / "logs",
        "qc_raw":  root_out / "qc_raw",
        "trim":    root_out / "trimmed_reads",
        "qc_trim": root_out / "qc_trimmed",
        "aln":     root_out / "alignment",
        "vc":      root_out / "variant_calling",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    master.info("  Output  : %s", root_out)
    master.info("  Step logs : %s/", dirs['logs'])

    t0 = datetime.now()
    try:
        # QC raw
        if not args.skip_qc:
            step_qc_raw(args, dirs, r1, r2, master)
        else:
            master.info("  Skipping raw QC (--skip-qc)")

        # Trim + trimmed QC
        if not args.skip_trim:
            p1, p2 = step_trim(args, dirs, base, r1, r2, master)
            if not args.skip_qc:
                step_qc_trimmed(args, dirs, p1, p2, master)
        else:
            master.info("  Skipping trimming (--skip-trim)")
            p1, p2 = r1, r2

        # Alignment
        if not args.skip_alignment:
            sam        = step_align(args, dirs, base, p1, p2, master)
            sorted_bam = step_sort(args, dirs, base, sam, master)
            rg_bam     = step_add_rg(args, dirs, base, sorted_bam, master)
            final_bam  = step_mark_dups(args, dirs, base, rg_bam, master)
        else:
            final_bam = f"{dirs['aln']}/{base}.bam"
            master.info("  Skipping alignment (--skip-alignment); expecting: %s", final_bam)

        # Variant calling
        if not args.skip_vc:
            gvcf     = step_haplotype_caller(args, dirs, base, final_bam, master)
            vcf      = step_genotype_gvcfs(args, dirs, base, gvcf, master)
            snps_vcf = step_select_snps(args, dirs, base, vcf, master)
            filt_vcf = step_filter_snps(args, dirs, base, snps_vcf, master)
            pass_vcf = step_pass_snps(args, dirs, base, filt_vcf, master)
            master.info("  Final PASS SNPs : %s", pass_vcf)
        else:
            master.info("  Skipping variant calling (--skip-vc)")

    except RuntimeError as exc:
        master.error("  FAILED : %s  (%s)", base, exc)
        elapsed = datetime.now() - t0
        return False, elapsed

    elapsed = datetime.now() - t0
    master.info("╚══ DONE  : %s  (elapsed: %s)", base, str(elapsed).split(".")[0])
    return True, elapsed


# ─────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────
def main() -> None:
    print(BANNER)
    parser = build_parser()
    args   = parser.parse_args()

    # ── Resolve sample list ──────────────────────────────────────────
    samples: list[dict] = []

    n_modes = sum([
        bool(args.input1 or args.input2),
        bool(args.sample_sheet),
        bool(args.samples),
    ])
    if n_modes == 0:
        parser.error(
            "Provide at least one input mode:\n"
            "  -1/-2  (single sample)\n"
            "  --sample-sheet FILE  (headerless file: col1=R1, col2=R2)\n"
            "  --samples R1:R2 ...  (inline pairs)"
        )
    if n_modes > 1:
        parser.error("Use only ONE input mode: -1/-2, --sample-sheet, or --samples")

    # Mode A: single sample via -1 / -2
    if args.input1 or args.input2:
        if not args.input1 or not args.input2:
            parser.error("Both -1/--input1 and -2/--input2 are required together")
        sid = derive_sample_id(args.input1)
        samples = [{"sample_id": sid, "read1": args.input1, "read2": args.input2}]

    # Mode B: sample sheet (headerless, col1=R1, col2=R2)
    elif args.sample_sheet:
        samples = parse_sample_sheet(args.sample_sheet)

    # Mode C: inline --samples R1:R2 ...
    elif args.samples:
        for token in args.samples:
            parts = token.split(":")
            if len(parts) != 2:
                parser.error(
                    f"--samples format must be R1_path:R2_path, got: {token!r}\n"
                    f"  Example: /data/nili_ravi_1.fastq.gz:/data/nili_ravi_2.fastq.gz"
                )
            r1, r2 = parts[0].strip(), parts[1].strip()
            samples.append({"sample_id": derive_sample_id(r1), "read1": r1, "read2": r2})

    if not samples:
        sys.exit("[ERROR] No samples found — check your input.")

    # ── Setup master logger (global log dir) ─────────────────────────
    global_log_dir = Path(args.outdir) / "_pipeline_logs"
    master = setup_master_logging(global_log_dir, args.log_level)

    master.info("Pipeline v3.0.0 started")
    master.info("Reference genome   : %s", args.genome)
    master.info("Output root        : %s", args.outdir)
    master.info("Threads (all steps): %d", args.threads)
    master.info("Java heap          : %s", args.java_mem)
    master.info("Total samples      : %d", len(samples))
    for i, s in enumerate(samples, 1):
        master.info("  [%d] %s  |  R1: %s  |  R2: %s",
                    i, s["sample_id"], s["read1"], s["read2"])

    if args.dry_run:
        master.warning("DRY-RUN mode — commands will NOT be executed")

    # Validate genome
    if not Path(args.genome).exists() and not args.dry_run:
        master.error("Genome not found: %s", args.genome)
        sys.exit(1)

    # ── Run all samples ───────────────────────────────────────────────
    total    = len(samples)
    passed   = []
    failed   = []
    timings  = {}          # sample_id → timedelta
    t_global = datetime.now()

    for i, sample in enumerate(samples, 1):
        master.info("")
        master.info("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        master.info("Processing sample %d / %d", i, total)

        sa              = sample_args(args, sample)
        success, elapsed = run_sample(sa, sample, master)
        timings[sample["sample_id"]] = elapsed

        if success:
            passed.append(sample["sample_id"])
        else:
            failed.append(sample["sample_id"])
            if args.stop_on_error:
                master.error("--stop-on-error set; aborting after first failure")
                break

    # ── Timing summary ────────────────────────────────────────────────
    total_elapsed = datetime.now() - t_global

    def fmt(td) -> str:
        """Format a timedelta as  Xh Ym Zs  (omitting leading zero units)."""
        secs  = int(td.total_seconds())
        h, r  = divmod(secs, 3600)
        m, s  = divmod(r, 60)
        if h:
            return f"{h}h {m}m {s}s"
        if m:
            return f"{m}m {s}s"
        return f"{s}s"

    SEP  = "╔" + "═" * 58 + "╗"
    MID  = "╠" + "═" * 58 + "╣"
    END  = "╚" + "═" * 58 + "╝"
    ROW  = "║  {:<28}  {:>10}  {:>12}  ║"

    master.info("")
    master.info(SEP)
    master.info("║%s║", "  PIPELINE TIMING SUMMARY".center(58))
    master.info(MID)
    master.info(ROW.format("Sample", "Status", "Time taken"))
    master.info("║" + "─" * 58 + "║")
    for sid, td in timings.items():
        status = "✓ PASS" if sid in passed else "✗ FAIL"
        master.info(ROW.format(sid, status, fmt(td)))
    master.info(MID)
    master.info(ROW.format("TOTAL  (%d sample/s)" % total, "", fmt(total_elapsed)))
    master.info(END)
    master.info("")
    master.info("  Samples passed : %d / %d", len(passed), total)
    if failed:
        master.error("  Samples failed : %s", ", ".join(failed))
    master.info("  Master log     : %s", global_log_dir / "pipeline_master.log")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
