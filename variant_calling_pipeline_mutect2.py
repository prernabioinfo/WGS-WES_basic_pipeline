#!/usr/bin/env python3
"""
======================================================================
  Variant Calling Pipeline — WGS/WES  (Mutect2 Tumor-Only + PoN)
  Automates: QC → Trimming → Alignment → PoN Creation → Mutect2
  Supports:
    • Normal-only  → builds Panel of Normals (PoN)
    • Tumor-only   → calls somatic variants using PoN
    • Full         → PoN build + tumor calling in one run
======================================================================

Author  : Prerna
Version : 4.0.0

SAMPLE SHEET FORMAT — TUMOR/NORMAL (TSV or CSV)
------------------------------------------------
No header row required. Each line = one sample:

  Column 1 : sample type  → "normal"  or  "tumor"
  Column 2 : path to forward reads (R1)
  Column 3 : path to reverse reads (R2)

Example (samples.tsv):
  normal  /data/normal1_R1.fastq.gz  /data/normal1_R2.fastq.gz
  normal  /data/normal2_R1.fastq.gz  /data/normal2_R2.fastq.gz
  tumor   /data/tumor1_R1.fastq.gz   /data/tumor1_R2.fastq.gz

Lines starting with '#' are treated as comments and ignored.
Delimiter is auto-detected (tab or comma or whitespace).

PIPELINE MODES
--------------
  --mode pon-only    : Align normals → Mutect2 (normal-only) per sample
                       → GenomicsDBImport → CreateSomaticPanelOfNormals
  --mode tumor-only  : Align tumors  → Mutect2 (tumor-only) → FilterMutectCalls
  --mode full        : Both phases in sequence (PoN first, then tumor calling)
"""

import argparse
import logging
import os
import shutil
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
    "trimmomatic_jar": "/home/qsbrp2022/prerna/apps/Trimmomatic-0.39/trimmomatic-0.39.jar",
    "adapters":        "/home/qsbrp2022/prerna/apps/Trimmomatic-0.39/adapters/TruSeq2-PE.fa",
    "picard_jar":      "/home/qsbrp2022/prerna/apps/picard.jar",
    "bwa":             "/home/qsbrp2022/prerna/apps/bwa/bwa",
    "gatk":            "/home/qsbrp2022/prerna/apps/gatk/gatk",
    "rgid":            "SAMPLE",
    "rglb":            "lib1",
    "rgpl":            "ILLUMINA",
    "rgpu":            "unit1",
    "rgsm":            "SAMPLE",
    "sliding_window":  "4:20",
    "min_len":         36,
    "phred":           "33",
    # Mutect2 / PoN
    "pon_db":          "pon_db",
    "pon_vcf":         "panel_of_normals.vcf.gz",
    "mode":            "full",
}

BANNER = r"""
╔══════════════════════════════════════════════════════════════════╗
║     W G S  V A R I A N T  C A L L I N G  P I P E L I N E      ║
║     BWA-MEM  ▸  Picard  ▸  Mutect2  (Tumor-Only + PoN)         ║
╚══════════════════════════════════════════════════════════════════╝
"""


# ─────────────────────────────────────────────────────────────────────
# SAMPLE SHEET PARSER
# ─────────────────────────────────────────────────────────────────────
def derive_sample_id(r1_path: str) -> str:
    name = Path(r1_path).name
    for suffix in ("_1.fastq.gz", "_1.fastq", "_R1.fastq.gz", "_R1.fastq",
                   "_1.fq.gz", "_1.fq", "_R1.fq.gz", "_R1.fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name.split(".")[0]).name


def parse_sample_sheet(path: str) -> tuple[list[dict], list[dict]]:
    """
    Parse a 3-column sample sheet: type | R1 | R2
    Returns (normals, tumors).
    """
    p = Path(path)
    if not p.exists():
        sys.exit(f"[ERROR] Sample sheet not found: {path}")

    lines = [ln for ln in p.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")]
    if not lines:
        sys.exit(f"[ERROR] Sample sheet is empty: {path}")

    first = lines[0]
    if "\t" in first:
        delim = "\t"
    elif "," in first:
        delim = ","
    else:
        delim = None

    normals: list[dict] = []
    tumors:  list[dict] = []

    for i, line in enumerate(lines, 1):
        cols = line.split(delim) if delim else line.split()
        cols = [c.strip() for c in cols if c.strip()]
        if len(cols) < 3:
            sys.exit(
                f"[ERROR] Line {i} has fewer than 3 columns:\n"
                f"        {line!r}\n"
                f"        Expected: <normal|tumor>  <R1_path>  <R2_path>"
            )
        stype, r1, r2 = cols[0].lower(), cols[1], cols[2]
        if stype not in ("normal", "tumor"):
            sys.exit(
                f"[ERROR] Line {i}: type must be 'normal' or 'tumor', got: {cols[0]!r}"
            )
        sid   = derive_sample_id(r1)
        entry = {"sample_id": sid, "read1": r1, "read2": r2}
        (normals if stype == "normal" else tumors).append(entry)

    for gname, group in [("normal", normals), ("tumor", tumors)]:
        seen: dict[str, int] = {}
        for s in group:
            seen[s["sample_id"]] = seen.get(s["sample_id"], 0) + 1
        dupes = [sid for sid, cnt in seen.items() if cnt > 1]
        if dupes:
            print(f"[WARNING] Duplicate {gname} IDs: {dupes}", file=sys.stderr)

    return normals, tumors


# ─────────────────────────────────────────────────────────────────────
# ARGUMENT PARSER
# ─────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant_pipeline_mutect2.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Full run — build PoN from normals then call tumors:
    python variant_pipeline_mutect2.py \\
        --sample-sheet samples.tsv -g ref.fna -o /results/ \\
        --intervals intervals.interval_list

  # Build PoN only:
    python variant_pipeline_mutect2.py \\
        --sample-sheet normals.tsv -g ref.fna -o /results/ \\
        --intervals intervals.interval_list --mode pon-only

  # Tumor calling with an existing PoN:
    python variant_pipeline_mutect2.py \\
        --sample-sheet tumors.tsv -g ref.fna -o /results/ \\
        --mode tumor-only --pon-vcf /results/panel_of_normals.vcf.gz

  # Dry run (log commands only):
    python variant_pipeline_mutect2.py \\
        --sample-sheet samples.tsv -g ref.fna -o /results/ \\
        --intervals intervals.interval_list --dry-run
""",
    )

    inp = parser.add_argument_group("Input")
    inp.add_argument("--sample-sheet", required=True, metavar="FILE",
                     help="Headerless file: col1=type(normal|tumor), col2=R1, col3=R2")

    req = parser.add_argument_group("Required")
    req.add_argument("-g", "--genome", required=True, metavar="FASTA",
                     help="Reference genome FASTA")
    req.add_argument("-o", "--outdir", required=True, metavar="DIR",
                     help="Root output directory")

    mut = parser.add_argument_group("Mutect2 / Panel of Normals")
    mut.add_argument("--mode", choices=["full", "pon-only", "tumor-only"],
                     default=DEFAULTS["mode"],
                     help=(
                         "full = PoN + tumor calling | "
                         "pon-only = build PoN only | "
                         "tumor-only = call tumors with existing --pon-vcf  "
                         f"[default: {DEFAULTS['mode']}]"
                     ))
    mut.add_argument("--intervals", metavar="FILE",
                     help="Interval list for GenomicsDBImport / Mutect2 (required for PoN)")
    mut.add_argument("--pon-db", default=DEFAULTS["pon_db"], metavar="DIR",
                     help=f"GenomicsDB workspace dir name (under --outdir)  [default: {DEFAULTS['pon_db']}]")
    mut.add_argument("--pon-vcf", default=None, metavar="VCF",
                     help=(
                         "PoN VCF output path (pon-only/full) or "
                         "path to existing PoN VCF (tumor-only)  "
                         "[default: <outdir>/panel_of_normals.vcf.gz]"
                     ))
    mut.add_argument("--germline-resource", default=None, metavar="VCF",
                     help="Germline resource VCF (e.g. gnomAD) for Mutect2 (optional)")
    mut.add_argument("--af-only-gnomad", action="store_true",
                     help="Pass --af-only-gnomad to Mutect2 (use with gnomAD AF-only VCF)")

    tools = parser.add_argument_group("Tool paths")
    tools.add_argument("--bwa",         default=DEFAULTS["bwa"])
    tools.add_argument("--trimmomatic", default=DEFAULTS["trimmomatic_jar"], metavar="JAR")
    tools.add_argument("--adapters",    default=DEFAULTS["adapters"])
    tools.add_argument("--picard",      default=DEFAULTS["picard_jar"], metavar="JAR")
    tools.add_argument("--gatk",        default=DEFAULTS["gatk"])

    perf = parser.add_argument_group("Performance")
    perf.add_argument("-t", "--threads", type=int, default=DEFAULTS["threads"],
                      help=f"Threads  [default: {DEFAULTS['threads']}]")
    perf.add_argument("--java-mem", default=DEFAULTS["java_mem"],
                      help=f"Java heap  [default: {DEFAULTS['java_mem']}]")

    rg = parser.add_argument_group("Default read group metadata")
    rg.add_argument("--rgid", default=DEFAULTS["rgid"])
    rg.add_argument("--rglb", default=DEFAULTS["rglb"])
    rg.add_argument("--rgpl", default=DEFAULTS["rgpl"])
    rg.add_argument("--rgpu", default=DEFAULTS["rgpu"])
    rg.add_argument("--rgsm", default=DEFAULTS["rgsm"])

    trim = parser.add_argument_group("Trimmomatic parameters")
    trim.add_argument("--sliding-window", default=DEFAULTS["sliding_window"])
    trim.add_argument("--min-len", type=int, default=DEFAULTS["min_len"])
    trim.add_argument("--phred", default=DEFAULTS["phred"], choices=["33", "64"])

    ctrl = parser.add_argument_group("Pipeline control")
    ctrl.add_argument("--skip-qc",            action="store_true")
    ctrl.add_argument("--skip-trim",          action="store_true")
    ctrl.add_argument("--skip-alignment",     action="store_true")
    ctrl.add_argument("--dry-run",            action="store_true")
    ctrl.add_argument("--keep-intermediates", action="store_true")
    ctrl.add_argument("--stop-on-error",      action="store_true")

    misc = parser.add_argument_group("Miscellaneous")
    misc.add_argument("--log-level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    misc.add_argument("-v", "--version", action="version", version="%(prog)s 4.0.0")

    return parser


# ─────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────
def setup_master_logging(log_dir: Path, level: str) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    fmt  = "%(asctime)s  [%(levelname)s]  %(message)s"
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(getattr(logging, level))
    sh.setFormatter(logging.Formatter(fmt))
    root.addHandler(sh)
    fh = logging.FileHandler(log_dir / "pipeline_master.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(fmt))
    root.addHandler(fh)
    return logging.getLogger("pipeline")


def step_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(log_path.stem)
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
        raise RuntimeError(f"Step failed — see log: {_log_path(logger)}")


def _log_path(logger: logging.Logger) -> str:
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            return h.baseFilename
    return "<unknown>"


# ─────────────────────────────────────────────────────────────────────
# PER-SAMPLE ARGS
# ─────────────────────────────────────────────────────────────────────
def sample_args(global_args, sample: dict):
    sa      = deepcopy(global_args)
    sid     = sample["sample_id"]
    sa.rgid = sample.get("rgid") or sid
    sa.rglb = sample.get("rglb") or global_args.rglb
    sa.rgpl = sample.get("rgpl") or global_args.rgpl
    sa.rgpu = sample.get("rgpu") or global_args.rgpu
    sa.rgsm = sample.get("rgsm") or sid
    return sa


# ─────────────────────────────────────────────────────────────────────
# ALIGNMENT STEPS  01–07  (shared by normals and tumors)
# ─────────────────────────────────────────────────────────────────────

def step_qc_raw(args, dirs, r1, r2, master):
    master.info("  ── STEP 01 : Raw QC (FastQC + seqkit)")
    log = step_logger(dirs["logs"] / "step01_qc_raw.log")
    run(f"seqkit stats {r1} {r2} > {dirs['qc_raw']}/seqkit_stats.txt", log, args.dry_run)
    run(f"fastqc -t {args.threads} {r1} {r2} -o {dirs['qc_raw']}", log, args.dry_run)


def step_trim(args, dirs, base, r1, r2, master):
    master.info("  ── STEP 02 : Trimmomatic")
    log = step_logger(dirs["logs"] / "step02_trimmomatic.log")
    p1  = f"{dirs['trim']}/{base}_R1_paired.fastq.gz"
    u1  = f"{dirs['trim']}/{base}_R1_unpaired.fastq.gz"
    p2  = f"{dirs['trim']}/{base}_R2_paired.fastq.gz"
    u2  = f"{dirs['trim']}/{base}_R2_unpaired.fastq.gz"
    run(
        f"java -Xmx{args.java_mem} -jar {args.trimmomatic} PE "
        f"-threads {args.threads} -{args.phred} "
        f"{r1} {r2} {p1} {u1} {p2} {u2} "
        f"ILLUMINACLIP:{args.adapters}:2:30:10 "
        f"SLIDINGWINDOW:{args.sliding_window} MINLEN:{args.min_len}",
        log, args.dry_run,
    )
    return p1, p2


def step_qc_trimmed(args, dirs, p1, p2, master):
    master.info("  ── STEP 03 : QC on trimmed reads")
    log = step_logger(dirs["logs"] / "step03_qc_trimmed.log")
    run(f"seqkit stats {p1} {p2} > {dirs['qc_trim']}/seqkit_stats.txt", log, args.dry_run)
    run(f"fastqc -t {args.threads} {p1} {p2} -o {dirs['qc_trim']}", log, args.dry_run)


def step_align(args, dirs, base, r1, r2, master):
    master.info("  ── STEP 04 : BWA-MEM alignment")
    log = step_logger(dirs["logs"] / "step04_bwa_mem.log")
    sam = f"{dirs['aln']}/{base}.sam"
    run(f"{args.bwa} mem -t {args.threads} {args.genome} {r1} {r2} -o {sam}",
        log, args.dry_run)
    return sam


def step_sort(args, dirs, base, sam, master):
    master.info("  ── STEP 05 : Picard SortSam")
    log        = step_logger(dirs["logs"] / "step05_picard_sort.log")
    sorted_bam = f"{dirs['aln']}/{base}.sorted.bam"
    tmp        = f"{dirs['aln']}/tmp"
    os.makedirs(tmp, exist_ok=True)
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} SortSam "
        f"I={sam} SORT_ORDER=coordinate O={sorted_bam} TMP_DIR={tmp}",
        log, args.dry_run,
    )
    if not args.keep_intermediates and not args.dry_run:
        Path(sam).unlink(missing_ok=True)
    return sorted_bam


def step_add_rg(args, dirs, base, sorted_bam, master):
    master.info("  ── STEP 06 : Picard AddOrReplaceReadGroups")
    log    = step_logger(dirs["logs"] / "step06_picard_rg.log")
    rg_bam = f"{dirs['aln']}/{base}.rg.bam"
    tmp    = f"{dirs['aln']}/tmp"
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} AddOrReplaceReadGroups "
        f"I={sorted_bam} O={rg_bam} "
        f"RGID={args.rgid} RGLB={args.rglb} RGPL={args.rgpl} "
        f"RGPU={args.rgpu} RGSM={args.rgsm} TMP_DIR={tmp}",
        log, args.dry_run,
    )
    if not args.keep_intermediates and not args.dry_run:
        Path(sorted_bam).unlink(missing_ok=True)
    return rg_bam


def step_mark_dups(args, dirs, base, rg_bam, master):
    master.info("  ── STEP 07 : Picard MarkDuplicates + samtools index")
    log       = step_logger(dirs["logs"] / "step07_picard_markdup.log")
    dedup_bam = f"{dirs['aln']}/{base}.bam"
    metrics   = f"{dirs['aln']}/{base}.dup_metrics.txt"
    tmp       = f"{dirs['aln']}/tmp"
    run(
        f"java -Xmx{args.java_mem} -jar {args.picard} MarkDuplicates "
        f"I={rg_bam} O={dedup_bam} M={metrics} "
        f"REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR={tmp}",
        log, args.dry_run,
    )
    run(f"samtools index {dedup_bam}", log, args.dry_run)
    if not args.keep_intermediates and not args.dry_run:
        Path(rg_bam).unlink(missing_ok=True)
    return dedup_bam


# ─────────────────────────────────────────────────────────────────────
# PON STEPS  08–10
# ─────────────────────────────────────────────────────────────────────

def step_mutect2_normal(args, dirs, base, bam, master) -> str:
    """
    STEP 08 — Mutect2 in normal-only mode.

    gatk Mutect2 -R reference.fasta -I normal.bam -max-mnp-distance 0 -O normal.vcf.gz
    """
    master.info("  ── STEP 08 : Mutect2 normal-only call (for PoN)")
    log = step_logger(dirs["logs"] / "step08_mutect2_normal.log")
    vcf = f"{dirs['vc']}/{base}_normal.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" Mutect2 "
        f"-R {args.genome} "
        f"-I {bam} "
        f"-max-mnp-distance 0 "
        f"-O {vcf}",
        log, args.dry_run,
    )
    master.info("     output → %s", vcf)
    master.info("     log    → %s", dirs['logs'] / 'step08_mutect2_normal.log')
    return vcf


def step_genomicsdb_import(args, normal_vcfs: list, pon_db_path: str,
                            global_log_dir: Path, master) -> None:
    """
    STEP 09 — GenomicsDBImport: consolidate all normal VCFs into a workspace.

    gatk GenomicsDBImport -R reference.fasta \\
        --genomicsdb-workspace-path pon_db \\
        -V normal1.vcf.gz -V normal2.vcf.gz ...
    """
    master.info("  ── STEP 09 : GenomicsDBImport — consolidate normal VCFs")
    log = step_logger(global_log_dir / "step09_genomicsdb_import.log")

    # GATK refuses to overwrite an existing workspace — remove it first
    if Path(pon_db_path).exists() and not args.dry_run:
        master.warning("  Removing existing GenomicsDB workspace: %s", pon_db_path)
        shutil.rmtree(pon_db_path)

    v_flags = " ".join(f"-V {v}" for v in normal_vcfs)

    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" GenomicsDBImport "
        f"-R {args.genome} "
        f"--genomicsdb-workspace-path {pon_db_path} "
        f"{v_flags}",
        log, args.dry_run,
    )
    master.info("     GenomicsDB workspace → %s", pon_db_path)
    master.info("     log → %s", global_log_dir / 'step09_genomicsdb_import.log')


def step_create_somatic_pon(args, pon_db_path: str, pon_vcf_path: str,
                             global_log_dir: Path, master) -> str:
    """
    STEP 10 — CreateSomaticPanelOfNormals.

    gatk CreateSomaticPanelOfNormals -R reference.fasta \\
        -V gendb://pon_db -O panel_of_normals.vcf.gz
    """
    master.info("  ── STEP 10 : CreateSomaticPanelOfNormals")
    log = step_logger(global_log_dir / "step10_create_somatic_pon.log")
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" CreateSomaticPanelOfNormals "
        f"-R {args.genome} "
        f"-V gendb://{pon_db_path} "
        f"-O {pon_vcf_path}",
        log, args.dry_run,
    )
    master.info("     PoN VCF → %s", pon_vcf_path)
    master.info("     log     → %s", global_log_dir / 'step10_create_somatic_pon.log')
    return pon_vcf_path


# ─────────────────────────────────────────────────────────────────────
# TUMOR CALLING STEPS  11–12
# ─────────────────────────────────────────────────────────────────────

def step_mutect2_tumor(args, dirs, base, bam, pon_vcf, master) -> str:
    """
    STEP 11 — Mutect2 tumor-only somatic calling with PoN.

    gatk Mutect2 -R reference.fasta -I tumor.bam \\
        --panel-of-normals panel_of_normals.vcf.gz \\
        [--germline-resource gnomad.vcf.gz [--af-only-gnomad]] \\
        [-L intervals.interval_list] \\
        -O tumor_unfiltered.vcf.gz
    """
    master.info("  ── STEP 11 : Mutect2 tumor-only somatic calling")
    log     = step_logger(dirs["logs"] / "step11_mutect2_tumor.log")
    raw_vcf = f"{dirs['vc']}/{base}_somatic_unfiltered.vcf.gz"

    pon_flag      = f"--panel-of-normals {pon_vcf}" if pon_vcf else ""
    germline_flag = (f"--germline-resource {args.germline_resource}"
                     if args.germline_resource else "")
    af_only_flag  = "--af-only-gnomad" if args.af_only_gnomad else ""
    interval_flag = f"-L {args.intervals}" if args.intervals else ""

    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" Mutect2 "
        f"-R {args.genome} "
        f"-I {bam} "
        f"{pon_flag} "
        f"{germline_flag} "
        f"{af_only_flag} "
        f"{interval_flag} "
        f"-O {raw_vcf}",
        log, args.dry_run,
    )
    master.info("     raw somatic VCF → %s", raw_vcf)
    master.info("     log             → %s", dirs['logs'] / 'step11_mutect2_tumor.log')
    return raw_vcf


def step_filter_mutect_calls(args, dirs, base, raw_vcf, master) -> str:
    """
    STEP 12 — FilterMutectCalls: produce the final filtered somatic VCF.

    gatk FilterMutectCalls -R reference.fasta \\
        -V tumor_unfiltered.vcf.gz -O tumor_filtered.vcf.gz
    """
    master.info("  ── STEP 12 : FilterMutectCalls")
    log      = step_logger(dirs["logs"] / "step12_filter_mutect_calls.log")
    filt_vcf = f"{dirs['vc']}/{base}_somatic_filtered.vcf.gz"
    run(
        f"{args.gatk} --java-options \"-Xmx{args.java_mem}\" FilterMutectCalls "
        f"-R {args.genome} "
        f"-V {raw_vcf} "
        f"-O {filt_vcf}",
        log, args.dry_run,
    )
    master.info("     filtered somatic VCF → %s", filt_vcf)
    master.info("     log                  → %s", dirs['logs'] / 'step12_filter_mutect_calls.log')
    return filt_vcf


# ─────────────────────────────────────────────────────────────────────
# ALIGNMENT PHASE RUNNER  (shared by normal + tumor samples)
# ─────────────────────────────────────────────────────────────────────

def run_alignment_phase(args, sample: dict, master):
    """Steps 01–07 for a single sample. Returns (final_bam, dirs)."""
    base = sample["sample_id"]
    r1   = sample["read1"]
    r2   = sample["read2"]

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

    master.info("  Output    : %s", root_out)
    master.info("  Step logs : %s/", dirs['logs'])

    for label, fpath in [("R1", r1), ("R2", r2)]:
        if not Path(fpath).exists() and not args.dry_run:
            raise RuntimeError(f"{label} not found: {fpath}")

    if not args.skip_qc:
        step_qc_raw(args, dirs, r1, r2, master)
    else:
        master.info("  Skipping raw QC")

    if not args.skip_trim:
        p1, p2 = step_trim(args, dirs, base, r1, r2, master)
        if not args.skip_qc:
            step_qc_trimmed(args, dirs, p1, p2, master)
    else:
        master.info("  Skipping trimming")
        p1, p2 = r1, r2

    if not args.skip_alignment:
        sam       = step_align(args, dirs, base, p1, p2, master)
        srt_bam   = step_sort(args, dirs, base, sam, master)
        rg_bam    = step_add_rg(args, dirs, base, srt_bam, master)
        final_bam = step_mark_dups(args, dirs, base, rg_bam, master)
    else:
        final_bam = str(dirs['aln'] / f"{base}.bam")
        master.info("  Skipping alignment; expecting BAM: %s", final_bam)

    return final_bam, dirs


# ─────────────────────────────────────────────────────────────────────
# PHASE ORCHESTRATORS
# ─────────────────────────────────────────────────────────────────────

def run_pon_phase(args, normals: list, pon_vcf_path: str,
                  global_log_dir: Path, master) -> str:
    """PHASE 1 — align normals, call each with Mutect2, build PoN."""
    if not normals:
        master.error("No normal samples — cannot build PoN.")
        sys.exit(1)

    master.info("")
    master.info("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    master.info("  PHASE 1 : Panel of Normals (PoN) Construction")
    master.info("  Normal samples : %d", len(normals))
    master.info("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

    normal_vcfs: list = []

    for i, sample in enumerate(normals, 1):
        base = sample["sample_id"]
        master.info("")
        master.info("╔══ NORMAL %d/%d : %s ══", i, len(normals), base)
        t0 = datetime.now()
        try:
            sa              = sample_args(args, sample)
            final_bam, dirs = run_alignment_phase(sa, sample, master)
            normal_vcf      = step_mutect2_normal(sa, dirs, base, final_bam, master)
            normal_vcfs.append(normal_vcf)
        except RuntimeError as exc:
            master.error("  FAILED : %s  (%s)", base, exc)
            if args.stop_on_error:
                sys.exit(1)
            continue
        master.info("╚══ DONE  : %s  (%s)", base,
                    str(datetime.now() - t0).split(".")[0])

    if not normal_vcfs:
        master.error("All normal samples failed — cannot build PoN.")
        sys.exit(1)

    master.info("")
    master.info("  %d normal VCF(s) collected:", len(normal_vcfs))
    for v in normal_vcfs:
        master.info("    %s", v)

    pon_db_path = str(Path(args.outdir) / args.pon_db)
    step_genomicsdb_import(args, normal_vcfs, pon_db_path, global_log_dir, master)
    step_create_somatic_pon(args, pon_db_path, pon_vcf_path, global_log_dir, master)

    master.info("")
    master.info("  ✓  PoN ready : %s", pon_vcf_path)
    return pon_vcf_path


def run_tumor_phase(args, tumors: list, pon_vcf: str, master):
    """PHASE 2 — align tumors, run Mutect2 tumor-only, filter calls."""
    if not tumors:
        master.warning("No tumor samples — skipping tumor calling phase.")
        return [], []

    master.info("")
    master.info("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    master.info("  PHASE 2 : Tumor-Only Somatic Variant Calling")
    master.info("  Tumor samples : %d", len(tumors))
    master.info("  PoN VCF       : %s", pon_vcf or "None (not recommended)")
    master.info("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

    passed: list = []
    failed: list = []

    for i, sample in enumerate(tumors, 1):
        base = sample["sample_id"]
        master.info("")
        master.info("╔══ TUMOR %d/%d : %s ══", i, len(tumors), base)
        t0 = datetime.now()
        try:
            sa              = sample_args(args, sample)
            final_bam, dirs = run_alignment_phase(sa, sample, master)
            raw_vcf         = step_mutect2_tumor(sa, dirs, base, final_bam, pon_vcf, master)
            filt_vcf        = step_filter_mutect_calls(sa, dirs, base, raw_vcf, master)
            master.info("  Final somatic VCF : %s", filt_vcf)
            passed.append(base)
        except RuntimeError as exc:
            master.error("  FAILED : %s  (%s)", base, exc)
            failed.append(base)
            if args.stop_on_error:
                break
        master.info("╚══ DONE  : %s  (%s)", base,
                    str(datetime.now() - t0).split(".")[0])

    return passed, failed


# ─────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────
def main() -> None:
    print(BANNER)
    parser = build_parser()
    args   = parser.parse_args()

    # Resolve PoN VCF path default
    if args.pon_vcf is None:
        args.pon_vcf = str(Path(args.outdir) / DEFAULTS["pon_vcf"])

    # Parse sample sheet
    normals, tumors = parse_sample_sheet(args.sample_sheet)

    # Master logger
    global_log_dir = Path(args.outdir) / "_pipeline_logs"
    master = setup_master_logging(global_log_dir, args.log_level)

    master.info("Pipeline v4.0.0  (mode: %s)", args.mode)
    master.info("Reference genome  : %s", args.genome)
    master.info("Output root       : %s", args.outdir)
    master.info("Threads           : %d", args.threads)
    master.info("Java heap         : %s", args.java_mem)
    master.info("Intervals         : %s", args.intervals or "not set")
    master.info("PoN VCF path      : %s", args.pon_vcf)
    master.info("Normal samples    : %d", len(normals))
    master.info("Tumor  samples    : %d", len(tumors))
    if args.dry_run:
        master.warning("DRY-RUN — commands will NOT be executed")

    if not Path(args.genome).exists() and not args.dry_run:
        master.error("Genome not found: %s", args.genome)
        sys.exit(1)

    # Mode validation
    if args.mode in ("full", "pon-only") and not normals:
        master.error("Mode '%s' requires normal samples but none were found.", args.mode)
        sys.exit(1)

    if args.mode == "tumor-only":
        pon_ok = args.pon_vcf and (args.dry_run or Path(args.pon_vcf).exists())
        if not pon_ok:
            master.warning(
                "Mode 'tumor-only' but PoN VCF not found at '%s'. "
                "Running without PoN (strongly not recommended).",
                args.pon_vcf,
            )

    t_global = datetime.now()
    pon_vcf  = args.pon_vcf
    t_passed: list = []
    t_failed: list = []

    if args.mode in ("full", "pon-only"):
        pon_vcf = run_pon_phase(args, normals, args.pon_vcf, global_log_dir, master)

    if args.mode in ("full", "tumor-only"):
        t_passed, t_failed = run_tumor_phase(args, tumors, pon_vcf, master)

    # ── Final summary ─────────────────────────────────────────────────
    elapsed = datetime.now() - t_global
    SEP = "╔" + "═" * 60 + "╗"
    MID = "╠" + "═" * 60 + "╣"
    END = "╚" + "═" * 60 + "╝"

    master.info("")
    master.info(SEP)
    master.info("║%s║", "  PIPELINE SUMMARY".center(60))
    master.info(MID)
    master.info("║  Mode           : %-42s ║", args.mode)
    if args.mode in ("full", "pon-only"):
        master.info("║  Normals used   : %-42s ║", len(normals))
        master.info("║  PoN VCF        : %-42s ║", str(pon_vcf)[-42:])
    if args.mode in ("full", "tumor-only"):
        master.info("║  Tumors passed  : %-42s ║",
                    f"{len(t_passed)} / {len(tumors)}")
        if t_failed:
            master.info("║  Tumors FAILED  : %-42s ║", ", ".join(t_failed))
    master.info(MID)
    master.info("║  Total elapsed  : %-42s ║", str(elapsed).split(".")[0])
    master.info(END)
    master.info("")
    master.info("  Master log : %s", global_log_dir / "pipeline_master.log")

    if t_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
