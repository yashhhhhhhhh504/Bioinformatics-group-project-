import argparse
from pathlib import Path
import re
from cyvcf2 import VCF
import pandas as pd
import numpy as np
COARSE_FUNC_MAP = {
    "DEL": "sequence_loss (coarse)",
    "INS": "sequence_gain (coarse)",
    "DUP": "copy_gain (coarse)",
    "INV": "inversion (coarse)",
}
SEVERITY_RANK = {
    "transcript_ablation": 1,
    "exon_loss_variant": 2,
    "frameshift_variant": 3,
    "stop_gained": 4,
    "start_lost": 5,
    "splice_acceptor_variant": 6,
    "splice_donor_variant": 7,
    "splice_region_variant": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "missense_variant": 11,
    "protein_altering_variant": 12,
    "coding_sequence_variant": 13,
    "incomplete_terminal_codon_variant": 14,
    "synonymous_variant": 15,
    "5_prime_UTR_variant": 20,
    "3_prime_UTR_variant": 21,
    "upstream_gene_variant": 30,
    "downstream_gene_variant": 31,
    "intron_variant": 40,
    "regulatory_region_variant": 45,
    "TF_binding_site_variant": 46,
    "intergenic_variant": 50,
}
def parse_args():
    ap = argparse.ArgumentParser(
        description="Summarize trio SV combinations (ALT/REF/NA), locations, type and function directly from VCF/BCF."
    )
    ap.add_argument("--vcf", required=True, help="Input VCF/BCF (.vcf/.vcf.gz/.bcf)")
    ap.add_argument("--samples", nargs=3, required=True, metavar=("S1","S2","S3"),
                    help="Three sample IDs (e.g., HG00512 HG00513 HG00514)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (path or basename)")
    ap.add_argument("--max-examples", type=int, default=5, help="Number of example loci per row (default 5)")
    return ap.parse_args()

# ---------- Helpers for genotypes & lengths ----------

def gt_state_from_cyvcf2(gt_tuple):
    """
    cyvcf2 genotype is [a, b, phased_flag, ...]; -1 means missing.
    Return ALT/REF/NA (ALT if any non-zero allele; REF if 0/0; NA otherwise).
    """
    if gt_tuple is None or len(gt_tuple) < 2:
        return "NA"
    a, b = gt_tuple[0], gt_tuple[1]
    if a == -1 or b == -1:
        return "NA"
    if a == 0 and b == 0:
        return "REF"
    # any non-zero allele -> ALT
    if a > 0 or b > 0:
        return "ALT"
    return "NA"

def get_info_int(info, key, fallback=None):
    v = info.get(key, None)
    if v is None:
        return fallback
    try:
        return int(v)
    except Exception:
        try:
            return int(float(v))
        except Exception:
            return fallback

def estimate_svlen(svtype, pos, end, inslen):
    # Prefer INSLEN for INS; else end - pos; else NaN
    if svtype == "INS" and inslen is not None:
        return abs(inslen)
    if pos is not None and end is not None:
        return max(0, abs(end - pos))
    return np.nan

# ---------- VEP/SnpEff parsing for function ----------

def parse_csq_header(vcf):
    """
    Return tuple (has_csq, csq_fields) where csq_fields is the ordered list from the CSQ header.
    """
    raw = vcf.raw_header or ""
    m = re.search(r'##INFO=<ID=CSQ,.*?Format: ([^">]+)>', raw)
    if not m:
        return False, []
    fields = m.group(1).strip().split("|")
    return True, fields

def parse_ann_header(vcf):
    """
    Return tuple (has_ann, ann_fields) where ann_fields is the ordered list from the ANN header (SnpEff).
    """
    raw = vcf.raw_header or ""
    m = re.search(r'##INFO=<ID=ANN,.*?Functional annotations: \'?([^\'">]+)\'?>', raw)
    if not m:
        # Some SnpEff headers use "Annotation_Fields: " instead
        m = re.search(r'##INFO=<ID=ANN,.*?Annotation_Fields: ([^">]+)>', raw)
    if not m:
        return False, []
    fields = [f.strip() for f in m.group(1).split("|")]
    return True, fields

def most_severe_consequence_from_csq(csq_value, csq_fields):
    """
    Parse CSQ string(s) and return the most severe SO consequence term (by SEVERITY_RANK).
    csq_value may be a comma-separated string of transcripts; each transcript is pipe-delimited by csq_fields order.
    """
    if not csq_value:
        return None
    best_term, best_rank = None, 10**9
    transcripts = str(csq_value).split(",")
    try:
        idx_consequence = csq_fields.index("Consequence")
    except ValueError:
        return None
    for tr in transcripts:
        parts = tr.split("|")
        if idx_consequence >= len(parts):
            continue
        # Consequence can be "term1&term2"
        for term in parts[idx_consequence].split("&"):
            term = term.strip()
            rank = SEVERITY_RANK.get(term, 9999)
            if rank < best_rank:
                best_term, best_rank = term, rank
    return best_term

def most_severe_consequence_from_ann(ann_value, ann_fields):
    """
    Parse ANN string(s) and return the most severe SO consequence term (by SEVERITY_RANK).
    ANN is comma-separated; each part is pipe-delimited by ann_fields order.
    Commonly, the "Annotation" field carries the consequence term(s).
    """
    if not ann_value:
        return None
    best_term, best_rank = None, 10**9
    try:
        idx_annotation = ann_fields.index("Annotation")
    except ValueError:
        # Some SnpEff headers use "Annotation" or "Putative_Impact", but we need consequence term(s).
        return None
    transcripts = str(ann_value).split(",")
    for tr in transcripts:
        parts = tr.split("|")
        if idx_annotation >= len(parts):
            continue
        for term in parts[idx_annotation].split("&"):
            term = term.strip()
            rank = SEVERITY_RANK.get(term, 9999)
            if rank < best_rank:
                best_term, best_rank = term, rank
    return best_term

# ---------- Summarization helpers ----------

def top_chroms(series, k=3):
    vc = series.value_counts()
    items = [f"{chrom}({int(cnt)})" for chrom, cnt in vc.head(k).items()]
    return ",".join(items)

def example_loci(df_sub, n=5):
    rows = []
    for _, r in df_sub.head(n).iterrows():
        rows.append(f"{r['CHROM']}:{r['POS']}-{r['END']}")
    return "; ".join(rows)

def main():
    args = parse_args()
    vcf = VCF(args.vcf)
    samples = list(vcf.samples)

    # Trio sample indices
    try:
        s1, s2, s3 = args.samples
        i1, i2, i3 = samples.index(s1), samples.index(s2), samples.index(s3)
    except ValueError as e:
        raise SystemExit(f"[Error] Sample not found in VCF header: {e}\nAvailable: {samples}")

    # Detect annotation headers (VEP CSQ / SnpEff ANN)
    has_csq, csq_fields = parse_csq_header(vcf)
    has_ann, ann_fields = parse_ann_header(vcf)

    rows = []
    for var in vcf:
        chrom = var.CHROM
        pos   = var.POS  # 1-based
        info  = var.INFO
        svtype = info.get("SVTYPE", None)
        end    = get_info_int(info, "END", getattr(var, "end", None))
        inslen = get_info_int(info, "INSLEN", None)

        # Genotype states
        gts = var.genotypes  # list of [a,b,phased,...]
        s1_state = gt_state_from_cyvcf2(gts[i1]) if i1 < len(gts) else "NA"
        s2_state = gt_state_from_cyvcf2(gts[i2]) if i2 < len(gts) else "NA"
        s3_state = gt_state_from_cyvcf2(gts[i3]) if i3 < len(gts) else "NA"

        svlen_est = estimate_svlen(svtype, pos, end, inslen)

        # Function term: prefer CSQ, then ANN, else coarse SV-type mapping
        func_term = None
        if has_csq and "CSQ" in var.INFO:
            func_term = most_severe_consequence_from_csq(var.INFO.get("CSQ"), csq_fields)
        if func_term is None and has_ann and "ANN" in var.INFO:
            func_term = most_severe_consequence_from_ann(var.INFO.get("ANN"), ann_fields)
        if func_term is None:
            func_term = COARSE_FUNC_MAP.get(str(svtype), "unknown")

        rows.append({
            "CHROM": chrom,
            "POS": pos,
            "END": end,
            "SVTYPE": svtype,
            "INSLEN": inslen if inslen is not None else np.nan,
            "SVLEN_EST": svlen_est,
            f"{s1}_STATE": s1_state,
            f"{s2}_STATE": s2_state,
            f"{s3}_STATE": s3_state,
            "FUNCTION": func_term,  # final function term
        })

    df = pd.DataFrame(rows)

    # Combination label: e.g. "HG00512=ALT|HG00513=REF|HG00514=REF"
    df["COMBINATION"] = (
        f"{s1}=" + df[f"{s1}_STATE"] + "|" +
        f"{s2}=" + df[f"{s2}_STATE"] + "|" +
        f"{s3}=" + df[f"{s3}_STATE"]
    )

    # --- Main summary: COMBINATION × SVTYPE × FUNCTION ---
    grp_keys = ["COMBINATION", "SVTYPE", "FUNCTION"]
    grp = df.groupby(grp_keys, dropna=False)

    def med_len(x):
        xx = pd.to_numeric(x, errors="coerce")
        return float(np.nanmedian(xx)) if xx.notna().any() else np.nan

    summary = grp.agg(
        COUNT=("SVTYPE", "size"),
        MEDIAN_LEN_BP=("SVLEN_EST", med_len),
        TOP_CHROMS=("CHROM", top_chroms),
    ).reset_index()

    # Example loci
    examples = []
    for keys, idx in grp.groups.items():
        sub = df.loc[idx, ["CHROM","POS","END"]]
        examples.append(example_loci(sub, n=args.max_examples))
    summary["EXAMPLE_LOCI"] = examples

    # Per-chromosome breakdown
    by_chrom = (
        df.groupby(["COMBINATION","CHROM","SVTYPE","FUNCTION"], dropna=False)
          .size().reset_index(name="COUNT")
    )

    out1 = f"{args.out_prefix}.combination_summary.csv"
    out2 = f"{args.out_prefix}.combination_by_chrom.csv"
    summary.to_csv(out1, index=False)
    by_chrom.to_csv(out2, index=False)

    # Minimal self-checks to guard correctness
    total_rows = len(df)
    summed = summary["COUNT"].sum()
    if summed != total_rows:
        print(f"[WARN] Grouped count ({summed}) != total variants ({total_rows}). "
              "This can happen if there are missing SVTYPE/FUNCTION in some rows.")
    print(f"Wrote:\n  {out1}\n  {out2}")
    print(f"Variants processed: {total_rows}; Summary rows: {len(summary)}")

if __name__ == "__main__":
    main()