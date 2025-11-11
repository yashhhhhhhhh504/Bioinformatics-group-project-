#!/usr/bin/env gawk -f
# Usage:
#   gawk -f trio_summary_v2.awk \
#        -v s1=HG00512 -v s2=HG00513 -v s3=HG00514 \
#        -v out1=DellyVariation.combination_summary.tsv \
#        -v out2=DellyVariation.combination_by_chrom.tsv \
#        -v maxExamples=5 \
#        DellyVariation.vcf
#
# Outputs:
#   out1: COMBINATION, SVTYPE, FUNCTION, COUNT, MEDIAN_LEN_BP, TOP_CHROMS, EXAMPLE_LOCI
#   out2: COMBINATION, CHROM, SVTYPE, FUNCTION, COUNT

BEGIN {
  if (s1=="" || s2=="" || s3=="") {
    print "ERROR: pass -v s1=... -v s2=... -v s3=..." > "/dev/stderr"; exit 1
  }
  if (out1=="") out1="combination_summary.tsv"
  if (out2=="") out2="combination_by_chrom.tsv"
  if (maxExamples=="") maxExamples=5
  FS = "\t"; OFS = "\t"
}

# Skip meta headers
/^##/ { next }

# Header: find columns
/^#CHROM/ {
  for (i=1; i<=NF; i++) hdr[$i] = i
  fmtCol = hdr["FORMAT"]
  s1col = hdr[s1]; s2col = hdr[s2]; s3col = hdr[s3]
  if (!fmtCol || !s1col || !s2col || !s3col) {
    print "ERROR: could not locate FORMAT and/or sample columns in header" > "/dev/stderr"; exit 1
  }
  next
}

# Data lines
!/^#/ {
  CHROM = $1; POS=$2; INFO=$8; FORMAT=$9

  ENDv    = infoget(INFO,"END")
  INSLENv = infoget(INFO,"INSLEN")
  SVTYPE  = infoget(INFO,"SVTYPE"); if (SVTYPE=="") SVTYPE="."

  # length estimate
  len = ""
  if (SVTYPE=="INS" && INSLENv!="")      len = abs(INSLENv+0)
  else if (ENDv!="") { len = (ENDv+0)-(POS+0); if (len<0) len=-len }

  # coarse function
  if (SVTYPE=="DEL")      FUNCTION="sequence_loss (coarse)"
  else if (SVTYPE=="INS") FUNCTION="sequence_gain (coarse)"
  else if (SVTYPE=="DUP") FUNCTION="copy_gain (coarse)"
  else if (SVTYPE=="INV") FUNCTION="inversion (coarse)"
  else                    FUNCTION="unknown"

  # FORMAT field positions
  nfmt = split(FORMAT, F, ":"); delete pos
  for (i=1;i<=nfmt;i++) pos[F[i]] = i

  # per-sample ALT/REF/NA
  s1state = sample_state($s1col, pos)
  s2state = sample_state($s2col, pos)
  s3state = sample_state($s3col, pos)

  combo = s1 "=" s1state "|" s2 "=" s2state "|" s3 "=" s3state

  # aggregate: COMBINATION × SVTYPE × FUNCTION
  key = combo SUBSEP SVTYPE SUBSEP FUNCTION
  count[key]++
  if (len != "") {
    lenCount[key]++
    lenList[key, lenCount[key]] = len+0
  }
  if (exCount[key] < maxExamples) {
    locus = CHROM ":" POS "-" (ENDv=="" ? POS : ENDv)
    exCount[key]++
    exList[key, exCount[key]] = locus
  }
  chrCounts[key SUBSEP CHROM]++

  # by-chrom table
  bychrom[combo SUBSEP CHROM SUBSEP SVTYPE SUBSEP FUNCTION]++
  next
}

END {
  # by-chrom output
  print "COMBINATION","CHROM","SVTYPE","FUNCTION","COUNT" > out2
  for (k in bychrom) {
    split(k, p, SUBSEP)
    print p[1], p[2], p[3], p[4], bychrom[k] >> out2
  }

  # combination summary output
  print "COMBINATION","SVTYPE","FUNCTION","COUNT","MEDIAN_LEN_BP","TOP_CHROMS","EXAMPLE_LOCI" > out1
  for (k in count) {
    split(k, p, SUBSEP)
    combo = p[1]
    svt   = p[2]
    funcv = p[3]

    # median
    n = lenCount[k] + 0
    med = "."
    if (n > 0) {
      delete tmp
      for (i=1; i<=n; i++) tmp[i] = lenList[k, i]
      asort(tmp)  # numeric asc (gawk)
      if (n % 2 == 1) med = tmp[(n+1)/2]
      else            med = (tmp[n/2] + tmp[n/2 + 1]) / 2.0
    }

    # top 3 chromosomes
    top1c=0; top2c=0; top3c=0; top1=""; top2=""; top3=""
    for (cc in chrCounts) {
      split(cc, q, SUBSEP)
      if (q[1] != k) continue
      cchr = q[2]; val = chrCounts[cc] + 0
      if (val > top1c) { top3c=top2c; top3=top2; top2c=top1c; top2=top1; top1c=val; top1=cchr }
      else if (val > top2c) { top3c=top2c; top3=top2; top2c=val; top2=cchr }
      else if (val > top3c) { top3c=val; top3=cchr }
    }
    top = "."
    if (top1c>0) {
      top = top1 "(" top1c ")"
      if (top2c>0) top = top "," top2 "(" top2c ")"
      if (top3c>0) top = top "," top3 "(" top3c ")"
    }

    # example loci
    exs = "."
    if (exCount[k] > 0) {
      exs = ""
      for (i=1; i<=exCount[k]; i++) {
        if (i>1) exs = exs "; "
        exs = exs exList[k, i]
      }
    }

    print combo, svt, funcv, count[k], med, top, exs >> out1
  }
}

# ---------- helpers ----------
function infoget(s, key,    m,t) {
  if (match(s, "(^|;)" key "=[^;]+")) {
    t = substr(s, RSTART, RLENGTH); sub(/^[^=]*=/,"",t); return t
  }
  return ""
}
function sample_state(field, pos,    a,gt) {
  n = split(field, a, ":")
  gt = (pos["GT"] ? a[pos["GT"]] : "")
  gsub(/\|/, "/", gt)
  if (gt=="" || gt=="." || gt=="./.") return "NA"
  if (gt ~ /1/) return "ALT"
  if (gt == "0/0") return "REF"
  return "NA"
}
function abs(x){ return (x<0)?-x:x }