#!/bin/bash
# Script to clean VCF, filter to trio, summarize variants, and identify de novo mutations
# Steps: 
#   1) Clean VCF (remove LowQual, filter by quality, remove missing genotypes)
#   2) Extract trio samples (HG00512, HG00513, HG00514)
#   3) Summarize variant combinations (location, type, function)
#   4) Find de novo mutations (pedigree-agnostic)
# Output: Single cleaned trio VCF (DellyVariation.trio.clean.vcf.gz)
# De novo detection: Pedigree-agnostic - checks all three samples as potential probands
#                    Identifies variants where exactly one sample has variant and other two are homozygous reference
INPUT="DellyVariation.vcf"
OUTPUT_VCF="DellyVariation.trio.clean.vcf"
DENOVO_VCF="DellyVariation.denovo.vcf"
DENOVO_LIST="DellyVariation.denovo.txt"
SUMMARY_REPORT="DellyVariation.trio.summary.txt"
TEMP_FILE="DellyVariation.temp.vcf"
SAMPLE1="HG00512"
SAMPLE2="HG00513"
SAMPLE3="HG00514"
echo "=========================================="
echo "VCF Cleaning, Trio Filtering and De Novo Detection"
echo "=========================================="
echo ""
echo "Trio samples (pedigree unknown):"
echo "  Sample 1: $SAMPLE1"
echo "  Sample 2: $SAMPLE2"
echo "  Sample 3: $SAMPLE3"
echo ""
echo "Note: Pedigree relationships unknown - checking all three samples as potential probands"
echo ""
echo "Output: Single cleaned trio VCF ($OUTPUT_VCF.gz)"
echo ""
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file '$INPUT' not found!"
    exit 1
fi
# ============================================================================
# PART 1: Clean VCF - Remove low-quality variants
# ============================================================================
echo "=========================================="
echo "PART 1: Cleaning VCF"
echo "=========================================="
echo ""
echo "Step 1: Analyzing original VCF file..."
ORIGINAL_COUNT=$(bcftools view -H $INPUT | wc -l)
LOWQUAL_COUNT=$(bcftools view -H $INPUT | awk '$7 == "LowQual"' | wc -l)
PASS_COUNT=$(bcftools view -H $INPUT | awk '$7 == "PASS"' | wc -l)
echo "  - Total variants: $ORIGINAL_COUNT"
echo "  - PASS variants: $PASS_COUNT"
echo "  - LowQual variants: $LOWQUAL_COUNT"
echo ""
echo "Step 2: Filtering PASS variants only..."
bcftools view -f PASS $INPUT -o $TEMP_FILE
AFTER_PASS=$(bcftools view -H $TEMP_FILE | wc -l)
echo "  - Variants after PASS filter: $AFTER_PASS"
echo "  - Removed: $((ORIGINAL_COUNT - AFTER_PASS)) variants"
echo ""

# Step 2: Filter by minimum read support (DV + RV >= 3)
echo "Step 3: Filtering by minimum read support (DV + RV >= 3)..."
bcftools filter -i '(sum(FMT/DV) + sum(FMT/RV) >= 3)' $TEMP_FILE -o ${TEMP_FILE}.2
mv ${TEMP_FILE}.2 $TEMP_FILE
AFTER_READS=$(bcftools view -H $TEMP_FILE | wc -l)
echo "  - Variants after read support filter: $AFTER_READS"
echo "  - Removed: $((AFTER_PASS - AFTER_READS)) variants"
echo ""

# Step 3: Filter by minimum PE/SR support or MAPQ (PE >= 3 OR SR >= 3 OR MAPQ >= 20)
echo "Step 4: Filtering by PE/SR support or MAPQ..."
bcftools filter -i '(INFO/PE >= 3 || INFO/SR >= 3 || INFO/MAPQ >= 20)' $TEMP_FILE -o ${TEMP_FILE}.clean
mv ${TEMP_FILE}.clean $TEMP_FILE

CLEANED_COUNT=$(bcftools view -H $TEMP_FILE | wc -l)
echo "  - Variants after quality filter: $CLEANED_COUNT"
echo "  - Removed: $((AFTER_READS - CLEANED_COUNT)) variants"
echo ""

echo "Cleaning summary:"
echo "  - Original variants:  $ORIGINAL_COUNT"
echo "  - After cleaning:     $CLEANED_COUNT"
echo "  - Total removed:      $((ORIGINAL_COUNT - CLEANED_COUNT)) ($(echo "scale=2; ($ORIGINAL_COUNT - $CLEANED_COUNT) * 100 / $ORIGINAL_COUNT" | bc)%)"
echo ""
echo "Filtering criteria applied:"
echo "  1. FILTER = PASS (removed LowQual)"
echo "  2. Minimum read support: sum(DV + RV) >= 3"
echo "  3. Quality support: PE >= 3 OR SR >= 3 OR MAPQ >= 20"
echo "  4. Remove variants with missing genotypes (./. or .|.) in any sample"
echo ""

# ============================================================================
# PART 2: Extract trio samples directly into final output
# ============================================================================

echo "=========================================="
echo "PART 2: Extracting Trio Samples to Final Output"
echo "=========================================="
echo ""

# Step 1: Extract only the three samples from cleaned temp VCF directly to final output
echo "Step 5: Extracting trio samples ($SAMPLE1, $SAMPLE2, $SAMPLE3) from cleaned VCF..."
bcftools view -s "$SAMPLE1,$SAMPLE2,$SAMPLE3" $TEMP_FILE -o ${TEMP_FILE}.trio

# Clean up temp file
rm -f $TEMP_FILE

TRIO_COUNT_BEFORE=$(bcftools view -H ${TEMP_FILE}.trio | wc -l)
echo "  - Variants in trio VCF (before missing genotype filter): $TRIO_COUNT_BEFORE"
echo ""

# Step 2: Remove variants with missing genotypes (./. or .|.)
echo "Step 6: Removing variants with missing genotypes (./. or .|.)..."
# Filter: keep only variants where all three samples have non-missing genotypes
# Extract GT field from FORMAT columns and filter using awk
bcftools view -h ${TEMP_FILE}.trio > $OUTPUT_VCF
bcftools view -H ${TEMP_FILE}.trio | \
awk 'BEGIN {FS="\t"; OFS="\t"} {
    # Extract GT fields (first field in FORMAT columns 10, 11, 12)
    split($10, gt1, ":")
    split($11, gt2, ":")
    split($12, gt3, ":")
    # Check if any GT is missing (./. or .|.)
    # Use regex to match ./., .|., or any variant with missing alleles
    if (gt1[1] !~ /\.\/\./ && gt1[1] !~ /\.\|\./ && 
        gt2[1] !~ /\.\/\./ && gt2[1] !~ /\.\|\./ && 
        gt3[1] !~ /\.\/\./ && gt3[1] !~ /\.\|\./) {
        print
    }
}' >> $OUTPUT_VCF

TRIO_COUNT=$(bcftools view -H $OUTPUT_VCF | wc -l)
MISSING_REMOVED=$((TRIO_COUNT_BEFORE - TRIO_COUNT))
echo "  - Variants after removing missing genotypes: $TRIO_COUNT"
echo "  - Removed variants with missing genotypes: $MISSING_REMOVED"
echo ""
rm -f ${TEMP_FILE}.trio

# ============================================================================
# PART 3: Summarize Variant Combinations (Location, Type, Function)
# ============================================================================

echo "=========================================="
echo "PART 3: Summarizing Variant Combinations"
echo "=========================================="
echo ""

echo "Step 7: Analyzing variant combinations, locations, types, and functions..."

# Create comprehensive summary report
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/PRECISE\t%INFO/PE\t%INFO/SR\t%INFO/END\t[%GT]\t[%GT]\t[%GT]\n' \
    $OUTPUT_VCF | \
awk -v s1="$SAMPLE1" -v s2="$SAMPLE2" -v s3="$SAMPLE3" -v summary_file="$SUMMARY_REPORT" -v input_file="$INPUT" \
'BEGIN {
    FS="\t"; OFS="\t"
    total=0
    
    # Pattern counters
    all_00=0; all_var=0; two_var=0; one_12=0; one_13=0; one_14=0; other_pat=0
    
    # SV type counters
    del_total=0; del_all_var=0
    ins_total=0; ins_all_var=0
    dup_total=0; dup_all_var=0
    
    # Chromosome counters
    for(i=1;i<=22;i++) chr_count["chr"i]=0
    chr_count["chrX"]=0; chr_count["chrY"]=0; chr_count["other"]=0
    
    # Precision counters
    precise_count=0; imprecise_count=0
    
    print "CHROM\tPOS\tEND\tID\tSVTYPE\tPRECISE\t" s1 "_GT\t" s2 "_GT\t" s3 "_GT\tPATTERN" > (summary_file ".detailed")
}
{
    total++
    gt1=$9; gt2=$10; gt3=$11
    svtype=$4; chr=$1; end=$8
    precise=$5
    
    # Chromosome counting
    chr_key=chr
    if(chr !~ /^chr[0-9XY]/) chr_key="other"
    chr_count[chr_key]++
    
    # Precision counting
    if(precise == "1") precise_count++
    else imprecise_count++
    
    # Pattern detection
    has_var1=(gt1 ~ /[1-9]|1\/1|1\|1|0\/1|0\|1|1\|0/)
    has_var2=(gt2 ~ /[1-9]|1\/1|1\|1|0\/1|0\|1|1\|0/)
    has_var3=(gt3 ~ /[1-9]|1\/1|1\|1|0\/1|0\|1|1\|0/)
    
    pattern=""
    if (gt1 == "0/0" && gt2 == "0/0" && gt3 == "0/0") {
        pattern="All_homozygous_ref"
        all_00++
    }
    else if (has_var1 && has_var2 && has_var3) {
        pattern="All_three_with_variant"
        all_var++
        if(svtype=="DEL") del_all_var++
        else if(svtype=="INS") ins_all_var++
        else if(svtype=="DUP") dup_all_var++
    }
    else if ((has_var1 && has_var2 && gt3 == "0/0") || 
             (has_var1 && has_var3 && gt2 == "0/0") || 
             (has_var2 && has_var3 && gt1 == "0/0")) {
        pattern="Two_with_variant"
        two_var++
    }
    else if (has_var1 && gt2 == "0/0" && gt3 == "0/0") {
        pattern="Only_" s1
        one_12++
    }
    else if (has_var2 && gt1 == "0/0" && gt3 == "0/0") {
        pattern="Only_" s2
        one_13++
    }
    else if (has_var3 && gt1 == "0/0" && gt2 == "0/0") {
        pattern="Only_" s3
        one_14++
    }
    else {
        pattern="Other"
        other_pat++
    }
    
    # SV type counting
    if(svtype=="DEL") {del_total++; if(!has_var1 && !has_var2 && !has_var3) del_all_var--}
    else if(svtype=="INS") {ins_total++; if(!has_var1 && !has_var2 && !has_var3) ins_all_var--}
    else if(svtype=="DUP") {dup_total++; if(!has_var1 && !has_var2 && !has_var3) dup_all_var--}
    
    # Print to detailed report (first 100 variants as examples)
    if(total <= 100) print chr, $2, end, $3, svtype, precise, gt1, gt2, gt3, pattern >> (summary_file ".detailed")
}
END {
    # Generate summary report
    print "==========================================" > summary_file
    print "TRIO VARIANT COMBINATION SUMMARY REPORT" >> summary_file
    print "==========================================" >> summary_file
    print "" >> summary_file
    cmd = "date"; cmd | getline date_str; close(cmd)
    print "Generated: " date_str >> summary_file
    print "Input VCF: " input_file >> summary_file
    print "Samples: " s1 ", " s2 ", " s3 >> summary_file
    print "" >> summary_file
    print "TOTAL VARIANTS: " total >> summary_file
    print "" >> summary_file
    print "--- GENOTYPE PATTERNS ---" >> summary_file
    printf "All homozygous reference (0/0): %d (%.2f%%)\n", all_00, (all_00/total)*100 >> summary_file
    printf "All three with variant: %d (%.2f%%)\n", all_var, (all_var/total)*100 >> summary_file
    printf "Two samples with variant: %d (%.2f%%)\n", two_var, (two_var/total)*100 >> summary_file
    printf "Only %s with variant: %d (%.2f%%)\n", s1, one_12, (one_12/total)*100 >> summary_file
    printf "Only %s with variant: %d (%.2f%%)\n", s2, one_13, (one_13/total)*100 >> summary_file
    printf "Only %s with variant: %d (%.2f%%)\n", s3, one_14, (one_14/total)*100 >> summary_file
    printf "Other patterns: %d (%.2f%%)\n", other_pat, (other_pat/total)*100 >> summary_file
    print "" >> summary_file
    print "--- BY STRUCTURAL VARIANT TYPE ---" >> summary_file
    print "" >> summary_file
    print "DELETIONS (DEL) - Function: Removes DNA sequence" >> summary_file
    printf "  Total: %d (%.1f%%)\n", del_total, (del_total/total)*100 >> summary_file
    printf "  All with variant: %d (%.1f%% of DEL)\n", del_all_var, (del_all_var/del_total)*100 >> summary_file
    print "" >> summary_file
    print "INSERTIONS (INS) - Function: Adds DNA sequence" >> summary_file
    printf "  Total: %d (%.1f%%)\n", ins_total, (ins_total/total)*100 >> summary_file
    printf "  All with variant: %d (%.1f%% of INS)\n", ins_all_var, (ins_all_var/ins_total)*100 >> summary_file
    print "" >> summary_file
    print "DUPLICATIONS (DUP) - Function: Copies DNA segment" >> summary_file
    printf "  Total: %d (%.1f%%)\n", dup_total, (dup_total/total)*100 >> summary_file
    printf "  All with variant: %d (%.1f%% of DUP)\n", dup_all_var, (dup_all_var/dup_total)*100 >> summary_file
    print "" >> summary_file
    print "--- CHROMOSOMAL DISTRIBUTION (LOCATION) ---" >> summary_file
    print "" >> summary_file
    for (c in chr_count) if (chr_count[c] > 0) printf "  %s: %d variants (%.2f%%)\n", c, chr_count[c], (chr_count[c]/total)*100 >> summary_file
    print "" >> summary_file
    print "--- PRECISION ---" >> summary_file
    printf "  Precise: %d (%.1f%%)\n", precise_count, (precise_count/total)*100 >> summary_file
    printf "  Imprecise: %d (%.1f%%)\n", imprecise_count, (imprecise_count/total)*100 >> summary_file
    print "" >> summary_file
    print "--- KEY FINDINGS ---" >> summary_file
    print "" >> summary_file
    if (one_12 + one_13 + one_14 == 0) {
        print "  * No de novo mutations detected (no variants with single sample variant)" >> summary_file
    } else {
        print "  * Potential de novo mutations found:" >> summary_file
        printf "    - %s: %d variants\n", s1, one_12 >> summary_file
        printf "    - %s: %d variants\n", s2, one_13 >> summary_file
        printf "    - %s: %d variants\n", s3, one_14 >> summary_file
    }
    print "" >> summary_file
    print "==========================================" >> summary_file
    
    # Print summary to stdout
    print "  Summary statistics generated:"
    print "    - Total variants: " total
    print "    - Genotype patterns: " (all_00 + all_var + two_var + one_12 + one_13 + one_14 + other_pat) " categories"
    print "    - SV types: DEL (" del_total "), INS (" ins_total "), DUP (" dup_total ")"
    print "    - Chromosomes: " length(chr_count) " chromosomes with variants"
    print "    - Summary saved to: " "summary_file"
}' 

# Display summary
cat $SUMMARY_REPORT | head -50

echo ""
echo "  - Detailed variant list (first 100) saved to: ${SUMMARY_REPORT}.detailed"
echo ""

# ============================================================================
# PART 4: Identify de novo mutations
# ============================================================================

echo "=========================================="
echo "PART 4: De Novo Mutation Detection"
echo "=========================================="
echo ""

# De novo criteria (pedigree-agnostic):
# For each variant, check if exactly one sample has the variant (0/1 or 1/1)
# and the other two samples are homozygous reference (0/0)
# This identifies potential de novo mutations without assuming which sample is the child

echo "Step 8: Identifying potential de novo mutations..."
echo "  Criteria (checking all three samples as potential probands):"
echo "    - Exactly one sample has variant (0/1 or 1/1)"
echo "    - The other two samples are homozygous reference (0/0)"
echo "    - This identifies variants that could be de novo in any of the three samples"
echo ""

# Create a script to filter de novo variants
# Using bcftools filter with expression to check genotypes
# Format: GT[0] is first allele, GT[1] is second allele
# 0/0 means both are 0 (reference)
# 0/1 means heterozygous
# 1/1 means homozygous alternate

# Extract de novo variants using bcftools query + awk
# This approach is more reliable for genotype checking
echo "  Filtering de novo variants..."

# First, create a list of potential de novo variants using bcftools query and awk
# Check all three combinations: any sample with variant while other two are homozygous reference
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t[%GT]\t[%GT]\t[%GT]\tPE=%INFO/PE;SR=%INFO/SR;MAPQ=%INFO/MAPQ\n' \
    $OUTPUT_VCF | \
awk -v s1="'"$SAMPLE1"'" -v s2="'"$SAMPLE2"'" -v s3="'"$SAMPLE3"'" 'BEGIN {
        FS="\t"
        OFS="\t"
        denovo_count = 0
        print "CHROM\tPOS\tID\tREF\tALT\tFILTER\tSVTYPE\t" s1 "_GT\t" s2 "_GT\t" s3 "_GT\tPROBAND\tQUALITY_METRICS" > "'"$DENOVO_LIST"'.tmp"
      } 
      {
        gt1 = $8
        gt2 = $9
        gt3 = $10
        
        # Function to check if genotype has variant
        function has_variant(gt) {
          return (gt == "0/1" || gt == "1/1" || gt == "0|1" || gt == "1|1" || gt == "1|0")
        }
        
        # Function to check if genotype is homozygous reference
        function is_ref(gt) {
          return (gt == "0/0" || gt == "0|0")
        }
        
        # Check all three scenarios: which sample could be the proband with de novo mutation?
        proband = ""
        
        # Scenario 1: Sample 1 has variant, samples 2 and 3 are homozygous reference
        if (has_variant(gt1) && is_ref(gt2) && is_ref(gt3)) {
          proband = s1
          print $1, $2, $3, $4, $5, $6, $7, gt1, gt2, gt3, proband, $11 >> "'"$DENOVO_LIST"'.tmp"
          denovo_count++
        }
        # Scenario 2: Sample 2 has variant, samples 1 and 3 are homozygous reference
        else if (has_variant(gt2) && is_ref(gt1) && is_ref(gt3)) {
          proband = s2
          print $1, $2, $3, $4, $5, $6, $7, gt1, gt2, gt3, proband, $11 >> "'"$DENOVO_LIST"'.tmp"
          denovo_count++
        }
        # Scenario 3: Sample 3 has variant, samples 1 and 2 are homozygous reference
        else if (has_variant(gt3) && is_ref(gt1) && is_ref(gt2)) {
          proband = s3
          print $1, $2, $3, $4, $5, $6, $7, gt1, gt2, gt3, proband, $11 >> "'"$DENOVO_LIST"'.tmp"
          denovo_count++
        }
      }
      END {
        print denovo_count
      }' > denovo_count.tmp

DENOVO_COUNT=$(cat denovo_count.tmp)
rm -f denovo_count.tmp

echo "  - De novo variants identified: $DENOVO_COUNT"
echo ""

# Step 2: Extract de novo variants from trio VCF
echo "Step 9: Extracting de novo variants to VCF..."

if [ "$DENOVO_COUNT" -gt 0 ]; then
    # Create a regions file for de novo positions
    tail -n +2 "${DENOVO_LIST}.tmp" | awk '{print $1":"$2}' | sort -u > denovo_regions.txt
    
    # Extract variants at these positions from the final trio VCF
    bcftools view -R denovo_regions.txt $OUTPUT_VCF -o $DENOVO_VCF
    
    # Move the de novo list to final location
    mv "${DENOVO_LIST}.tmp" $DENOVO_LIST
    
    # Clean up
    rm -f denovo_regions.txt
    
    echo "  - De novo variants VCF saved to: $DENOVO_VCF"
    echo "  - De novo variants report saved to: $DENOVO_LIST"
else
    echo "  - No potential de novo mutations found."
    mv "${DENOVO_LIST}.tmp" $DENOVO_LIST
    echo -e "No potential de novo mutations found in the trio.\n" >> $DENOVO_LIST
    echo -e "De novo criteria (pedigree-agnostic):\n- Exactly one sample has variant (0/1 or 1/1)\n- The other two samples are homozygous reference (0/0)\n- Checked all three samples ($SAMPLE1, $SAMPLE2, $SAMPLE3) as potential probands" >> $DENOVO_LIST
    rm -f $DENOVO_VCF
    touch $DENOVO_VCF
fi

echo ""

# ============================================================================
# PART 5: Compress and index output files
# ============================================================================

echo ""
echo "=========================================="
echo "PART 5: Compressing and Indexing"
echo "=========================================="
echo ""

echo "Step 10: Compressing and indexing output files..."
bgzip -f $OUTPUT_VCF
bcftools index ${OUTPUT_VCF}.gz
echo "  - Final cleaned trio VCF: ${OUTPUT_VCF}.gz"

if [ "$DENOVO_COUNT" -gt 0 ]; then
    bgzip -f $DENOVO_VCF
    bcftools index ${DENOVO_VCF}.gz
    echo "  - De novo VCF: ${DENOVO_VCF}.gz"
else
    echo "  - De novo VCF: (none found)"
fi

echo ""

# Display summary
echo ""
echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo "Input file:             $INPUT"
echo "Final output VCF:       ${OUTPUT_VCF}.gz (cleaned + trio only)"
echo "Variant summary:        $SUMMARY_REPORT"
echo "De novo VCF:            ${DENOVO_VCF}.gz"
echo "De novo report:         $DENOVO_LIST"
echo ""
echo "Variant counts:"
echo "  Original variants:        $ORIGINAL_COUNT"
echo "  After cleaning:           $CLEANED_COUNT ($(echo "scale=2; ($CLEANED_COUNT) * 100 / $ORIGINAL_COUNT" | bc)%)"
echo "  Final trio variants:      $TRIO_COUNT"
echo "  Potential de novo:        $DENOVO_COUNT"
echo ""
echo "Note: De novo mutations identified assuming unknown pedigree."
echo "      Each variant shows which sample could be the proband (child)."
echo ""

if [ "$DENOVO_COUNT" -gt 0 ]; then
    echo "Potential de novo rate: $(echo "scale=4; $DENOVO_COUNT * 1000000 / $TRIO_COUNT" | bc) per million variants"
    echo ""
    echo "First few potential de novo variants:"
    head -6 $DENOVO_LIST | column -t
    echo ""
    echo "PROBAND column indicates which sample has the variant (potential child)"
fi

echo ""
echo "✓ Done!"

