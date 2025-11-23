#!/usr/bin/env python3
"""
De novo mutation analysis script
"""
import pandas as pd
import re

# File paths
input_vcf = '/data/bioinf-m2_group-project/HG00512-514_PASS_variants.vcf'
output_file = '/data/bioinf-m2_group-project/denovo_mutations.tsv'

# Simple function to parse the format string
def get_data(fmt_str, sample_str):
    keys = fmt_str.split(':')
    vals = sample_str.split(':')
    data = {}
    # Map keys to values
    for i in range(len(keys)):
        if i < len(vals):
            data[keys[i]] = vals[i]
    return data

def main():
    results = []
    
    print("Starting analysis...")
    print(f"Reading from: {input_vcf}")
    
    # Filter thresholds
    MIN_GQ = 20
    MIN_RC = 20
    
    with open(input_vcf, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines
            if line.startswith('#'):
                continue
                
            cols = line.split('\t')
            
            # Basic info
            chrom = cols[0]
            pos = cols[1]
            rsid = cols[2]
            ref = cols[3]
            alt = cols[4]
            flt = cols[6]
            info_txt = cols[7]
            fmt = cols[8]
            
            # Only keep PASS variants
            if flt != 'PASS':
                continue

            # Parse data for the trio
            # Column indices: 9=Father, 10=Mother, 11=Child
            dad_str = cols[9]
            mom_str = cols[10]
            kid_str = cols[11]
            
            dad = get_data(fmt, dad_str)
            mom = get_data(fmt, mom_str)
            kid = get_data(fmt, kid_str)
            
            # Check Genotypes (GT)
            # Parents must be 0/0, child should be 0/1 or 1/1
            d_gt = dad.get('GT', './.')
            m_gt = mom.get('GT', './.')
            k_gt = kid.get('GT', './.')
            
            if d_gt == '0/0' and m_gt == '0/0' and (k_gt == '0/1' or k_gt == '1/1'):
                
                # Quality Control
                try:
                    # Check Genotype Quality (GQ)
                    if int(dad.get('GQ', 0)) < MIN_GQ: continue
                    if int(mom.get('GQ', 0)) < MIN_GQ: continue
                    if int(kid.get('GQ', 0)) < MIN_GQ: continue
                    
                    # Check Read Count (RC)
                    if int(dad.get('RC', 0)) < MIN_RC: continue
                    if int(mom.get('RC', 0)) < MIN_RC: continue
                    if int(kid.get('RC', 0)) < MIN_RC: continue
                    
                    # Check Allele Balance (AB) for the child
                    # For het calls, we expect AB around 0.5
                    k_ab = 0
                    k_ref = 0
                    k_alt = 0
                    
                    # Try to find split reads first (RR/RV)
                    if 'RR' in kid and 'RV' in kid:
                        k_ref = int(kid['RR'])
                        k_alt = int(kid['RV'])
                    # If not, try paired reads (DR/DV)
                    elif 'DR' in kid and 'DV' in kid:
                        k_ref = int(kid['DR'])
                        k_alt = int(kid['DV'])
                        
                    total_reads = k_ref + k_alt
                    
                    if total_reads > 0:
                        k_ab = k_alt / total_reads
                    else:
                        # No reads to support calculation
                        continue
                        
                    # Filter based on AB if heterozygous
                    if k_gt == '0/1':
                        if k_ab < 0.2 or k_ab > 0.8:
                            continue

                    # Extract SV TYPE from INFO field
                    svtype = "Unknown"
                    if 'SVTYPE=' in info_txt:
                        # Use regex to find the type
                        match = re.search(r'SVTYPE=([^;]+)', info_txt)
                        if match:
                            svtype = match.group(1)

                    # Store the valid result
                    res = {
                        'CHROM': chrom,
                        'POS': pos,
                        'ID': rsid,
                        'REF': ref,
                        'ALT': alt,
                        'SV_TYPE': svtype,
                        'Father_GT': d_gt,
                        'Mother_GT': m_gt,
                        'Child_GT': k_gt,
                        'Father_GQ': dad.get('GQ'),
                        'Mother_GQ': mom.get('GQ'),
                        'Child_GQ': kid.get('GQ'),
                        'Father_RC': dad.get('RC'),
                        'Mother_RC': mom.get('RC'),
                        'Child_RC': kid.get('RC'),
                        'Child_Allele_Balance': k_ab,
                        'INFO': info_txt
                    }
                    
                    # Keep raw counts just in case
                    res['Child_RR'] = kid.get('RR', 0)
                    res['Child_RV'] = kid.get('RV', 0)
                    res['Child_DR'] = kid.get('DR', 0)
                    res['Child_DV'] = kid.get('DV', 0)
                    
                    results.append(res)
                    
                except ValueError:
                    # Skip if any conversion fails
                    continue

    # Save to file
    df = pd.DataFrame(results)
    if not df.empty:
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Done! Found {len(df)} mutations.")
        print(f"Saved results to {output_file}")
        
        # Print some stats
        print("\nCounts by SV Type:")
        print(df['SV_TYPE'].value_counts())
        
        print("\nTop Chromosomes:")
        print(df['CHROM'].value_counts().head())
        
        # Check average quality
        avg_gq = df['Child_GQ'].astype(float).mean()
        print(f"\nAverage Child GQ: {avg_gq:.2f}")
        
    else:
        print("No mutations found with these filters.")

if __name__ == "__main__":
    # Run the script
    print("--- De Novo Mutation Analysis ---")
    main()
    # print("Analysis complete.")