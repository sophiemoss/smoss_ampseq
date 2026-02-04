#####______ Doing an additional filter due to multiplexing, to check whether it make a big difference to the data ______#####

# In the normal amplicon script, this happens:

    run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.genotyped.vcf.gz | "
            "bcftools view --threads 20 -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o tmp.vcf.gz" % vars(args))
    
    log(f"Finished filtering for FMT/DP>10")

# Now we need to edit the genotypes using a python script which is NOT FOR PHASED DATA, IT DOES NOT DISTINGUISH BETWEEN 0/1 and 1/0.
# could have the below as a separate python script but kept here for ease to follow this workflow.

import pysam

infile = "tmp.vcf.gz"
outfile = "combined.genotypes_redone_25_75.vcf.gz"

vcf_in = pysam.VariantFile(infile)
vcf_out = pysam.VariantFile(outfile, 'w', header=vcf_in.header)

for record in vcf_in.fetch():
    for sample in record.samples:
        sample_data = record.samples[sample]
        if sample_data.get('GT') == (None, None):
            continue  # keep as ./.
        
        ad = sample_data.get('AD')
        if ad and len(ad) >= 2:
            ref_count = ad[0]
            alt_count = ad[1]
            total = ref_count + alt_count
            if total == 0:
                continue  # ambiguous, leave as is

            ref_frac = ref_count / total
            alt_frac = alt_count / total

            if ref_frac >= 0.75:
                sample_data['GT'] = (0, 0)
            elif alt_frac >= 0.75:
                sample_data['GT'] = (1, 1)
            else:
                sample_data['GT'] = (0, 1)
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()

## Now we have our edited VCF with edited genotypes, we carry on with the other steps from the amplicon script:

    run_cmd("bcftools view --threads 20 -v snps combined.genotypes_redone_25_75.vcf.gz | bcftools csq -p a -f /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa -g /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.56.gff3 -Oz -o snps_regenotyped_25_75.vcf.gz" % vars(args))

    run_cmd("(echo -e 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tDP\tAD\tTBCSQ' && bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%DP\t%AD\t%TBCSQ\n]' snps_regenotyped_25_75.vcf.gz) > combined_regenotyped_25_75_filtered_formatted.snps.trans.txt")   

## Run the 3.compute_SNP_allele_frequencies.py script with the appropriate files
## Take the output combined_regenotyped_25_75_filtered_formatted.snps.trans.txt and snp_frequencies_filtered_by_coverage_for_regenotyped_25_75.tsv 
## to look at the data.