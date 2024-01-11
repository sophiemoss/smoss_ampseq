# anopheles_ampseq

This repository is for processing amplicon fastq files of anopheles mosquitoes to identify variants including SNPs and INDELs

This combines code from amp_seq with additional code for processing VCF files

First, the raw fastq files need to be demultiplexed.

Then, the amplicon pipeline needs to be run to create a combined VCF using the variant caller of choice (in my case both gatk and freebayes)

Then, this combined VCF should be filtered and annotated using snpeff (sophie_makecombinedvcf_snpeff.py)

Finally, the txt files of annotated SNPs and INDELs can be filtered further using sophie_combined_vcf_filtering.py

The parameters and metadata need to change depending on the sample set.


