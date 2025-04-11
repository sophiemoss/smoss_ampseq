## need to change chromosomes in reference to AgamP4_2L etc in the reference genome

# input GCF_000005575.2_AgamP3_genomic.fasta
# output GCF_000005575.2_AgamP3_genomic_renamedchrs.fasta

awk '{gsub(/NT_078265.2/, "AgamP4_2L"); gsub(/NC_004818.2/, "AgamP4_X"); gsub(/NT_078266.2/, "AgamP4_2R"); gsub(/NT_078267.5/, "AgamP4_3L"); gsub(/NT_078268.4/, "AgamP4_3R"); gsub(/NC_002084.1/, "AgamP4_MT"); print;}' GCF_000005575.2_AgamP3_genomic.fasta > GCF_000005575.2_AgamP3_genomic_renamedchrs.fasta

