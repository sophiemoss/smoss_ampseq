# Extract coverage for a particular sequencing pool

for i in `ls *_region_coverage.txt`; do cat $i | awk -v sn=$i '{print sn"\t"$4"\t"$5}'>> allsamples_region_coverage.txt; done

# Read into R
# Transform
# Export to excel
# Colour for regions with coverage under 30 - these need to be amplified again (may need to re-do PCR, or if there was a band, send to genewiz again)

## Can look for depth at a specific position 
samtools depth bu1020.bam | egrep 'AgamP4_2L' | egrep '2391228'
