## After running the fastq2matrix pipeline we can extract information from these files to gain an idea of their quality.
## We are interested in:
## Number of reads
## Percentage of reads mapping to the reference genome
## Coverage
## Number of SNPs

# To extract the number of reads
# For multiple samples, feed in a list of all the samples and use xargs

#code below is saying show the first line (head -n 1) and extract the first item (print $1)
#it creates a list (total_reads) of the number for all samples in the order of the sample list

ls *.mkdup.bamstats | xargs -i -P1 sh -c 'head -1 {} | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > total_reads
#ls *.bqsr.bamstats | xargs -i -P1 sh -c 'head -1 {} | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > total_reads

# Number of reads mapping to reference genome
# Same principle, but extracting the first item on the 7th line

#for some files
ls *.mkdup.bamstats| xargs -i -P1 sh -c 'head -7 {} | tail -1 | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > mapped_reads

# Look at coverage as a %
# Step 1 - Create a list of all the sample names, you just want the prefix

ls *.mkdup.bam | sed -r 's/.mkdup.bam//g' > covlist

# Step 2 - Create a text file with the bash script that you want to run in parallel

#vim coverage_script.sh
#i
#samtools depth $1.mkdup.bam -a > $1.coverage 

# Step 3 - Call the command using parallel - to create coverage files
# Open the sample list, then using that list use parallel to go through the list,
# one by one (-j 1) and execute the bash script coverage_script.sh. The --bar is to show a progress bar.

 cat covlist | parallel -j 25 --bar bash coverage_script.sh {}

 # Step 4 - Extract info from coverage files
 # Now we are interested how many postions in the genome have a coverage of >5 or >10
# Use following script to go through each cov file and extract the number of positions that have a coverage of 5 or above - easily modified to 10
# For mosquitoes could use 10 but Holly is using 5, which is what Susana suggested
# This creates a list in a text file which can easily be copied into a table - divide by genome length for proportion

ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=5)print}'"'"' | wc -l' > coverage_5
ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=10)print}'"'"' | wc -l' > coverage_10
ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=20)print}'"'"' | wc -l' > coverage_20

## Code to look at what the coverage is of wgs data by looking at the bam file and create coverage plots across genome

## Check coverage as an average read depth across entire genome

#samtools depth -aa bamfile.mkdup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'

#example:
#samtools depth -aa bu1002_Combined.mkdup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
#for multiple samples:

ls *.mkdup.bam | parallel --keep-order -j 60 --bar 'samtools depth -aa {} | awk '\''{sum+=$3} END { print sum/NR}'\''' > average_read_depth


