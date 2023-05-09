
## rename all of the fastq files to be lowercase

rename 'y/A-Z/a-z/' *

## move fastq files to folders dependent on species

while read i; do echo "$i"; done <gambiae_samples.txt

while read i; do mv ${i}* "an_gambiae_fastq"; done <gambiae_samples.txt

while read i; do mv ${i}* "an_coluzzii_fastq"; done <coluzzii_samples.txt

while read i; do mv ${i}* "an_melas_fastq"; done <melas_samples.txt

while read i; do mv ${i}* "an_hybrid_fastq"; done <hybrid_samples.txt
