### Read Trimming ###
java -jar path/to/trimmomatic.jar PE .R1.fq.gz input.R2.fq.gz output.R1paired.fq.gz sample.R1unpaired.fq.gz pool.R2paired.fq.gz pool.R2unpaired.fq.gz ILLUMINACLIP:/path/to/Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:20:10:1:true TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:30 TOPHRED33

### Mapping ###
bwa index path/to/reference.fasta
bwa mem -M -t 8 path/to/reference.fasta pool.R1paired.fq.gz pool.R2paired.fq.gz | samtools view  -Sb - | samtools sort - pool.bam && samtools index pool.bam

### Duplicate removal ###
java -jar path/to/picard.jar MarkDuplicates I=pool.bam O=pool_nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=pool_nodup.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 && samtools index pool_nodup.bam

### Estimate coverage ### 
samtools depth -aa pool_nodup.bam | awk '{sum+=$3; cnt++} END {if (cnt > 0) print "Average coverage:", sum/cnt}'

### Preparing data for popoolation1 ###
#Note: This step will create one pileup file for the respective pool
samtools mpileup -B -f -aa path/to/reference.fasta -Q30 -q30 pool_nodup.bam | grep -v mpileup | awk '$4 >9' > file_per_pool.pileup

### Calculate pi and Tajima’s D (popoolation1) ###
#Note: adjust the maximum coverage (2.5x the average coverage) and pool size (2x pool size for diploid samples)
perl path/to/Variance-sliding.pl --input file_per_pool.pileup --output file.pi --measure pi --window-size 50000 --step-size 50000 --min-covered-fraction 0.2 --min-count 2 --min-coverage 10 --max-coverage 75 --min-qual 30 --pool-size 96 --fastq-type sanger
perl path/to/Variance-sliding.pl --input file_per_pool.pileup --output file.D --measure D --window-size 50000 --step-size 50000 --min-covered-fraction 0.2 --min-count 2 --min-coverage 10 --max-coverage 75 --min-qual 30 --pool-size 96 --fastq-type sanger

### Preparing data for popoolation2 ###
#Note: This step will create one single file for both pools
samtools mpileup -B -f -aa path/to/reference.fasta -Q30 -q30 pool1_nodup.bam pool2_nodup.bam | grep -v mpileup | awk '$4 >4 && $7 >4' > pool1_pool2.mpileup
java -ea -Xmx7g -jar path/to/mpileup2sync.jar --input pool1_pool2.mpileup --output pool1_pool2.sync --fastq-type sanger --min-qual 30 --threads 15

### FST and Allele frequency differences (popoolation2)
#Note: adjust the maximum coverage (2.5x the average coverage) and pool size (number of individuals of both pools)
perl path/to/fst-sliding.pl --input pool1_pool2.sync --output pool1_pool2.fst --min-count 3 --min-coverage 10 --max-coverage 93,93 --min-covered-fraction 0.2 --window-size 50000 --step-size 12500 --pool-size 48:48
perl path/to/snp-frequency-diff.pl --input pool1_pool2.sync --output-prefix pool1_pool2 --min-count 3 --min-coverage 10 --max-coverage 93,93
