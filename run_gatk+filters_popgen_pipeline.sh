
#all the previous steps, including adapter removal and mapping, followed the script provided in A standard pipeline for processing short-read sequencing data from Littorina snails V.3, by James Reeve et al.


#Personally, I recomend parallelizing this to process a few chromosomes (if not all) at a time

conda activate gatk

REF="path/to/ref.fasta"
BAMS="path/to/folder/with/bams"
RG_ADDED="path/to/empty/folder"
FINAL_VCF_DIR="output/folder"
LOCAL_TMP="/tmp/vcfs"

mkdir -p "$RG_ADDED"
mkdir -p "$LOCAL_TMP"

i=0

for file in "$BAMS"/*_nodup.bam
do
echo "Procssing: $file"

sample=$(basename "$file" _nodup.bam)
out_RG="$RG_ADDED/${sample}_addedRG.bam"
out_vcf="$LOCAL_TMP/${sample}_${CHR}.g.vcf.gz"

lib_name="LIB-${sample}-01"
((i++))
unit_name="unit$i"

    picard AddOrReplaceReadGroups \
     I="$file" \
    O="$out_RG" \
   RGLB="$lib_name" \
  RGPL=ILLUMINA \
 RGPU="$unit_name" \
RGSM="$sample"

    samtools index "$out_RG"

gatk --java-options "-Xmx30g" HaplotypeCaller \
-R "$REF" \
-I "$out_RG" \
-O "$out_vcf" \
-L "$GATK_REGION" \
-ERC GVCF \
--min-base-quality-score 13 \
--minimum-mapping-quality 30 \
--mapping-quality-threshold-for-genotyping 30 \
--pcr-indel-model NONE

mv "$LOCAL_TMP/${sample}_${CHR}.g.vcf.gz" "$FINAL_VCF_DIR/"

done

echo "Calling individual done"


for file in "$BAMS"/*_nodup.bam; do
sample=$(basename "$file" _nodup.bam)

gvcf_path="$FINAL_VCF_DIR/${sample}_${CHR}.g.vcf.gz"

echo "$gvcf_path" >> "$FINAL_VCF_DIR/gvcfs_${CHR}.list"

if [ ! -f "${gvcf_path}.tbi" ]; then
gatk IndexFeatureFile -I "$gvcf_path"
fi
done

gatk --java-options "-Xmx35g" CombineGVCFs \
-R "$REF" \
-V "$FINAL_VCF_DIR/gvcfs_${CHR}.list" \
-O "$FINAL_VCF_DIR/${CHR}.g.vcf.gz"

gatk IndexFeatureFile -I "$FINAL_VCF_DIR/${CHR}.g.vcf.gz"

gatk --java-options "-Xmx35g" GenotypeGVCFs \
-R "$REF" \
-V "$FINAL_VCF_DIR/${CHR}.g.vcf.gz" \
-O "$FINAL_VCF_DIR/${CHR}.vcf.gz"



gatk --java-options "-Xmx16g" GatherVcfs \
-V "path/to/cromosome1.vcf.gz" \
...
# add all chromosomes
-O "final_unfilt.vcf.gz"

gatk IndexFeatureFile -I "final_unfilt.vcf.gz"



#--------filters

bcftools view -Ou -m 2 -M 2 -v snps final_unfilt |
  bcftools filter -Oz -g 5:indel,other \
-o filt_SNP.vcf.gz

bcftools filter -Oz -e 'INFO/DP>2.5xAverage Depth' \
-o filt_COV.vcf.gz filt_SNP.vcf.gz

bcftools filter -Oz -e 'MQ<30' \
-o filt_MQ.vcf.gz filt_COV.vcf.gz

bcftools filter -Oz -e 'QUAL<30' \
-o filt_QUAL.vcf.gz filt_MQ.vcf.gz

bcftools filter -Oz -e 'FS>60.0 || SOR>3.0' \
-o filt_strand.vcf.gz filt_QUAL.vcf.gz

bcftools filter -Ou -S . -e 'FMT/GQ<20' filt_strand.vcf.gz |
  bcftools filter -Oz -S . -e 'FMT/DP<3' \
-o filt_soft.vcf.gz

bcftools filter -Oz -e 'F_MISSING>0.2' \
-o filt.vcf.gz filt_soft.vcf.gz

#remove individuals with missing >0.2

#adapt the chromosome name
bcftools view -Oz -r 'ENA|OZ249416|OZ249416.1' \
-o chromosome11.vcf.gz filt.vcf.gz

