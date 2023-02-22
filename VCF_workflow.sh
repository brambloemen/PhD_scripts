#! /bin/bash
ml load bcftools/1.9 trimmomatic bowtie2 samtools/1.9

bowtie2-build Reference/Bvel10075.fasta Reference/Bvel10075
samtools faidx Reference/Bvel10075.fasta
samtools dict -o Reference/Bvel10075.dict Reference/Bvel10075.fasta

mkdir Alignment
for file in $(ls ../Isolate_reads/*.fastq)
do
    acn=$(echo $file | sed -e s/.fastq// | sed 's#.*/##')
    echo $acn $file
    bowtie2 -x Reference/Bvel10075 -p 30 --rg-id $acn -U $file -S Alignment/$acn.sam 

    samtools view -@30 -b -h -o Alignment/$acn.bam Alignment/$acn.sam
    samtools sort -@30 -o Alignment/$acn.sort.bam Alignment/$acn.bam
    mv Alignment/$acn.sort.bam Alignment/$acn.bam
    samtools index -@30 Alignment/$acn.bam
    
# done
cd Alignment
bcftools mpileup --threads 20 -Ou -f ../Reference/Bvel10075.fasta $(ls *.bam) | bcftools call --threads 20 -vmO z -o Protease1_Bvel10075.vcf.gz

tabix -p vcf Protease1_Bvel10075.vcf.gz

bcftools stats --threads 30 -F ../Reference/Bvel10075.fasta -s - Protease1_Bvel10075.vcf.gz Protease1_Bvel10075.vcf.gz.stats
mkdir plots
plot-vcfstats -p plots/ Protease1_Bvel10075.vcf.gz.stats