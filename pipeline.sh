#!/bin/bash

#file paths
fqs="/localdisk/data/BPSM/ICA1/fastq"
genome="/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
genes="/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"

#doing the fastq jobs
gzip -dc $fqs/*.gz >> ${HOME}/BPSM/ICA1/fastq/
cd ${HOME}/BPSM/ICA1/fastq
fastq *.fq

#copy the sample details
cp $fqs/Tco2.fqfiles ${HOME}/BPSM/ICA1/

#check the quality of each gene
for file in ${HOME}/BPSM/ICA1/htmls/*.html
do
	#check if there is any failure in base sequence content
	if grep -q '\[FAIL\].*Per base sequence content' ${file}; then
    echo "${file} There might be adapter contamination in the sequences."
	fi
done


#using hisat2 for alignment of RNA seqences
cd ${HOME}/BPSM/ICA1
mkdir hisat2
##unzip the fasta file
gzip -dc $genome >> ${HOME}/BPSM/ICA1/hisat2/tcongo_genome.fasta
##constructing index for hisat2
cd ${HOME}/BPSM/ICA1/hisat2
hisat2-build tcongo_genome.fasta reference_index
##executing hisat2 alignment with loop
for file in ${HOME}/BPSM/ICA1/fastq/*.fq
do
	base_name=$(basename "$file" .fq)
	hisat2 -x reference_index -U ${file} -S ${HOME}/BPSM/ICA1/hisat2/${base_name}.sam
done

#converting sam files to bam files
for file in ${HOME}/BPSM/ICA1/hisat2/*.sam
do
	base_name=$(basename "$file" .sam)
	samtools view -b ${file} > ${HOME}/BPSM/ICA1/hisat2/${base_name}.bam
done

#sorting the bam files
for file in ${HOME}/BPSM/ICA1/hisat2/*.bam
do
	base_name=$(basename "$file" .bam)
	samtools sort ${file} -o sorted_${base_name}.bam
done

#counting the number of reads that aligned to the coding region of the genome
cd ${HOME}/BPSM/ICA1
mkdir counting
cd ${HOME}/BPSM/ICA1/counting
cp $genes ./genes.bed

#generate counting results for each gene
for file in ${HOME}/BPSM/ICA1/hisat2/sorted_*.bam
do
	base_name=$(basename "$file" .bam)
	bedtools coverage -a genes.bed -b ${file} -counts > ${base_name}.coveredoutput
done

#sum up to get counts for the whole genome
for file in ${HOME}/BPSM/ICA1/counting/*.coveredoutput
do
    base_name=$(basename "$file" .coveredoutput)
    realname="${base_name#sorted_}"
    awk '{counts += $NF} END {print counts}' "$file" | sed "s/^/${realname} /" >> counting.txt
done

#calculate mean counts for each gene
#cut off verbose info first
for file in ${HOME}/BPSM/ICA1/counting/*.coveredoutput
do
    base_name=$(basename "$file" .coveredoutput)
    realname="${base_name#sorted_}"
	awk -F'\t' '{print $4, $5, $6}' $file > $realname.count
done

#sum up counts by gene to get expression level
awk '
{
    sum[$2] += $NF;
    line[$2] = $0
}
END {
    for (id in sum) {
        sub(/[0-9]+$/, sum[id], line[id]);
        print line[id]
    }
}
' *.count > gene_expression_level.txt

#



