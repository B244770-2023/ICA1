#!/bin/bash

#file paths
fqs="/localdisk/data/BPSM/ICA1/fastq"
genome="/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
genes="/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"

cd ${HOME}/BPSM/ICA1
#copy and process the sample details(drop the suffix of End1 and End2)
cp $fqs/Tco2.fqfiles ${HOME}/BPSM/ICA1/
awk 'BEGIN {OFS="\t"} {gsub(/\.gz/, "", $6); gsub(/\.gz/, "", $7); print}' Tco2.fqfiles > Tco2.details

#divide fq files into groups(naming method: time-treatment-sampletype)
mkdir fastq
gzip -dc $fqs/*.gz >> ${HOME}/BPSM/ICA1/fastq/
cd ${HOME}/BPSM/ICA1/fastq
while IFS=$'\t' read -r SampleName SampleType Replicate Time Treatment End1 End2
do
    if [ "$SampleName" == "SampleName" ]; then
        continue
    fi
    foldername="${Time}_${Treatment}_${SampleType}"
    if [ ! -d "$foldername" ]; then
        mkdir "$foldername"
    fi
    mv "$End1" "$foldername/"
    mv "$End2" "$foldername/"
done < ${HOME}/BPSM/ICA1/Tco2.details


#doing the fastq jobs
cd ${HOME}/BPSM/ICA1/fastq
fastq *.fq

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


#!!those actions have to done by hand in the different directory for a few times!!#
##executing hisat2 alignment with loop
wdir="${HOME}/BPSM/ICA1/fastq/0_Uninduced_Clone1" #change directroy here for group-wise
for file in ${HOME}/BPSM/ICA1/fastq/*.fq 
do
	base_name=$(basename "$file" .fq)
	hisat2 -x reference_index -U ${file} -S ${HOME}/BPSM/ICA1/hisat2/${base_name}.sam
done

#converting sam files to bam files
for file in $wdir/*.sam
do
	base_name=$(basename "$file" .sam)
	samtools view -b ${file} > ${HOME}/BPSM/ICA1/hisat2/${base_name}.bam
done

#sorting the bam files
for file in wdir/*.bam
do
	base_name=$(basename "$file" .bam)
	samtools sort ${file} -o sorted_${base_name}.bam
done

#counting the number of reads that aligned to the coding region of the genome
cd $wdir
mkdir counting
cd $wdir/counting
cp $genes ./genes.bed

#generate counting results for each gene
for file in $wdir/sorted_*.bam
do
	base_name=$(basename "$file" .bam)
	bedtools coverage -a genes.bed -b ${file} -counts > ${base_name}.coveredoutput
done

#sum up to get counts for the whole genome
for file in $wdir/*.coveredoutput
do
    base_name=$(basename "$file" .coveredoutput)
    realname="${base_name#sorted_}"
    #output a file which calculates number of reads that align to the regions of coding genes in genome
    awk '{counts += $NF} END {print counts}' "$file" | sed "s/^/${realname} /" >> counting.txt
done

#cut off verbose info, left only gene names, descriptions and counts
for file in ${PWD}/BPSM/ICA1/counting/*.coveredoutput
do
    base_name=$(basename "$file" .coveredoutput)
    realname="${base_name#sorted_}"
	awk -F'\t' '{print $4, $5, $6}' $file > $realname.count
done

#calculate average expression level of genes by groups
awk '
{
    sum[$2] += $NF;
    count[$2]++;
    line[$2] = $0;
}
END {
    for (id in sum) {
        avg = sum[id] / count[id];
        sub(/[0-9]+$/, avg, line[id]);
        print line[id];
    }
}
' *.count > gene_average_expression_level.txt

#calculate fold changes(over time relative to uninduced controls)
#e.g. choosing induced clone1 at 24h and uninduced clone1 at 24h

awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/24_Induced_Clone1/gene_average_expression_level.txt $wir/24_Uninduced_Clone1/gene_average_expression_level.txt > $wir/24_Induced_Clone1/fold_change.txt
awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/48_Induced_Clone1/gene_average_expression_level.txt $wir/48_Uninduced_Clone1/gene_average_expression_level.txt > $wir/48_Induced_Clone1/fold_change.txt
awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/24_Induced_Clone2/gene_average_expression_level.txt $wir/24_Uninduced_Clone2/gene_average_expression_level.txt > $wir/24_Induced_Clone2/fold_change.txt
awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/48_Induced_Clone2/gene_average_expression_level.txt $wir/48_Uninduced_Clone2/gene_average_expression_level.txt > $wir/48_Induced_Clone2/fold_change.txt
awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/24_Induced_WT/gene_average_expression_level.txt $wir/24_Uninduced_WT/gene_average_expression_level.txt > $wir/24_Induced_WT/fold_change.txt
awk 'NR==FNR{a[NR]=$NF; next} {print $0, $NF/a[FNR]}' $wir/48_Induced_WT/gene_average_expression_level.txt $wir/48_Uninduced_WT/gene_average_expression_level.txt > $wir/48_Induced_WT/fold_change.txt
#...