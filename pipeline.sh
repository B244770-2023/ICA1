#doing the fastq jobs
gzip -dc /localdisk/data/BPSM/ICA1/fastq/*.gz >> ${HOME}/BPSM/ICA1/fastq/
cd ${HOME}/BPSM/ICA1/fastq
fastq *.fq

#doing the bowtie2 jobs
cd ${HOME}/BPSM/ICA1
mkdir bowtie2
#unzip the fasta file
gzip -dc /localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz >> ${HOME}/BPSM/ICA1/bowtie2/tcongo_genome.fasta
#constructing index for bowtie2
bowtie2-build tcongo_genome.fasta reference_index

#using hisat2 for alignment of RNA seqences
cd ${HOME}/BPSM/ICA1
mkdir hisat2
#unzip the fasta file
gzip -dc /localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz >> ${HOME}/BPSM/ICA1/hisat2/tcongo_genome.fasta
#constructing index for hisat2
hisat2-build tcongo_genome.fasta reference_index
#executing hisat2 alignment with loop
for file in ${HOME}/BPSM/ICA1/fastq
do
	base_name=${basename "$file" .fq}
	hisat2 -x reference_index -U ${file} -S ${HOME}/BPSM/ICA1/hisat2/${base_name}.sam
done

#converting to bam files
for file in ${HOME}/BPSM/ICA1/hisat2
do
	base_name=${basename "$file" .sam}
	samtools view -bS ${file} > ${HOME}/BPSM/ICA1/hisat2/${base_name}.bam
done

#counting the number of reads that aligned to the coding region
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ${HOME}/BPSM/ICA1/ttt.bed
cd ${HOME}/BPSM/ICA1


