for file in ${HOME}/BPSM/ICA1/fastq/*
do
	base_name=$(basename $file .fq)
	hisat2 -x reference_index -U ${file} -S ${HOME}/BPSM/ICA1/hisat2/${base_name}.sam
done


