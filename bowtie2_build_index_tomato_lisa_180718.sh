#PBS -S /bin/bash
#PBS -lnodes=1:ppn=16
#PBS -lwalltime=0:45:00
	# if the script does not finish in 8hrs, give more hours for calculations
	# (up to 120:00:00 for 5 days)

#bowtie2 version 2.2.4

module load bowtie/2.2.4

#copy to scratch
cp "${HOME}"/genomes/S_lycopersicum_chromosomes.3.00.fa "${TMPDIR}"/

#build index from the reference genome
bowtie2-build S_lycopersicum_chromosomes.3.00.fa tomato_inx_3.00
# bowtie-build [reference fasta] [base name for index file]

mv "${TMPDIR}"/* "${HOME}"/genomes/
