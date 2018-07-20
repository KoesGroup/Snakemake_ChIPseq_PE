#!/usr/bin/env python
#Author :  Tijs Bliek
#20/07/18

import os
import glob

lijst = glob.glob("*.bam")

for f in lijst:
    f = f.split(".bam")[0]
    # take out the proper mapped pairs (sam flags -f 2)
    command1  = "samtools view -b -hf 0x02 " + f + ".bam > " + f + "_mapped.bam"
    # bam to sam
    command2  = "samtools view -h -o " + f + "_mapped.sam " + f + "_mapped.bam"
    # only take the mapping up to 2 mismatches
    command3  = "samtools view -Sh " + f + "_mapped.sam | grep -e \"^@\" -e \"XM:i:[012][^0-9]\" > " + f + "_mapped_lowmiss.sam"
    # create file containing the sam-header
    command4  = "sed -n 1,15p " + f + "_mapped_lowmiss.sam > " + f + "_mapped_lowmiss_multi_header.sam"
    # only take the paires of wich both uniquely map, no secondairy map, indicated with XS:i: INT for mapping quality of the next best
    command5  = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss.sam | grep -v \"XS:i:\" | python WriteToSam.py > " + f + "_mapped_lowmiss_unique.txt"
    # concatanate header with uniquely mapping paires with no more then 2 mismatches
    command6  = "cat " + f + "_mapped_lowmiss_multi_header.sam " + f + "_mapped_lowmiss_unique.txt > " + f + "mapped_lowmiss_unique_ori.sam"
    # filter out the ones with highest mapping quality (=42)
    command7  = "samtools view -h -q 42 " + f + "_mapped_lowmiss_unique_ori.sam > " + f + "_mapped_lowmiss_unique.sam"
    # sam to bam
    command8  = "samtools view -bS -o " + f + "_mapped_lowmiss_unique.bam " + f + "_mapped_lowmiss_unique.sam"
    # sort the file to position
    command9  = "samtools sort " + f + "_mapped_lowmiss_unique.bam -o " + f + "_mapped_lowmiss_unique_sort.bam"
    # index the bam file, creates a bai file
    command10 = "samtools index " + f + "_mapped_lowmiss_unique_sort.bam"
    # back to the mapped.bam from command 1, now to take out the multi-mapped to sam (without the header)
    command11 = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss.sam | grep \"XS:i:\"| python WriteToSam.py > " + f + "_mapped_lowmiss_multi.txt"
    # add the header
    command12 = "cat " + f + "_mapped_lowmiss_multi_header.sam " + f + "_mapped_lowmiss_multi.txt > " + f + "_mapped_lowmiss_multi.sam"
    # take the ones with mapping quality of 10 or more
    command13 = "samtools view -h -q 10 " + f + "_mapped_lowmiss_multi.sam > " + f + "_mapped_lowmiss_multi_fq10.sam"
    # take the proper paired multimapping to sam (without header)
    command14 = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss_multi_fq10.sam | grep \"XS:i:\" | python WriteToSam.py > " + f + "_mapped_lowmiss_multi_fq10.txt"
    # create text file containing the number of lines
    command15 = "wc -l " + f + "_mapped_lowmiss_multi_fq10.txt >> " + f + ".stats 2>&1"
    # take the best mapping of the multi-mapping, in case of multiple best take one at random
    command16 = "Rscript multi_unique_extract_pairend.r " + f + "_mdim.stats " + f + "_mapped_lowmiss_multi_fq10.txt " + f + "_MU.RData " + f + "_mapped_lowmiss_multi_fq10_unique.txt"
    # add the header
    command17 = "cat " + f + "_mapped_lowmiss_multi_header.sam  " + f + "_mapped_lowmiss_multi_fq10_unique.txt > " + f + "_mapped_lowmiss_multi_fq10_unique.sam"
    # sam to bam
    command18 = "samtools view -bS -o " + f + "_mapped_lowmiss_multi_fq10_unique.bam " + f + "_mapped_lowmiss_multi_fq10_unique.sam"
    # sort the bam
    command19 = "samtools sort " + f + "_mapped_lowmiss_multi_fq10_unique.bam -o " + f + "_mapped_lowmiss_multi_fq10_unique_sort.bam"
    # index the bam to create bai
    command20 = "samtools index " + f + "_mapped_lowmiss_multi_fq10_unique_sort.bam"
    # join multimapping with uniquly mapping (from command9)
    command21 = "samtools merge  " + f + "_mapped_lowmiss_unique_both.bam " + f + "_mapped_lowmiss_unique_sort.bam  " + f + "_mapped_lowmiss_multi_fq10_unique_sort.bam"
    # sort the bam
    command22 = "samtools sort " + f + "_mapped_lowmiss_unique_both.bam -o " + f + "_mapped_lowmiss_unique_both_sort.bam"
    # index the bam
    command23 = "samtools index " + f + "_mapped_lowmiss_unique_both_sort.bam"

commands = [command1, command2, command3, command4, command5, command6, command7, command8, command9, command10,\
            command11, command12, command13, command14, command15, command16, command17, command18, command19, command20,\
            command21, command22, command23]
name = f + ".sh"

outFile = open(name, "w")

outFile.write("#!/bin/bash\n\n")
outFile.write("\n".join(commands))

outFile.close()

#os.system("chmod +x " + name)
#os.system("./" + name)
