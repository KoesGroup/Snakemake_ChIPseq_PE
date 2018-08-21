import os
import glob

lijst = glob.glob("*.bam")
n = 1
for f in lijst:
    f = f.split(".bam")[0]
    # take the proper mapped pairs (sam flags -f 2)
    command1  = "samtools view -b -hf 0x02 " + f + ".bam > " + f + "_mapped.bam"
    # bam to sam
    command2  = "samtools view -h -o " + f + "_mapped.sam " + f + "_mapped.bam"
    # only take the mappings up to 2 mismatches(for 75bp reads), 3 was choosen for 150bp reads.
    command3  = "samtools view -Sh " + f + "_mapped.sam | grep -e \"^@\" -e \"XM:i:[0123][^0-9]\" > " + f + "_mapped_lowmiss.sam"
    # remove sam file that is no longer needed
    command3r = "rm " + f + "_mapped.sam"
    # create file containing the sam-header, for tomato the header is the first 15 lines (1,15p), this needs to be changed for different species
    command4  = "sed -n 1,15p " + f + "_mapped_lowmiss.sam > " + f + "_mapped_lowmiss_multi_header.sam"
    # only take the paires of wich both uniquely map, no secondairy map, indicated with XS:i: INT for mapping quality of the next best
    command5  = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss.sam | grep -v \"XS:i:\" | python WriteToSam.py > " + f + "_mapped_lowmiss_unique.txt"
    # concatanate header with uniquely mapping paires with no more then 2 mismatches
    command6  = "cat " + f + "_mapped_lowmiss_multi_header.sam " + f + "_mapped_lowmiss_unique.txt > " + f + "_mapped_lowmiss_unique_ori.sam"
    # remove txt file that is no longer needed
    command6r = "rm " + f + "_mapped_lowmiss_unique.txt"
    # filter out the ones with highest mapping quality (=42)
    command7  = "samtools view -h -q 42 " + f + "_mapped_lowmiss_unique_ori.sam > " + f + "_mapped_lowmiss_unique.sam"
    # remove sam file that is no longer needed
    command7r = "rm " + f + "_mapped_lowmiss_unique_ori.sam"
    # sam to bam
    command8  = "samtools view -bS -o " + f + "_mapped_lowmiss_unique.bam " + f + "_mapped_lowmiss_unique.sam"
    # remove sam file that is no longer needed
    command8r = "rm " + f + "_mapped_lowmiss_unique.sam"
    # sort the file to position
    command9  = "samtools sort " + f + "_mapped_lowmiss_unique.bam -o " + f + "_mapped_lowmiss_unique_sort.bam"
    # index the bam file, creates a bai file, needed for visualisation in programs like IGV
    command10 = "samtools index " + f + "_mapped_lowmiss_unique_sort.bam"
    # back to the mapped.bam from command 1, now to take out the multi-mapped to sam (without the header)
    command11 = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss.sam | grep \"XS:i:\"| python WriteToSam.py > " + f + "_mapped_lowmiss_multi.txt"
    # remove sam file that is no longer needed
    command11r = "rm " + f + "_mapped_lowmiss.sam"
    # add the header
    command12 = "cat " + f + "_mapped_lowmiss_multi_header.sam " + f + "_mapped_lowmiss_multi.txt > " + f + "_mapped_lowmiss_multi.sam"
    # remove txt file that is no longer needed
    command12r = "rm " + f + "_mapped_lowmiss_multi.txt"
    # take the ones with mapping quality of 10 or more
    command13 = "samtools view -h -q 10 " + f + "_mapped_lowmiss_multi.sam > " + f + "_mapped_lowmiss_multi_fq10.sam"
    # remove sam file that is no longer needed
    command13r = "rm " + f + "_mapped_lowmiss_multi.sam"
    # take the proper paired multimapping to sam (without header)
    command14 = "samtools view -S -f 0x02 " + f + "_mapped_lowmiss_multi_fq10.sam | grep \"XS:i:\" | python WriteToSam.py > " + f + "_mapped_lowmiss_multi_fq10.txt"
    # remove sam file that is no longer needed
    command14r = "rm " + f + "_mapped_lowmiss_multi_fq10.sam"
    # create text file containing the number of lines
    command15 = "wc -l " + f + "_mapped_lowmiss_multi_fq10.txt >> " + f + "_mdim.stats 2>&1"
    # take the best mapping of the multi-mapping, in case of multiple best take one at random
    # 8 args are requiered, after the .r file, separated by spaces
    # args[1] is the summary stats file for multi_fq10.sam lines as created in command15
    # args[2] is the input sam file multi_fq10.txt as created in command14
    # args[3] is halfway saved the .RData file for multi-unique
    # args[4] is the final output sam file without header.
    # args[5] skip number, is the total lines of your header in sam file, for tomato 15, arabidopsis, 9
    # args[6] number of lines each time you want to read in R(this is trying to avoid big data problem.), for tomato 20000 is used
    # args[7] total columns in sam file, for single end: samfile has 20 columns, for pair-end: 21
    # args[8] columns for the tag, for single end: 15; pair-end: 16
    command16 = "Rscript multi_unique_extract_Sevgin.r " + f + "_mdim.stats " + f + "_mapped_lowmiss_multi_fq10.txt " + f + "_MU.RData " + f + "_mapped_lowmiss_multi_fq10_unique.txt 15 20000 21 16"
    # remove txt file that is no longer needed
    command16r = "rm " + f + "_mapped_lowmiss_multi_fq10.txt"
    # add the header
    command17 = "cat " + f + "_mapped_lowmiss_multi_header.sam " + f + "_mapped_lowmiss_multi_fq10_unique.txt > " + f + "_mapped_lowmiss_multi_fq10_unique.sam"
    # remove txt file that is no longer needed
    command17r = "rm " + f + "_mapped_lowmiss_multi_fq10_unique.txt"
    # sam to bam
    command18 = "samtools view -bS -o " + f + "_mapped_lowmiss_multi_fq10_unique.bam " + f + "_mapped_lowmiss_multi_fq10_unique.sam"
    # remove sam file that is no longer needed
    command18r = "rm " + f + "_mapped_lowmiss_multi_fq10_unique.sam"
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
    # remove _multy_header.sam, _MU.RData and  _mdim.stats as they are no longer needed.
    command23r = "rm " + f + "_mapped_lowmiss_multi_header.sam " + f + "_MU.RData " + f + "_mdim.stats"
    # this should result in 3 bam-files (might no longer be needed and can be removed), 3 sort.bam files and 3 bam.bai files.

    commands = [command1, command2, command3, command3r, command4, command5, command6, command6r, command7, command7r, command8,\
            command8r, command9, command10, command11, command11r, command12, command12r, command13, command13r, command14,\
            command14r, command15, command16, command16r, command17, command17r, command18, command18r, command19, command20,\
            command21, command22, command23, command23r]
    name = f + ".sh"

    outFile = open(name, "w")
    outFile.write("#!/bin/bash\n\n")
    outFile.write("\n".join(commands))
    outFile.close()

    os.system("chmod +x " + name)
    os.system("screen -dmS " + f + " ./" + name)
    print("sample " + f + ".bam is running in screen " + f)