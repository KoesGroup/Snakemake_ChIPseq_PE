import glob
import os

files = glob.glob('*.gz')
print(files)
for f in files:
    #need to add a split function to here the '.gz' at the end of the file name
    command1 = "seqtk sample -s100 " + f + " 10000 > " + "sub_" + f.split(".gz")[0]
    command2 = "echo " + f
    commands = [command1, command2]
    name = f + ".sh"
    outFile = open(name, "w")
    outFile.write("#!/bin/bash\n\n")
    outFile.write("\n".join(commands))
    outFile.close()

    os.system("chmod +x " + name)
    os.system("./" + name)
    print("sample " + f + " is running")
