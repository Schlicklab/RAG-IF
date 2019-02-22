# This script is designed to collect necessary input files automatically 
# given a prefix name 
#
#
#
import os.path
import sys
import os

prefix = sys.argv[1] #"Top236"

#pathDir1 = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/Swati_Test/BPSEQ/"
#pathDir2 = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/Swati_Test/SEQ/"
# change these two paths according to the target and your folder
pathDir1 = sys.argv[2] + "/BPSEQ/"
pathDir2 = sys.argv[2] + "/SEQ/"

if os.path.isfile(pathDir1+prefix+".bpseq"):
   os.system("cp "+pathDir1+prefix+".bpseq "+ " .")
else:
   print "Error: file missing.\nAbort."
   sys.exit()

if os.path.isfile(pathDir1+prefix+".nupack.bpseq"):
   os.system("cp "+pathDir1+prefix+".nupack.bpseq "+ " .")
else:
   print "Error: file missing.\nAbort."
   sys.exit()

if os.path.isfile(pathDir2+prefix+".seq"):
   os.system("cp "+pathDir2+prefix+".seq "+ " .")
else:
   print "Error: file missing.\nAbort."
   sys.exit()


if os.path.isfile(pathDir2+prefix+".rnafold"):
   os.system("cp "+pathDir2+prefix+".rnafold "+ " .")
else:
   print "Error: file missing.\nAbort."
   sys.exit()


if os.path.isfile(pathDir2+prefix+".forna"):
   os.system("cp "+pathDir2+prefix+".forna "+ prefix+".rnafold.forna")
else:
   print "Error: file missing.\nAbort."
   sys.exit()

print "Necessary files have been copied."









