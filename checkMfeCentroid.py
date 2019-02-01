# This script is designed to decide to use mfe result or centroid result 
# If both have the same graph ID, we prefer the mfe result 

# The chosen result is renamed as "rnafold"
import os.path
import os
import sys 


#inpstr = "Top49"
inpstr = sys.argv[1]


gID1 = ''
gID2 = ''
gID3 = ''

if os.path.isfile(inpstr+".tg_log"):
   for line in open(inpstr+".tg_log"):  
       if "Graph ID" in line:
           #print line.strip()
           gID1 = line.split()[-1]

#print gID1

if os.path.isfile(inpstr+".mfe.tg_log"):
   for line in open(inpstr+".mfe.tg_log"):  
       if "Graph ID" in line:
           #print line.strip()
           gID2 = line.split()[-1]

#print gID2

if os.path.isfile(inpstr+".centroid.tg_log"):
   for line in open(inpstr+".centroid.tg_log"):  
       if "Graph ID" in line:
           #print line.strip()
           gID3 = line.split()[-1]

#print gID3


mfePass = 'n'
centroidPass = 'n'


#/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script

if gID1 == gID2:
   #
   os.system("python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script/analyzer.py "+ inpstr +".mfe.tg_log "+inpstr +".tg_log >"+inpstr+".s1_summary1.a")
   #
   if os.path.isfile(inpstr+".s1_summary1.a"):
      for line in open( inpstr+".s1_summary1.a"):
          if "All vertices have been analyzed" in line:
             mfePass = 'y'
else:
    if gID1 == gID3:
       #
       os.system("python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script/analyzer.py "+ inpstr+".centroid.tg_log "+ inpstr + ".tg_log >"+ inpstr+".s1_summary1.b")
       #
       if os.path.isfile(inpstr+".s1_summary1.b"):
          for line in open(inpstr+".s1_summary1.b"):
              if "All vertices have been analyzed" in line:
                 centroidPass = 'y'

#
if mfePass == 'y':
   os.system("cp "+inpstr+".mfe.tg_log "+inpstr+".rnafold.tg_log")
   os.system("cp "+inpstr+".mfe.bpseq  "+inpstr+".rnafold.bpseq")
else:
   if centroidPass == 'y':
      os.system("cp "+inpstr+".centroid.tg_log "+inpstr+".rnafold.tg_log")
      os.system("cp "+inpstr+".centroid.bpseq " +inpstr+".rnafold.bpseq")


# clean up
os.system("rm -rf "+inpstr+".s1_summary1.a")
os.system("rm -rf "+inpstr+".s1_summary1.b")



