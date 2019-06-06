# This script is designed as a driver script that directly reads in
# (A) and (B), calls several other scripts  
#
# Necessary Input files:
#
# (A) 1. Top###.bpseq 
#     2. Top###.nupack.bpseq
#
# (B) 1. Top###.seq 
#     2. Top###.rnafold
#     3. Top###.rnafold.forna
#
# Input string: Top### as the prefix     
#
#

import os.path
import sys
import os 


# # # #
inpstr = sys.argv[1] #"Top94"
# # # #

#subroutinDir="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/Cases/subroutines/"
subroutinDir="/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/"

# added on 2018-09-19
s3fail_correct_mfe_script = subroutinDir+"correct_mfe.py"


if (not os.path.isfile(inpstr+".bpseq")) or\
   (not os.path.isfile(inpstr+".nupack.bpseq")) or\
   (not os.path.isfile(inpstr+".seq")) or\
   (not os.path.isfile(inpstr+".rnafold")) or\
   (not os.path.isfile(inpstr+".rnafold.forna")):
   #
   pass
   print "Necessary files missing.\nAbort."
   sys.exit()

# S1 Generate BPSEQ file for RNAfold result
os.system("python "+subroutinDir+"rnafoldforna2bpseq.py "+inpstr+".rnafold.forna"+" >/dev/null")

if os.path.isfile(inpstr+".mfe.bpseq") and os.path.isfile(inpstr+".centroid.bpseq"):
   print "S1 done."
else:
   print "S1 failed."
   sys.exit()


# S2 Run TreeGraph code
os.system(subroutinDir+"run-TreeGraph.sh "+ inpstr)
if os.path.isfile(inpstr+".tg_log") and os.path.isfile(inpstr+".nupack.tg_log") and\
   os.path.isfile(inpstr+".mfe.tg_log") and os.path.isfile(inpstr+".centroid.tg_log"):
   #
   print "S2 done."
else:
   print "S2 failed."
   sys.exit()

# S3 For RNAfold result, decide whether to use [mfe] result or [centroid] result
# This is done by checking whether all vertices have been correlated to the designed secondary structure  
os.system("python "+subroutinDir+"checkMfeCentroid_Swati.py "+inpstr)
if os.path.isfile(inpstr+".rnafold.tg_log"):
   print "S3 done."
else:
   print "S3 failed. trying to solve the problem... "
   
   # added 09-19-2018  
   os.system("python "+s3fail_correct_mfe_script)
   # check 'STOP_SIGNAL'
   if os.path.isfile('STOP_SIGNAL'):
      print "S3 failed." 
      sys.exit()


   #
   # we need to remove/rename several files in order to avoid conflicts 
   # 
   # new files we have introduced:
   # 1. heaven.txt - [rename]
   os.system("mv heaven.txt heaven.txt.s3_fail")

   # 2. design_mfe.s1_summary

   # 3. dm-ifmutate.log - [rename] 

   # 4. TopXXX.rnaInverseInp - [rename]
        # renamed already as 'TopXXX.rnaInverseInp.s3_fail' in the above script  

   # 5. mut-TopXXX.seq - [rename as 'TopXXX.seq']
   #                     Do not forget to backup the original TopXXX.seq file!!!
   #                     In the end, the backup file should be renamed back to 'TopXXX.seq'

   os.system("mv "+inpstr+".seq "+  inpstr+".seq.bak "    )
   os.system("mv "+"mut-"+inpstr+".seq " + inpstr+".seq " )

   # 6. '.gaif_conf' - [remove!]
   os.system("rm -rf .gaif_conf")

   #
   # files that will be used in the following steps
   # 1. TopXXX.nupack.tg_log <- bpseq <- seq   -- S4, S5 
   # 2. TopXXX.rnafold.tg_log <- bpseq <- seq  -- S4, S5 
   # 3. TopXXX.rnafold.bpseq <- seq            -- S6 
   #
   # above three files have been generated by the script above   
   #
   # 4. TopXXX.seq                             -- S7,8 - for GAIF
   # 5. TopXXX.tg_log                          -- S7,8 - only use the 'Graph ID' line 
   #
   #

   
   #sys.exit()




# S4 Pre-mutation analysis - NuPACK vs RNAfold 
#os.system("python "+"/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script/analyzer.py "+ inpstr+".nupack.tg_log " + inpstr+".rnafold.tg_log > "+inpstr+".s1_summary2 " )
os.system("python "+"/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/analyzer.py "+ inpstr+".nupack.tg_log " + inpstr+".rnafold.tg_log > "+inpstr+".s1_summary2 " )
if os.path.isfile(inpstr+".s1_summary2"):
   print "S4 done."
else:
   print "S4 failed."
   sys.exit()

# S5 Mutation targeting with S2 script -> has been corrected - 08/21/2018   
#os.system("python "+"/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S2-script/s2_script.py "+ inpstr+".s1_summary2 "+ inpstr +".nupack.tg_log "  + inpstr+".rnafold.tg_log > "+inpstr+"-ifmutate.log " )
os.system("python "+"/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/s2_script.py "+ inpstr+".s1_summary2 "+ inpstr +".nupack.tg_log "  + inpstr+".rnafold.tg_log > "+inpstr+"-ifmutate.log " )
if os.path.isfile(inpstr+"-ifmutate.log"):
   print "S5 done."
else:
   print "S5 failed."
   sys.exit()


# S6 Call Mutation Engine 
os.system("tail -n +2 "+inpstr+".rnafold.bpseq > tmp1")
os.system("mv tmp1 " + inpstr+".rnafold.bpseq")

#os.system("python "+"/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/mutationEngine/mutEngine.py "+ inpstr+".rnafold.bpseq "+ inpstr+"-ifmutate.log "+">/dev/null"  )
os.system("python "+"/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/mutEngine.py "+ inpstr+".rnafold.bpseq "+ inpstr+"-ifmutate.log "+">/dev/null"  )
if os.path.isfile(inpstr+".rnafold.rnaInverseInp"):
   print "S6 done."
else:
   print "S6 failed."
   sys.exit()

# S7 + S8
os.system("python "+subroutinDir+"runRNAinverse-gaif.py "+inpstr+".rnafold.rnaInverseInp "+inpstr+".seq ") # modified 


# If S3-fail branch was walked, the TopXXX.seq should be recovered for later analysis 
# 2018-09-19

if os.path.isfile(inpstr+".seq.bak"):
   os.system("mv "+inpstr+".seq.bak "+inpstr+".seq " ) 





