# This script is designed to categorize the candidates from a certain
# candidate type, e.g. IIIb, based on the number of common vertices shared 
# by NuPACK result and desired topology.
#
#
#
import os.path
import sys
import os

targetType = 'IIIb'

nupackResultDir = sys.argv[1] # S.J. to take command line input
prefix = "Top"
comparisonInpf = sys.argv[2]

# change nupackResultDir and comparisonInpf variables based on the target topology and corresponding folders

pool = []

with open(comparisonInpf) as pointer:
     for line in pointer:
         #
         if len(line) > 1:
            if "Candidate" in line:
               #
               if targetType in line:
                  #
                  label = str( line.split()[3] )
                  pool.append( label )



for i in range(len(pool)):
# Run TreeGraph code for NuPACK and target topology 
  # the bpseq files are available in {nupackResultDir}
  #
  if not os.path.isfile(nupackResultDir+prefix+pool[i]+".bpseq"):
     print "Error:",prefix+pool[i]+".bpseq not exist."
     sys.exit()

  if not os.path.isfile(nupackResultDir+prefix+pool[i]+".nupack.bpseq"):
     print "Error:",prefix+pool[i]+".nupack.bpseq not exist."
     sys.exit()


  #os.system( "python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py "+nupackResultDir+prefix+pool[i]+".bpseq "+ " > "+prefix+pool[i]+".tg_log"    )
  os.system( "python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/modified-treeGraph/treeGraphs.py "+nupackResultDir+prefix+pool[i]+".bpseq "+ " > "+prefix+pool[i]+".tg_log"    )


  #os.system( "python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py "+nupackResultDir+prefix+pool[i]+".nupack.bpseq "+ " > "+prefix+pool[i]+".nupack.tg_log"    )
  os.system( "python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/modified-treeGraph/treeGraphs.py "+nupackResultDir+prefix+pool[i]+".nupack.bpseq "+ " > "+prefix+pool[i]+".nupack.tg_log"    )

  #  python analyzer.py {current_topology}.log {target_topology}.log
  #
  #os.system( "python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script/analyzer.py "+prefix+pool[i]+".nupack.tg_log "+ prefix+pool[i]+".tg_log "+" > "+prefix+pool[i]+".summary1" )
  os.system( "python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/analyzer.py "+prefix+pool[i]+".nupack.tg_log "+ prefix+pool[i]+".tg_log "+" > "+prefix+pool[i]+".summary1" )

  os.system("rm -rf "+prefix+pool[i]+".nupack.tg_log "+prefix+pool[i]+".tg_log ")

  # count number of correlated vertices 
  # count "':"
  # count "->"
  nv = 0
  with open(prefix+pool[i]+".summary1") as pointer:
       for line in pointer:
           #
           if len(line) > 1:
              nv = nv + line.count("':") 
              nv = nv + line.count("->")

  #
  os.system("rm -rf "+prefix+pool[i]+".summary1")
  print prefix+pool[i],"nv=",nv



# Make folders and sub-folders 




