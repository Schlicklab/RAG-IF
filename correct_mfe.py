# This script is designed to mutate the sequence whose 
# RNAfold-mfe structure has the correct topology but 
# partially wrong correlation in vertices.  
#
# This part will be interfaced to the 'master(-gaif)' script   
# in the S3 fail branch.  
#

# issue: 


import os
import os.path
import sys
from subprocess import check_output

# get filenames 
for fname in os.listdir('.'):
    if fname.endswith('.seq'):
       #print fname
       prefix = fname.split(".")[0]
       break 


inpf1 = prefix+".tg_log"
inpf2 = prefix+".mfe.tg_log"
inpf3 = prefix+".bpseq"

#prefix = inpf3.split(".")[0]



analyze_path = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S1-script/"
s2_path = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/analyzer/S2-script/"
me_path = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/mutationEngine/"
sub_path = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/Cases/subroutines/" 
RNAinverseDir="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/GAIF/"


def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])



def getSS(ex_seq): # based on NUPACK

       f1 = open("tmpNupack.in","w")
       f1.write(ex_seq)
       f1.close()

       os.system("/opt/nupack3.2.2/bin/mfe -material rna tmpNupack")
       os.system("rm tmpNupack.in")

       target = ""
       with open("tmpNupack.mfe") as pointer:
            for line in pointer:
                if "((" in line:
                   target = line
                   break
                if ".." in line:
                   target = line
                   break
                if "))" in line:
                   target = line
                   break

       os.system("rm -rf tmpNupack.mfe")
       return target



def getSS_RNAfold( ex_seq ):

    f1 = open("tmpRNAfold.in",'w')
    f1.write(ex_seq)
    f1.close()

    os.system("/opt/viennaRNA2.4.8/bin/RNAfold -p -d2 --noLP < tmpRNAfold.in > tmpRNAfold.out" )
    os.system("rm -rf tmpRNAfold.in rna.ps dot.ps")

    # we don't need this fasta file
    #f1 = open("tmpRNAfold-mfe.fa","w")
    #f1.write(">test\n")
    #f1.write(ex_seq)
    #f1.close()

    # we only take the mfe structure
    counter = 0
    with open("tmpRNAfold.out") as pointer:
         for line in pointer:
             counter = counter + 1
             if counter == 2:
                tmpRNA_mfe_ss = line.split()[0]
    os.system("rm -rf tmpRNAfold.out")

    return tmpRNA_mfe_ss


def getTopo_RNAfold( ex_seq, flag ):

    fd2seqDir = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/"
    TGpath="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/"

    if flag == 1:
       ss_RNAfold = getSS_RNAfold( ex_seq )
    if flag == 2:
       ss_RNAfold = getSS( ex_seq ) # NUPACK 

    # prepare fasta and dotbracket files
    f1 = open("tmpRNAfold-fasta-1.fa",'w')
    f1.write(">test1\n")
    f1.write(ex_seq)
    f1.close()

    f1 = open("tmpRNAfold-fasta-1.dotbracket",'w')
    f1.write(ss_RNAfold)
    f1.close()

    # get bpseq file
    os.system('rm -rf tmpRNAfold-fasta-1.bpseq')
    os.system("python "+fd2seqDir+"dotfa2bpseq.py tmpRNAfold-fasta-1.fa tmpRNAfold-fasta-1.dotbracket"+">/dev/null ")
    os.system("rm -rf tmpRNAfold-fasta-1.fa tmpRNAfold-fasta-1.dotbracket")

    # run treeGraph
    rnafold_mfe_top = ""
    os.system("tail -n +2 tmpRNAfold-fasta-1.bpseq > tmp1 ")
    os.system("mv tmp1 tmpRNAfold-fasta-1.bpseq")
    os.system("python "+TGpath+"treeGraphs.py tmpRNAfold-fasta-1.bpseq > tmpRNAfold-fasta-1.tg_log")
    with open("tmpRNAfold-fasta-1.tg_log") as pointer:
         for line in pointer:
             if "Graph ID" in line:
                 rnafold_mfe_top = line.split()[-1].strip()

    os.system("rm -rf tmpRNAfold-fasta-1.bpseq")#  tmpRNAfold-fasta-1.tg_log ")

    #print "current topology is:", rnafold_mfe_top

    return rnafold_mfe_top

#
def getTopo_forFiles( ex_seq, flag, prefix ):

    fd2seqDir = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/"
    TGpath="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/"

    if flag == 1:
       ss_RNAfold = getSS_RNAfold( ex_seq )
    if flag == 2:
       ss_RNAfold = getSS( ex_seq ) # NUPACK 

    # prepare fasta and dotbracket files
    f1 = open("tmpRNAfold-fasta-1.fa",'w')
    f1.write(">test1\n")
    f1.write(ex_seq)
    f1.close()

    f1 = open("tmpRNAfold-fasta-1.dotbracket",'w')
    f1.write(ss_RNAfold)
    f1.close()

    # get bpseq file
    os.system('rm -rf tmpRNAfold-fasta-1.bpseq')
    os.system("python "+fd2seqDir+"dotfa2bpseq.py tmpRNAfold-fasta-1.fa tmpRNAfold-fasta-1.dotbracket"+">/dev/null ")
    os.system("rm -rf tmpRNAfold-fasta-1.fa tmpRNAfold-fasta-1.dotbracket")

    # run treeGraph
    rnafold_mfe_top = ""
    os.system("cp tmpRNAfold-fasta-1.bpseq tmpRNAfold-fasta-1.bpseq.bak")
    os.system("tail -n +2 tmpRNAfold-fasta-1.bpseq > tmp1 ")
    os.system("mv tmp1 tmpRNAfold-fasta-1.bpseq")
    os.system("python "+TGpath+"treeGraphs.py tmpRNAfold-fasta-1.bpseq > tmpRNAfold-fasta-1.tg_log")

    #with open("tmpRNAfold-fasta-1.tg_log") as pointer:
    #     for line in pointer:
    #         if "Graph ID" in line:
    #             rnafold_mfe_top = line.split()[-1].strip()

    if flag == 1:
       os.system("rm -rf tmpRNAfold-fasta-1.bpseq")
       os.system("mv tmpRNAfold-fasta-1.bpseq.bak "+prefix+".rnafold.bpseq")
       os.system("mv tmpRNAfold-fasta-1.tg_log "+prefix+".rnafold.tg_log")

    if flag == 2:
       os.system("rm -rf tmpRNAfold-fasta-1.bpseq") #  tmpRNAfold-fasta-1.tg_log ")
       os.system("mv tmpRNAfold-fasta-1.tg_log "+prefix+".nupack.tg_log")

    return 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# main 

# get the designed seq.
f1 = open(prefix+".seq",'r')
seq_org = f1.readline().strip()
#print seq_org
#sys.exit()


# correlate design_tg and mfe_tg
os.system("python "+analyze_path+"analyzer.py "+ inpf2 + " " + inpf1  + "> design_mfe.s1_summary" )

os.system("python "+s2_path+"s2_script.py "+"design_mfe.s1_summary "+inpf2+" "+inpf1+"> dm-ifmutate.log")
# 'yes' - to be mutated 
# 'no' - not to be mutated 

# run Mutation Engine
os.system("python "+me_path+"mutEngine.py  "+inpf3+"  dm-ifmutate.log")

#sys.exit()




ifsufficient = 0
nstepheaven = -2
while ifsufficient == 0:
   
   # prepare the '.gaif_conf' file
   nstepheaven = nstepheaven + 2
   os.system("rm -rf .gaif_conf")
   f1 = open(".gaif_conf","w")
   f1.write("nstepheaven = "+str(nstepheaven)+"\n")
   f1.close()
    

   # Run GAIF 
   os.system("python "+RNAinverseDir+"gaif.py "+prefix+".rnaInverseInp "+prefix+".seq "   )
   #sys.exit()
   #
   # python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/GAIF/gaif.py Top104.rnaInverseInp   Top104.seq 
   #

   # check whether the sequences in 'heaven.txt' are sufficient or not 
   if os.path.isfile("heaven.txt"):
      nseq = wc("heaven.txt")
      print nseq,"sequences are in heaven.txt"
      if nseq > 100:
         ifsufficient = 1
      else:
         continue
   else: # if heaven.txt does not even exist! 
      continue 

   # After this, we may have a lot of sequences in the 'heaven.txt' file, we need to find just a few 
   # successful sequences that have all vertices correlated with the designed one 
   # we need to check the following:
   # 1)  topology type (RNAfold)
   # 2)  whether all vertices are correlated (RNAfold and design) 
   # 3)* whether this has been already a successful one which can pass the verification of NuPACK [optional]  
   #
   # Current optimal way is to find the seuquence with minimal mutation to replace the originally designed    
   # sequence.

   seqs = []
   with open("heaven.txt") as f: # read in sequences generated by GAIF
        for line in f:
            if len(line) > 2:
               seq = line.split()[0].strip()
               seqs.append(seq)

   target_topo = ""
   with open(inpf1) as f: # get the target topology 
        for line in f:
            if "Graph ID" in line:
                target_topo = line.split()[-1].strip()
                #print target_topo
                #sys.exit()

   min_mut = 99
   min_mut_seq = ''
   for i in range(len(seqs)):
       ex_seq = seqs[i]
       current_topo = getTopo_RNAfold( ex_seq, 1 ) # "tmpRNAfold-fasta-1.tg_log" is left
       #print "debug",current_topo , target_topo
       if current_topo == target_topo:
          # correlate vertices 
          os.system("python "+analyze_path+"analyzer.py "+"tmpRNAfold-fasta-1.tg_log "+inpf1+ ">design_mfe.s2_summary")
          with open("design_mfe.s2_summary") as f:
              if "All vertices have been analyzed" in f.read():
                 # find the sequence with minimal mutations 
                 count = 0
                 for x,y in zip(seq_org, ex_seq):
                     if x != y:
                        count = count + 1
                 #print count
                 if count < min_mut:
                    min_mut = count
                    min_mut_seq = ex_seq
                    #print min_mut 



                 # check NuPACK topo. -- this does not work for 8_12, but might work for others 
                 #nupack_topo = getTopo_RNAfold( ex_seq, 2 )
                 #print nupack_topo

       #
       os.system("rm -rf tmpRNAfold-fasta-1.tg_log design_mfe.s2_summary ")
       #sys.exit()
          # 
   # - - - END OF FOR LOOP - - - - - -

   if min_mut == 99:
      print "No good sequences.\nStop!"
      f4 = open("STOP_SIGNAL","w")
      f4.write("stop at S3 fail branch")
      f4.close()
      sys.exit()

       #
   #    
   # dump the current sequence as a file
   print "find minimal mutations:", min_mut
   f2 = open("mut-"+prefix+".seq","w")
   f2.write(min_mut_seq) # +" "+str(min_mut))
   f2.close()
 
   # in the end, "mut-[prefix].seq" will be generated

   # rename and clean-up
   os.system("mv "+prefix+".rnaInverseInp "+ prefix+".rnaInverseInp.s3_fail ")


   #backup necessary files
   os.system("mv "+prefix+".nupack.tg_log "+prefix+".nupack.tg_log.bak ")
   os.system("mv "+prefix+".rnafold.tg_log "+prefix+".rnafold.tg_log.bak ")
   os.system("mv "+prefix+".rnafold.bpseq "+prefix+".rnafold.bpseq.bak ")

   # obtain the above three files 
   # 1. [prefix].rnafold.seq <- 
   #
   getTopo_forFiles( min_mut_seq, 1, prefix ) # rnafold
   getTopo_forFiles( min_mut_seq, 2, prefix ) # nupack

   #

















