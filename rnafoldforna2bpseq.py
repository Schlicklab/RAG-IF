# This script is designed to convert the RNAfold.forna to BPSEQ
# Y. Tao - 2018-07-09


#import StringIO
import os 
import sys

PathDir="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/"



#inpf1 = "Top49.rnafold.forna"
inpf1 = sys.argv[1]




prefix = inpf1.split(".")[0]


lcount = 0

with open(inpf1) as pointer:
     for line in pointer:
         #
         lcount = lcount + 1
         if lcount == 2: 
            fasta_info = line.strip()
         if lcount == 3:
            mfe_ss = line.strip()
         if lcount == 6:
            centroid_ss = line.strip()



#f = StringIO.StringIO()
#f.write(">"+prefix+" temporary fasta file\n")
#f.write(fasta_info)
#f.close()

fa_file_1 = prefix+"-tmp-mfe.fa"
fa_file_2 = prefix+"-tmp-centroid.fa"
dot_file_1 =prefix+"-tmp-mfe.dotbracket"
dot_file_2 =prefix+"-tmp-centroid.dotbracket" 




f = open(fa_file_1,"w")
f.write("> test\n")
f.write(fasta_info+"\n")
f.close()

f = open(fa_file_2,"w")
f.write("> test\n")
f.write(fasta_info+"\n")
f.close()

f = open(dot_file_1,"w")
f.write(mfe_ss+"\n")
f.close()

f = open(dot_file_2,"w")
f.write(centroid_ss+"\n")
f.close()



# python dotfa2bpseq.py {INPUT}.fa {INPUT}.dotbracket


os.system("python "+PathDir+"dotfa2bpseq.py "+fa_file_1+" "+dot_file_1)
os.system("rm -rf "+fa_file_1)
os.system("rm -rf "+dot_file_1)

os.system("python "+PathDir+"dotfa2bpseq.py "+fa_file_2+" "+dot_file_2)
os.system("rm -rf "+fa_file_2)
os.system("rm -rf "+dot_file_2)


os.system("mv "+prefix+"-tmp-mfe.bpseq "+prefix+".mfe.bpseq")
os.system("mv "+prefix+"-tmp-centroid.bpseq "+prefix+".centroid.bpseq")











