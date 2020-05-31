# This script is designed to order the final sequences in:
# 1. heaven.txt
# 2. *.survivors   
# according to the number of point mutations 
#
# also it should print out the mutation pattern 
#
# to be added: automatically remove redundant point mutation(s) and keep only
#              necessary point mutations
#
# * 2018-09-05
# issue: deleting useless mutations one by one has no problem, but sometimes if they are deleted  
#        together, then the secondary structure will be changed. 
#        This implies that there exists interactions between these mutations.  
#        --> if these mutations form a set, then we need to check its all possible subsets and find the largest   
#            set that can keep the secondary structure 
# another issue: above remedy might give a set which is not minimal !  
#
# * 2018-09-06
# issue: when analyzing the results, not all sequences in the '.survivors' are included   
#        --> comment out 
#
#
# issue: try to remove pairs introduced by two mutations, but also test whether the secondary structure 
#        stays the same or not     
#       
#
#
# * 2018-09-11
# issue 1: check inclusion among the unique mutation sets. 
#          If found, only keep the shorter mutation set.    
#
# issue 2: we did not check the topology by RNAfold for those unique mutation sets.
#
#
# * 2018-09-13
# issue: inclusion check not complete 
#
# * 2018-09-20
# issue: the "target_topo" should be specified via the command line with 'sys.argv' 
#
#

import os
import sys
from itertools import chain, combinations

# -----------------
def getSS_RNAfold( ex_seq ):
    
    f1 = open("tmpRNAfold.in",'w')
    f1.write(ex_seq)
    f1.close()

    os.system("/opt/viennaRNA2.3.5/bin/RNAfold -p -d2 --noLP < tmpRNAfold.in > tmpRNAfold.out" )
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


# ------------------- S.J. 03/19/2019 - function to get the NUPACK topology
def getTopo_Nupack(ex_seq):

    f1 = open("tmpNupack.fa","w")
    f1.write(">test\n")
    f1.write(ex_seq)
    f1.close()
                    
    target = getSS(ex_seq)
    f1 = open("tmpNupack.dotbracket","w")
    f1.write(target)
    f1.close()

    os.system("python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/dotfa2bpseq.py tmpNupack.fa tmpNupack.dotbracket"+">/dev/null")
    os.system("tail -n +2 tmpNupack.bpseq > tmp1")
    os.system("mv tmp1 tmpNupack.bpseq")
    os.system("python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/modified-treeGraph/treeGraphs.py tmpNupack.bpseq > tmpNupack.tg_log")

    nupack_topo=""
    with open("tmpNupack.tg_log") as pointer:
        for line in pointer:
            if "Graph ID" in line:
                nupack_topo = line.split()[-1].strip()

    os.system("rm -rf tmpNupack.*")
    return nupack_topo


# ----------------- S.J. 03/20/2019
def check_inclusion(pattern): # to check if this pattern includes any of the unique patterns already selected

    global unique_pattern

    for pat in unique_pattern: # for every pattern in unique pattern
        
        flag = True # assume inclusion
        for mut in pat: # for every mutation in this unique pattern
            if mut not in pattern: # if this mutation is not in the pattern
                flag = False
                break

        if flag: # if all mutations in pat are in pattern, therefore this is redundant
            return flag

    # code will come here if none of the unique patterns are included in the pattern
    return False


# -----------------
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


#--------------
def inputDotBracket(array): # find pairs from dot-bracket notations 

    bracketPositions = []
    results = []
    for i, item in enumerate(array):
        if i == 0 and item == ')':
            pass
            #print("Non sense ! Exit")
            break

        if item == '(':
            bracketPositions.append(i)
        elif item ==')':
            if len(bracketPositions) > 0:
                openingPosition = bracketPositions.pop()
                results.append([openingPosition+1,i+1 ]  )
                #print(openingPosition+1, '-->', i+1)
        else:
            pass
                #print('ERROR: Not a bracket. Word is: %s.' % item)
    return results


#--------------- S.J. 03/20/2019 - changing the name of the parameters
def get_seq(seq1, res_num, seq2 ):
    
    list1 = list(seq1)
    list2 = list(seq2)
    
    #print "Replacing " + str(res_num) + " from " + list1[ res_num-1 ] + " to " + list2[ res_num-1 ]

    list1[ res_num-1 ] = list2[ res_num-1 ]
    test_seq = ''.join(list1)

    return test_seq


#--------------
def all_subsets(ss): # find all subsets except the null set   
    # ss is also a list 
    a = chain(*map(lambda x: combinations(ss, x), range(1, len(ss)+1)))
    list_a = []
    for k in a:
        list_a.append( list(k) )
    
    return list_a


#---------------
# to add topology check[getTopo_RNAfold] - 2018/09/12    
def kick_one( pattern, test_seq, seq_org, ss_current, target_topo ): 

    print "Kicking..."

    istop = 0
    to_kick = []

    for i in range(len(pattern)):
        res_num = int( pattern[i].split('-')[0] )

        new_seq = get_seq( test_seq, res_num, seq_org )
        #SJ commented new_ss = getSS( new_seq )
        new_ss = getTopo_Nupack( new_seq )
        new_topo = getTopo_RNAfold( new_seq )

        #SJ commented if new_ss == ss_current and new_topo == target_topo:# modified
        if new_ss == target_topo and new_topo == target_topo:
           to_kick.append( pattern[i] )
           out_seq = new_seq+'.'[:-1]
           break

        if i == (len(pattern)-1):
           out_seq = test_seq+'.'[:-1] 
           istop = 1

    if len(to_kick) != 0:
       pattern.remove( to_kick[0] )

    return pattern, istop, out_seq


#---------------
def getTopo_RNAfold( ex_seq ):

    #fd2seqDir = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/" 
    fd2seqDir = "/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/" 
    TGpath="/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/modified-treeGraph/"

    ss_RNAfold = getSS_RNAfold( ex_seq )

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

    os.system("rm -rf tmpRNAfold-fasta-1.bpseq  tmpRNAfold-fasta-1.tg_log ")
    #print "current topology is:", rnafold_mfe_top
    return rnafold_mfe_top

# -------------------------
# we need to add topology check for RNAfold secondary structure - 2018/09/12
def pat_analyzer( pattern, ss_current, ss_org, seq_sur, seq_org, target_topo ):

    pairs_current = inputDotBracket( ss_current ) # find pairs from dot-brackets
    pairs_org     = inputDotBracket( ss_org )

    pattern_2 = []
    info_2 = []

    useful_i = []
    
    print pattern, "a"
    #for s in range(0,len(seq_org)):
    #    if seq_org[s] != seq_sur[s]:
    #        print s+1,seq_org[s],"-",seq_sur[s] # this should be exactly the same as the min pattern being printed above

    info_1 = []
    # Analysis 1: find paired partner
    for i in range(len( pattern )):
       
        #print pattern 
        res_num = int( pattern[i].split('-')[0] )

        pair_num_org = 0 # the residue number of its paired partner
        pair_num_current = 0

        for j in range(len(pairs_org)):
            if pairs_org[j][0] == res_num:
               pair_num_org = pairs_org[j][1]
            if pairs_org[j][1] == res_num:
               pair_num_org = pairs_org[j][0]

        for j in range(len(pairs_current)):
            if pairs_current[j][0] == res_num:
               pair_num_current = pairs_current[j][1]
            if pairs_current[j][1] == res_num:
               pair_num_current = pairs_current[j][0]
        print res_num, pair_num_org, pair_num_current

        info_1.append( [res_num, pair_num_org, pair_num_current] )

        #if pair_num_org != pair_num_current: # exclude loop->loop and pair->pair(the same one)
        #   useful_i.append( i )
        #   pattern_2.append(  pattern[i] )
        #   info_2.append( [ res_num, pair_num_org, pair_num_current ] )



    # Analysis 2: try to find the mininal set of mutation 'pattern' 
    #             we need to loop over all mutations in this pattern, remove it then check  
    #             whether the secondary structure will change or not. If SS is not changed, 
    #             then safely remove this point mutation. 
    #             This should repeat until all point mutations will lead to change in SS. 
    #               
    # input: pattern, ss_current, seq_sur, seq_org

    allchange = 0
    min_pattern = pattern[:] # need to be updated
    new_seq_sur = (seq_sur + '.')[:-1] # need to be updated
    while allchange == 0:
          delList = []
          keep_num_list = []
          for i in range(len( min_pattern )):
              res_num = int( min_pattern[i].split('-')[0] )
              #print "Removing the mutation"
              test_seq = get_seq( new_seq_sur, res_num, seq_org ) # delete this point mutation
              # SJ commented - test_ss = getSS( test_seq ) # test the secondary structure afterwards
              test_ss = getTopo_Nupack(test_seq) # S.J. 03/19/2018 - NUPACK topology
              
              test_topo = getTopo_RNAfold( test_seq ) # get the RNAfold topology

              # SJ commented if test_ss == ss_current and test_topo == target_topo:
              if test_ss == target_topo and test_topo == target_topo:
                 delList.append( min_pattern[i] )# collects the mutation that can be safely deleted
              else:
                 keep_num_list.append( res_num ) # necessary point mutations

          if len(delList) == 0:
             allchange = 1
             break # break the while loop

          # update the 'new_seq_sur'
          new_seq_sur = (seq_org+'.')[:-1]
          for i in range(len(keep_num_list)):
              #print "Adding back the mutation"
              new_seq_sur = get_seq( new_seq_sur, keep_num_list[i], seq_sur)
          
          # update the 'min_pattern'
          for i in range(len(delList)):
              min_pattern.remove( delList[i] )

    # - - - - - -  - - - - - - -  - - - - - -  - - - -  -  - -

    # double check the secondary structure of 'new_seq_sur' 
    
    print new_seq_sur
    for s in range(0,len(seq_org)):
        if seq_org[s] != new_seq_sur[s]:
            print s+1,seq_org[s],"-",new_seq_sur[s] # this should be exactly the same as the min pattern being printed next
    
    # SJ commented check_ss = getSS( new_seq_sur )
    check_ss = getTopo_Nupack( new_seq_sur )
    # and topology...
    check_topo = getTopo_RNAfold( new_seq_sur )
    print check_ss, " ", check_topo

    #SJ commented if check_ss == ss_current and check_topo == target_topo:
    if check_ss == target_topo and check_topo == target_topo:
        print "remove ",str(len(pattern)-len(min_pattern)), "point mutations"
        print min_pattern #, "b"
        pattern_2 = min_pattern
        #print ""
    else:
        print "try to remove ",str(len(pattern)-len(min_pattern)), "point mutations"
        print min_pattern
        print "failed"
        
        # remedy solution -> add the smallest subset of 'deleted' into min_pattern  
        deleted = list( set(pattern) - set(min_pattern) )
        subsets_del = all_subsets( deleted )
        test_min_pattern = min_pattern[:] # added on 2018-09-22
        for i in range(len(subsets_del)):
            test_min_pattern = min_pattern + subsets_del[i]
            test_seq = (seq_org+'.')[:-1]
            for j in range(len(test_min_pattern)):
                res_num = int( test_min_pattern[j].split('-')[0] )
                test_seq = get_seq( test_seq, res_num, seq_sur )
             # SJ commented test_ss = getSS( test_seq )
            test_ss = getTopo_Nupack(test_seq)
            test_topo = getTopo_RNAfold( test_seq )
            # SJ commented if test_ss == ss_current and test_topo == target_topo:
            if test_ss == target_topo and test_topo == target_topo:
               print test_min_pattern #, "b"
               #pattern_2 = test_min_pattern # wait a minute...
               break # the for loop of 'i'
    
        # to solve another issue... 
        #    input: test_min_pattern 
        #
        #
        istop = 0
        while istop == 0:
            # kick_one has been modified with regard to topology
            # issue: 2018-09-22 
            # 'test_min_pattern' referenced before assignment 
            test_min_pattern, istop, test_seq = kick_one( test_min_pattern, test_seq, seq_org, ss_current, target_topo )

        print test_min_pattern #, 'b'
        pattern_2 = test_min_pattern


    # Analysis 3: find pairs introduced by two mutations, and test whether removing them will change  
    #             the secondary structure or not.  
    # 
    # input: info_1, pattern_2, ss_current, ... 
    #
    #
    pairs_int = []
    for i in range(len(info_1)):
        first_n = info_1[i][0]
        third_n = info_1[i][2]
        for j in range(len(info_1)):
            if j > i:
               first_m = info_1[j][0]
               third_m = info_1[j][2]
               if first_n == third_m and first_m == third_n:
                  print "found pair!"
                  pairs_int.append( [first_n, first_m] )

    num_list = []
    for i in range(len(pattern_2)):
        res_num = int( pattern_2[i].split('-')[0] )
        num_list.append( res_num )

    for i in range(len(pairs_int)):
        a = pairs_int[i][0]
        b = pairs_int[i][1]
 
        if a in num_list and b in num_list:
             print "pairs found in pattern_2"
             
             # try to remove these two mutations and test the secondary structure 
             test_seq2 = (seq_org+".")[:-1]
             for j in range(len(pattern_2)):
                 res_num = int( pattern_2[j].split('-')[0] )
                 if ( res_num != a ) and ( res_num != b ):
                    test_seq2 = get_seq( test_seq2, res_num, seq_sur )

             # SJ commented test_ss2 = getSS( test_seq2 )
             test_ss2 = getTopo_Nupack( test_seq2 )

             # get topology 
             test_topo2 = getTopo_RNAfold( test_seq2 )

             # SJ commented - if test_ss2 == ss_current and test_topo2 == target_topo:
             if test_ss2 == target_topo and test_topo2 == target_topo:
                to_remove = []
                for k in range(len(pattern_2)):
                    res_num = int( pattern_2[k].split('-')[0] )
                    if res_num == a or res_num == b:
                       to_remove.append( pattern_2[k] )
                for k in to_remove:
                    pattern_2.remove( k )
                #print pattern_2, 'b'
             # update pattern_2
    print pattern_2, 'b'

    return pattern_2


###########
#  MAIN START #
###########

if len(sys.argv) != 2:
    print "Please double check input command line."
    print "usage: python [this_script.py] [target_topo]"
    sys.exit()

target_topo = sys.argv[1]

for file in os.listdir("."):
    if file.endswith(".seq"):
        inpf1 = file

for file in os.listdir("."):
    if file.endswith(".survivors"):
        inpf2 = file

# read in original seq.
seq_orig = ""
with open(inpf1) as f:
     seq_org = f.readline().strip()
ss_org = getSS(seq_org)

# order ".survivors" first
seq_sur = []
dif_sur = []

with open(inpf2) as f: # read in sequences in ".survivors" and calculate the number of different residues  
     for line in f:
         if len(line) > 2:
            flag = line.split()[-1] 
            if flag == "C":
               continue
         
            temp_seq = line.split()[0].strip()
            count = 0
            temp_pattern = []
            seq_sur.append(temp_seq)
            
            for x, y in zip( seq_org, line.split()[0].strip() ):
                if x != y:
                    count = count + 1

            dif_sur.append( count )
#print dif_sur

ord_seq_sur = [x for _,x in sorted(zip(dif_sur,seq_sur))]
ord_dif_sur = sorted( dif_sur ) # number of mutations 

#print ord_seq_sur
#print sorted(dif_sur)

#for i in range(len(ord_seq_sur)):
    #print ord_seq_sur[i], ord_dif_sur[i]
    #pass

# print out mutation pattern
#print "Print mutation patterns:"
#for i in range(len(ord_seq_sur)):
#    print ord_pat_sur[i]
  #if ord_dif_sur[i] <= max_mut: # comment out 2018-09-06
  #pattern = []
  #seq_sur = ord_seq_sur[i]
  #for j in range(len(seq_org)):
  #      if seq_org[j] != seq_sur[j]:
  #         pattern.append( str(j+1)+"-"+seq_sur[j] )

#print " "

#print len(ord_seq_sur)


# S.J. 03/20/2019 - new code to check patterns and then take unique sequences - based on old code
unique_pattern = []
for i in range(len(ord_seq_sur)):

    cur_pat = []
    cur_seq = ord_seq_sur[i]
    for j in range(len(seq_org)):
        if seq_org[j] != cur_seq[j]:
            cur_pat.append( str(j+1)+"-"+cur_seq[j] ) # calculating the pattern
    
    flag = check_inclusion(cur_pat) # checking if this pattern includes any of the other unique patterns already selected
    if flag: # this includes one of the already selected unique patterns, so this is redundant, so skip this sequence
        continue

    print "\nOptimizing Sequence " + str(i+1)
    cur_ss = getSS(cur_seq)
    # come here if this sequence needs to be analyzed
    opt_pat = []
    opt_pat  =  pat_analyzer(cur_pat, cur_ss, ss_org, cur_seq, seq_org, target_topo)
    opt_pat = sorted(opt_pat)
    unique_pattern.append(opt_pat) # add it to the unique patterns, it will be different from others already selected, otherwise this sequnce would have been skipped anyways

print "unique patterns: --b"
for k in unique_pattern:
    print k, 'b'


# count the unique secondary structure (NuPACK secondary structures)
#if 1==1: # switch

#   ss_list = []
#   for i in range(len(ord_seq_sur)):
     #if ord_dif_sur[i] <= max_mut: # comment out 2018-09-06
#       ex_seq = ord_seq_sur[i]

#       f1 = open("tmpNupack.in","w")
#       f1.write(ex_seq)
#       f1.close()

#       os.system("/opt/nupack3.2.2/bin/mfe -material rna tmpNupack")
#       os.system("rm tmpNupack.in")

#       target = ""
#       with open("tmpNupack.mfe") as pointer:
#            for line in pointer:
#                if "((" in line:
#                   target = line
#                   break
#                if ".." in line:
#                   target = line
#                   break
#                if "))" in line:
#                   target = line
#                   break
                
#       os.system("rm -rf tmpNupack.mfe")
       #print target
#       if target not in ss_list:
#          ss_list.append( target )

   # check each s.s. type in ss_list, put sequences of the same type together
#   skipme = []
#   for i in range(len(ss_list)):
#       print "secondary structure: -- b"
#       ex_ss = ss_list[i]
#       print ex_ss.strip()
#       seq_type = [] # collect seq.
#       num_list = [] # collect number of mutations
#       j_list = [] # collect the index of seq. that belong to a specific secondary structure
#       for j in range(len(ord_seq_sur)):
         #if ord_dif_sur[j] <= max_mut:  # comment out 2018-09-06  
#           if j in skipme:
#              continue # in order to save time
#           this_ss = getSS( ord_seq_sur[j] )
#           if this_ss == ex_ss:
#             seq_type.append( ord_seq_sur[j] )
#              skipme.append( j )
#              num_list.append( ord_dif_sur[j] )
#              j_list.append( j )
#           pass
#       print len(seq_type), num_list              #, j_list

       # do analysis
#       unique_pattern = []
#       for k in range(len(j_list)):

           # get mutation pattern 
#           pattern = []
#           seq_sur = ord_seq_sur[  j_list[k]  ]
#           for l in range(len(seq_sur)):
#               if seq_sur[l] != seq_org[l]: # compare current seq. and original seq.
#                  pattern.append( str(l+1)+"-"+seq_sur[l] )
           #print pattern

           # analyze this pattern - remove useless patterns  
#           pattern_2  =  pat_analyzer( pattern, ex_ss, ss_org, seq_sur, seq_org, target_topo )
#           pattern_2 = sorted( pattern_2 )
#           if pattern_2 not in unique_pattern:
#              unique_pattern.append( pattern_2 )

       # remove 'inclusion' in unique patterns
       # - added on 2018-09-12
#       toRemove = [ ]
#       if len(unique_pattern) > 1:
#          for j in range(len(unique_pattern)):
#              for k in range(len(unique_pattern)):
#                  if k > j:
#                     p1 = unique_pattern[j]
#                     p2 = unique_pattern[k]
#                     ifall = 1 # assume inclusion first
#                     for l in p1:
#                        if l not in p2:
#                           ifall = 0
#                     if ifall == 1:
#                        if p2 not in toRemove:
#                           toRemove.append( p2 )
                     
#                     ifall_b = 1 # added 2018-09-13
#                     for l in p2:
#                         if l not in p1:
#                            ifall_b = 0
#                     if ifall_b == 1:
#                        if p1 not in toRemove:
#                           toRemove.append( p1 )
                        
#          if len(toRemove) != 0:
#             for j in range(len(toRemove)): # problematic!!!
#                 unique_pattern.remove( toRemove[j] )

#       print "unique patterns: --b"
#       for k in unique_pattern:
#           print k, 'b'





