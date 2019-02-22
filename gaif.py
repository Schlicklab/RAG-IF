# -*- coding: UTF-8 -*-

# Genetic Alogorithm Inverse Folding 
#
# Y. Tao - 2018-08-13
#
# usage: 
# python gaif.py XXXXX.rnaInverseInp [YYYYY.template]
#

# 08/20/2018 - add functionality of using template as the initial guess 


import random
import os
import os.path
import sys
import time
#import numpy 
#import multiprocessing as mp
from functools import partial
import multiprocessing



class Chrom:
    chrom = []
    fitness = 0.0

    def __init__(self, ngene):
        self.chrom = ['N']*ngene 

    def assign(self):
        lib = ['A','U','C','G']
        ngene = len(self.chrom)
        for i in range(ngene):
            self.chrom[i] = random.choice( lib )

#---------------

def eachFit(chrom, seq, ss_target, Nindex,k):
    #
    #print chrom
    fullseq = list(seq)
    for i in range(len(Nindex)):
        fullseq[ Nindex[i] ] = chrom[i]
    fullseq = ''.join(fullseq)
    #print fullseq

    jobID = str(random.randint(10000,99999))
    jobID = jobID+str(time.time()).split('.')[1]
    list1= ['a','b','c','d','e','f','g','h','i','j']
    list2= [1,2,3,4,5,6,7,8,9,0]
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))

    if k == 1 : # Run RNAfold
       f1 = open("tmpRNAfold"+jobID+".in","w")
       f1.write(fullseq)
       f1.close()

       os.system("/opt/viennaRNA2.3.5/bin/RNAfold -p -d2 --noLP < tmpRNAfold"+jobID+".in > tmpRNAfold"+jobID+".out 2>/dev/null " )
       #passMe = os.path.isfile("rna.ps") and os.path.isfile("dot.ps")
       passMe = os.path.isfile("tmpRNAfold"+jobID+".out")
       if passMe:
          passMe = ( os.stat("tmpRNAfold"+jobID+".out").st_size != 0 )
       os.system("rm -rf tmpRNAfold"+jobID+".in rna.ps dot.ps")
      
       
       if passMe:
          with open("tmpRNAfold"+jobID+".out") as f:
               f.readline()
               ss_new = f.readline().split()[0]
          os.system("rm -rf tmpRNAfold"+jobID+".out")
       else:
          ss_new = '.'*len(ss_target)
          os.system("rm -rf tmpRNAfold"+jobID+".out")

          #print ss_new     
          #print ss_target

       # get fitness 
       fitness = 0
       for x,y in zip( ss_new, ss_target ):
           if x == y:
              fitness = fitness + 1
    #print fitness

    return fitness 



def fit4par(idx, pop, seq, ss, Nindex, k):
    pop[idx].fitness = eachFit(pop[idx].chrom, seq, ss, Nindex, k)
    #print pop[idx].fitness

    return pop[idx] 

def calcFit(pop,seq,ss,Nindex,k,nproc):
    # k decides which engine to use 

 if nproc < 1:
    nproc = 1
    
 if nproc > 1:
    whole = range(len(pop)) 
    pool = multiprocessing.Pool(nproc)
    partial_fit4par = partial(fit4par, pop=pop, seq=seq, ss=ss, Nindex=Nindex, k=k)
    result = pool.map( partial_fit4par, whole )
    pool.close()
    pool.join()

    #print len(result), type(result[0]),result[0].fitness
    #print result[0].chrom
    #sys.exit()

    ################################################################
    #
    # Serial version
 if nproc == 1:   
    for ic in range(len(pop)):
        pop[ic].fitness = eachFit(pop[ic].chrom, seq, ss, Nindex, k)
    return pop
    ################################################################
    # fit4par(pop[c], seq, ss, Nindex, k)
    
 return result #return pop


def findBest( pop ):
    
    best = ['i', 0]
    for ic in range(len(pop)):
        #print len(pop[ic]),len(pop)
        #sys.exit()
        if best[1] < pop[ic].fitness:
           best = [ic, pop[ic].fitness] # modified 
    return best


def meanFit( pop ):
    N = len(pop)
    sumFit = 0
    for c in pop:
        sumFit = sumFit + c.fitness # modified

    mean = sumFit/(N+0.0)
    return mean





def select( pop, nreplace ):
   
    N = len(pop)
    if nreplace >= 0.3*N:
       nreplace = int(0.15*N)

    # find the worst ones 
    worstList = []
    for i in range(0, nreplace):
        worst = 9999
        worstL = ''
        for ic in range(len(pop)): # modified
            if pop[ic].fitness < worst:
               if ic not in worstList:
                  worst = pop[ic].fitness
                  worstL = ic
        worstList.append( worstL )
    #
    #print worstList

    # find the best ones
    bestList = []
    for i in range(0, nreplace):
        best = 0
        bestL = ''
        for ic in range(len(pop)): # modified
            #print pop[ic].fitness
            if pop[ic].fitness > best:
               if ic not in bestList:
                  best = pop[ic].fitness
                  bestL = ic
        bestList.append( bestL )
    #
    #print bestList
    #print worstList
    #sys.exit()

    # replace worst chromosomes with best chromosomes
    for i in range(len(bestList)):
        #print pop[ worstList[i] ].chrom
        pop[ worstList[i] ].chrom = pop[ bestList[i] ].chrom[:]
        pop[ worstList[i] ].fitness = pop[ bestList[i] ].fitness # fitness should also be updated                    
        #print pop[ worstList[i] ].chrom
    return pop


def xover1( pop, nresi, acr, nx, waiveList):

    if nx > 0.5*nresi:
       nx = int(0.1*nresi)
    if nx == 0:
       nx = 1

    for i in range(len(pop)):
        for j in range(len(pop)):
            if i != j:
              if acr > random.random():  
                for k in range(nx):
                   #
                   acrNode = random.randrange(0,nresi)
                   if (i not in waiveList) and (j not in waiveList): 
                      tmp = pop[i].chrom[acrNode]
                      pop[i].chrom[acrNode] = pop[j].chrom[acrNode]
                      pop[j].chrom[acrNode] = tmp
                   if (i  in waiveList) and (j not in waiveList): 
                      tmp = pop[j].chrom[acrNode]
                      pop[j].chrom[acrNode] = pop[i].chrom[acrNode]
                      #pop[j].chrom[acrNode] = tmp
                   if (i not in waiveList) and (j in waiveList): 
                      tmp = pop[i].chrom[acrNode]
                      pop[i].chrom[acrNode] = pop[j].chrom[acrNode]
                      #pop[j].chrom[acrNode] = tmp
                   

    return pop


def getOrder( pop ):

    fitList = []
    cList = []

    for ic in range(len(pop)): # modified 
        cList.append( ic )
        fitList.append( pop[ic].fitness )

    X = cList
    Y = fitList
 
    Z = [x for _,x in sorted(zip(Y,X),reverse=True)]

    return Z


def getWaive( pop, nwaive):

    orderC = getOrder( pop )
    waiveList = []
    for i in range(nwaive):
        waiveList.append(  orderC[i]  )

    return waiveList


def mutIt( r ):
    #
    lib = ['A','U','C','G']
    
    lib.remove( r )

    s = random.choice( lib )
    return s 




def mutation1( pop, nresi, orderC, mut_good, mut_ave, mut_bad, waiveList):
    # H   - M   - L 
    # 20% - 50% - 30%
    # 1   - 2   - 1 
    highList = []
    midList = []
    lowList = []
    N = len(pop)

    for i in range(N):
        if (i+1) <= 0.2*N:
           highList.append( orderC[i] )
        if (i+1) <= 0.7*N and (i+1) > 0.2*N:
           midList.append( orderC[i] )
        if (i+1) > 0.7*N:
           lowList.append( orderC[i] )
    #print highList
    #print midList
    #print lowList

    nloop = 1
    for ic in range(len(pop)):
        if ic in highList:
           mut = mut_good
           nloop = 1
        if ic in midList:
           mut = mut_ave
           nloop = 2
        if ic in lowList:
           mut = mut_bad
           nloop = 1
        if ic in waiveList:
           continue
        if mut > random.random():
           # do mutation
           for loop in range(nloop):
               mutNode = random.randrange(0, nresi)
               pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] ) 

    #
    return pop

def heaven( pop, fullfit, nstepheaven, heaven_rate, heavenList ):
    # modify:
    #   give 50% chance to use a lower standard (by 2)  
    #        40%                                 by 4
    #        20%                                 by 6
    #        15%                                 by 8
    # should give more chances to use lower standards (by 4, 6 and 8) 
    #
    for ic in range(len(pop)):
        diff = fullfit - pop[ic].fitness

        not_good = False
        if diff > nstepheaven:
           if diff > (nstepheaven+2):
              if diff > (nstepheaven+4):
                 if diff > (nstepheaven+6):
                    if diff > (nstepheaven+8): 
                       not_good = True
                    else: # diff <= (nstepheaven+8)
                       if 0.15 > random.random():
                          not_good = False
                       else:
                          not_good = True 

                 else: # diff <= (nstepheaven+6)
                    if 0.2 > random.random():
                       not_good = False     
                    else:
                       not_good = True

              else: # diff <= (nstepheaven+4)
                 if 0.4 > random.random(): 
                    not_good = False
                 else:
                    not_good = True

           else: # diff <= (nstepheaven+2)
              if 0.5 > random.random():
                 not_good = False
              else:
                 not_good = True 


        if not_good:     #diff > nstepheaven:
           continue
        else:
           # must enter heavenList 
           str1 = ''.join( pop[ic].chrom )
           if [ str1, pop[ic].fitness ] not in heavenList:
              heavenList.append( [ str1, pop[ic].fitness ] )
           # but not necessarily to be mutated    
           if heaven_rate > random.random(): # mutation
              # mutate this chromosome once  
              mutNode = random.randrange(0, len(pop[ic].chrom))
              pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] )
#              mutNode = random.randrange(0, len(pop[ic].chrom))
#              pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] )



    #heavenList = [list(t) for t in set(tuple(element) for element in heavenList)] # remove duplicates 

    return pop, heavenList


def stuckMutation( pop, waiveList ):

    #
    #print "stuckMutation triggered..."

    nprotect = 2
    iprotect = 0
    for ic in waiveList:

        iprotect = iprotect + 1  # protect 'nprotect' seq. out of the waiveList from mutation 
        if iprotect <= nprotect:
           continue           
        #
        for loop in range(30): 
            mutNode = random.randrange(0, len(pop[ic].chrom) )
            pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] ) 

    return pop


def chkSeq(seq, seq_tmp, Nindex):
    
    # length 
    if len(seq) != len(seq_tmp):
       print 'Template sequence has inconsistent length.'
       sys.exit()

    # check info
    for i in range(len(seq)):
        if i not in Nindex:
           if seq[i] != seq_tmp[i]:
              print 'Template sequence has inconsistent info.'
              sys.exit()
    
    # get initial guess in correct format
    chrom = []
    for i in range(len(Nindex)):
        #
        chrom.append( seq_tmp[Nindex[i]] )

    #print chrom

    return chrom



def main(inpf,tmpf,nseq,nreplace,nwaive,niter,nproc,acr_ave,nx,mut_good,mut_ave,mut_bad,nstepheaven,heaven_rate,nprintheaven,nsurvivors,nstill_1,nstill_0):

   # set up timer 
   timer_a = time.time()
   total_time = 30*60 # 30 min is the upper limit of execution 

   # Read in info. of input file
   with open( inpf ) as f:
        ss = f.readline().strip()
        seq = f.readline().strip().upper()

   print "Length of Seq.", len(seq)

   Nindex = [ ] # collect the index of unknown residues 
   for r in range(len(seq)):
       if seq[r] == "N":
          Nindex.append( r )

   # Construct chromosomes
   population = [] #{} modified on 2018-08-14
   for i in range(nseq):
       population.append( Chrom(len(Nindex)) ) # modified
       population[i].assign() # modified 


   # Use template [optional] 
   if tmpf != "": # use template 
      with open(tmpf) as f:
           seq_tmp = f.readline().strip()
      #print seq_tmp

      # verify the consistency of 'seq_tmp' and 'seq', 
      # and also set up the chromosomes by using the initial guess 
      for i in range(nwaive): # make the first 'nwaive' chromosomes to take the initial guess
         population[i].chrom = chkSeq(seq, seq_tmp, Nindex)
      #sys.exit()


   # Calculate initial fitness
         #print population[101].chrom
   population = calcFit( population, seq, ss, Nindex, 1, nproc )
         #print population[101].chrom
 
   # Find current best chromosome
   bestFitList = []
   bestChrom = findBest( population )
   bestFitList.append( bestChrom[1] )
   #print bestFitList

   # Calculate average fitness
   #print len(population)
   meanFitList = []
   meanFitList.append( meanFit( population ) )

   # get the waive list - {is this truly waive list?} 
   waiveList = getWaive( population, nwaive)  


   # PRINT 
   #print "Mean:", int(meanFit( population )), "Max:",bestChrom[1]
   sys.stdout.write("Mean: "+str(int(meanFit( population )))+" "+"Max: "+str(bestChrom[1])+'\r' )
   sys.stdout.flush()
   print ""

   # prepare the heaven list
   heavenList = []


   # Start Iteration
   t = 0
   nochange = 0 
   nboost = 200 # was 0
                                 #for t in range(niter): # might be better to change into while loop 
   while t <= niter:

       # check timer 
       timer_b = time.time()
       if timer_b-timer_a > total_time:
          print "\nTime is up."
          sys.exit()

       # if not sufficient...
       if t == niter and len(heavenList) < nsurvivors:
          t = 0
          nstepheaven = nstepheaven + 2 # increase the Hamming distance to the target secondary structure 

       # step size - incremental  
       t = t + 1

       # 1. selection
       population = select( population, nreplace) # improved - 08/16/2018 

       # 2. cross-over
       population = xover1( population, len(Nindex), acr_ave, nx, waiveList)
       
       # inject one template sequence if possible - 09/06/2018 
       if tmpf != "":
          population[0].chrom = chkSeq(seq, seq_tmp, Nindex) 
          pass 

       ## update the fitness 
       population = calcFit( population, seq, ss, Nindex, 1, nproc )

       ## update waiveList ?? 
       waiveList = getWaive( population, nwaive)

       ## get the ordering of fitness 
       orderC = getOrder( population )

       # 3. mutation
       population = mutation1( population, len(Nindex), orderC, mut_good, mut_ave, mut_bad, waiveList)

       # 4. calculate fitness
       population = calcFit( population, seq, ss, Nindex, 1, nproc )

       ## update waiveList ?? 
       waiveList = getWaive( population, nwaive)

       # print out current iteration result
       bestChrom = findBest( population )
       bestFitList.append( bestChrom[1] ) 
       meanFitList.append( meanFit( population ) )
       #print "Mean:", int(meanFit( population )), "Max:",bestChrom[1], str(t+1)+"/"+str(niter)
       sys.stdout.write("Mean: "+str(int(meanFit( population )))+" "+"Max: "+str(bestChrom[1])+" "+str(t+1)+"/"+str(niter)+'\r' )
       sys.stdout.flush()

       # 5. lift eligible chromosomes into heaven - modified on 09/06/2018
       population, heavenList = heaven( population, len(seq), nstepheaven, heaven_rate, heavenList )


       # mutate chromosomes in waiveList if the heavenList does not change for a while  (stuck somewhere)
       if t == 1: # first iteration 
          old_heavenlen = len(heavenList)
                        #old_heavenList = heavenList[:]

       if t > 1:
          if old_heavenlen != len(heavenList):   #old_heavenList != heavenList:
             old_heavenlen = len(heavenList) # update the heavenList length 
                                
             #old_heavenList = heavenList[:]
             nochange = 0 # go back to 0 in time, it means 'heavenList' has been changed  
          else: # 

             if old_heavenlen == 0: # nothing in the 'heavenList'
                nochange = nochange + 1 # set up a counter  
                if nochange >= nstill_0: # was set to 150 steps 
                   # to do mutation
                                     #population = stuckMutation( population, waiveList )
                   # re-initialize the whole population (99%)
                   for tmpl in range(1,nseq): #range(nseq):
                       if 0.01 < random.random(): # eliminate 99%  -- ??? 
                                                  # '>' changed in to '<' 
                          population[tmpl].assign()
                   nochange = 0
                   t = nboost
                   nboost = nboost + 100 # was 100

                   #
                   print "\npopulation re-initialized", 't =',t, 'nboost =', nboost 
                   # 


             if old_heavenlen > 0: # there are some seq. in the 'heavenList'
                #
                nochange = nochange + 1
                # if nochange exceeds a limit, then trigger mutation, also set nochange back to 0
                if nochange >= nstill_1: # was set to 100 steps 
                   # mutation
                   if len(heavenList) > 10:
                      population = stuckMutation( population, waiveList )
                      print "\nstuckMutation activated."
                   else:
                      
                      print "\npopulation re-initialized" 
                      for tmpl in range(1,nseq): # modified
                          if 0.01 < random.random():  
                             population[tmpl].assign()
                      #t = 0 -> leave for next 'nstepheaven' value      
                   nochange = 0
                                             

       ## print out heaven results
       if len(heavenList) != 0:
          if (t+1)%nprintheaven == 0: 
             f1 = open('heaven.txt','w')
             fullseq = list(seq)
             for ih in range(len(heavenList)):
                 newseq = []
                 for j in range(len(fullseq)):
                     if j not in Nindex:
                        newseq.append( fullseq[j] )
                     else:
                        newseq.append( heavenList[ih][0][Nindex.index(j)] )
                 newseq = ''.join(newseq)
                 f1.write(newseq)
                 #for h in heavenList[ih][0]:
                 #    f1.write( str(h)+" " )
                 f1.write(" ")
                 f1.write( str(heavenList[ih][1])+"/"+str(len(fullseq)) +"\n" )
             f1.close()
             pass
             # stop running the program
             if len(heavenList) >= nsurvivors:
                print "\nSufficient chromosomes in heaven:",len(heavenList)
                sys.exit()

       #




       # 5. lift eligible chromosomes into heaven
       #population, heavenList = heaven( population, len(seq), nstepheaven, heaven_rate, heavenList )







   # Print Final Result -> This is actually not called and used.
   f1 = open('result.txt','w')
   fullseq = list(seq)
   for ic in range(len(population)):
       newseq = []
       for j in range(len(fullseq)):
           if j not in Nindex:
              newseq.append( fullseq[j] )
           else:
              newseq.append( population[ic].chrom[Nindex.index(j)] )  
       newseq = ''.join(newseq)
       f1.write(newseq)
       f1.write(" ")
       f1.write( str(population[ic].fitness)+"/"+str(len(fullseq)) +"\n" )
   f1.close()

   return 0

   




if __name__== "__main__":

   if len(sys.argv) == 1:
      print "incorrect excution..."
      sys.exit()
   
   inpf = sys.argv[1]
   if not os.path.isfile(inpf):
      print "input file not exist..."
      sys.exit()

   tmpf = ""
   if len(sys.argv) == 3:
      tmpf = sys.argv[2] 
      if not os.path.isfile(tmpf):
         print "template file not exist..."
         sys.exit()



   # set up
   #inpf = "7_4-Top243.rnafold.rnaInverseInp"
   nseq = 500 # number of starting candidates (population size)
   nreplace = 50 # kill the last 'nreplace', replace them with the first 'nreplace'
   nwaive = 10 # waive any changes for 'nwaive' best chromosomes 
   
   #tmpf = "Top243.seq" 


   nstepheaven = 0 # if a chromosome has a distance to target less than 'nstepheaven', then pick it out
   heaven_rate = 0.75 # probability to be mutated after being lifted to heaven
   nprintheaven = 5 # print out heaven results every 'nprintheaven' steps if any
   nstill_1 = 100 # if the heavenList is not changed for up to 'nstill' iterations, mutate chromosomes in waiveList
   nstill_0 = 150 

   nsurvivors = 800 # 500 # when the heaven has 'nsurvivors' chromosomes, stop the program 

   niter = 500 # number of iteration

   nproc = 5 # number of CPU processors 

   mut_good = 0.30 # probability of mutation
   mut_bad =  0.75
   mut_ave =  0.50


   acr_good = 0.25 # probability of cross-over
   acr_bad =  0.25
   acr_ave =  0.25 
   nx = 4 # number of cross over genes   


   # Read optional external configuration file provided by user
   #   support variables: nstepheaven
   #   hidden input file: ".gaif_conf" in current directory
   #   This is used in 8_12 case, where S3 fails in the 'master-gaif.py'. 
   #
   if os.path.isfile(".gaif_conf"):
      with open(".gaif_conf") as f:
           for line in f:
               # set 'nstepheaven'
               if "nstepheaven" in line:
                  val = line.strip().split("=")[1]
                  nstepheaven = int(val) # if the sequence yield is too low, then set this value to a larger number 
                  print ".gaif_conf found"
                  print "nstepheaven set to", nstepheaven
   


   main(inpf,tmpf,nseq,nreplace,nwaive,niter,nproc,acr_ave,nx,mut_good,mut_ave,mut_bad,nstepheaven,heaven_rate,nprintheaven,nsurvivors,nstill_1,nstill_0)


#




