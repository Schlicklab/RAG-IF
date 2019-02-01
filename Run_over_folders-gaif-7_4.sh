#!/bin/bash 
# This bash script runs over folders to do the automatic mutations.

for folder in Top*;do
   #echo $folder
   
   # 1. Enter the folder
   cd $folder
   
   # 2. Run
   echo "Running "${folder} 
   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/Cases/subroutines/collect_input-7_4.py ${folder}

   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/Cases/subroutines/master-gaif.py ${folder} #  > FINAL_OUTPUT

   # 3. Exit this folder
   cd ..  

done 	

