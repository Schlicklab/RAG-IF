# RAG-IF

Usage details:

1) Install RNAfold and NUPACK on your system, the two 2D structure prediction packages used in RAG-IF.

2) Change the absolute paths for calling RNAfold, NUPACK, and any other scriprs code in all the code files of RAG-IF.


Step by step instructions of running the RAG-IF code are below. Example files and run are available in the directory ExampleRun/.

1) Create 2 folders one called BPSEQ, which will contain the BPSEQ files obtained from the design webserver (.bpseq extension); another folder called SEQ which will contain
the sequence files corresponding to the BPSEQ files (.seq extension). The .seq files are not provided by the RAGTOP webserver, these will have to be created.

2) Run "./run_forForna_NuPACK.sh ExampleRun/BPSEQ"
This will process the bpseq files and run NUPACK to determine the secondary structure for all input sequences. This will create .forna files
and .nupack.bpseq files in the BPSEQ directory.

3) Run "./run_forForna_RNAfold.sh ExampleRun/SEQ"
This will process the seq files and run RNAfold to determine the secondary structures for all input sequences. This will create .rnafold files and .forna files in the SEQ
directory.

4) Create files RNAfold_Rank_Topo.txt and Nupack_Rank_Topo.txt that contains the rank of the sequence (Top1.bpseq for example has rank 1) and their topology as predicted 
by RNAfold (mfe and centroid) and NUPACK. This information can be obtained from the results from the webserver.
   
5) Run "python check_can.py ExampleRun/ 7_4 > ExampleRun/candidate_results.txt". This will look at the two files created in Step 4 and categorize the sequences based 
on what the NUPACK and RNAfold predicted topologies are. The 7_4 is the target topology, therefore that will need to change depending on what your target is. Currently, 
we are interested in the category where the sequence has the target topology using RNAfold but not with NUPACK. That category is IIIb. Look at the code "check_can.py"
for category definitions.

6) Run "python filter_type.py ExampleRun/BPSEQ/ ExampleRun/candidate_results.txt > ExampleRun/result-filter.txt". This will calculate the number of common vertices between
the sequences in the category IIIb as predicted by NUPACK and the target topology.

7) Run "python make_folders.py ExampleRun". This will create folders based on number of common vertices and sequence rank to keep the results catalouged.

The next steps provide details on how to run RAG-IF for one sequence. If you want to run it on all sequences in a given folder, use the scripts 
"Run_over_folders-gaif.sh" and "run_min_mut_analysis.sh"

8) cd to the sequence folder you want to run RAG-IF on (ExampleRun/7/Top107).

9) Run "python <path to RAG-IF Code>/collect_input.py Top107 <path to ExampleRun>". This will copy all necessary files needed to run the RAG-IF GA to this directory.
Top107 is the name and rank of the sequence you want to run GA on. This will need to change according to your sequence.

10) Run "python <path to RAG-IF code>/master-gaif.py Top107". This will execute the GA. Top107 is the name and rank of the sequence you want to run GA on. 
This will need to change according to your sequence. This will create a directory Top107.survivors, which are the candidate sequences returned by GA that 
have the correct topology with both RNAfold and NUPACK. 

11) Run "python <path to RAG-IF code>/end_seq_order.py Top107.survivors > min_mut.analysis". This will run mutation optimization and give minimal mutations written
at the bottom of the file min_mut.analysis.
