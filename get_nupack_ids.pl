#!/usr/local/bin/perl
#ARGV[0] = diretcory where the nupack run files are, ARGV[1] directory where Nupack_Rank_Topo.txt should be written

# get the list of NUPACK bpseq output files
$command="ls ".$ARGV[0]."*.nupack.bpseq > ".$ARGV[0]."temp_list.txt";
system($command);

@ranks=();
@min_ids=();
$count = -1;

$list_file=$ARGV[0]."temp_list.txt";
open(LIST,"<$list_file");
while($line=<LIST>){ # for every rnafold output file

	$count++;
	$line=~ s/\n//g;
	#getting the rank / number of the struct file
	@cols=split(/\//,$line);
	@cols2=split(/\./,$cols[2]);
	print $cols2[0]."\n";
	$cols2[0] =~ s/Top//g;
	$ranks[$count] = $cols2[0];

	#runnning treeGraphs to get the RAG IDs
	$command="python2.7 modified-treeGraph/treeGraphs.py ".$line;
        $output_text=`$command`;
        @outlines=split(/\n/,$output_text);
        $topo_min="";
        for(my $i=0; $i<scalar(@outlines); $i++){
        	@cols5=split(/:/,$outlines[$i]);
                if($cols5[0] eq "Graph ID"){ # this is the graph id of the bpseq file

                        $topo_min = $cols5[1];
                        $topo_min =~ s/\s+//g;
                        $topo_min =~ s/\n//g;
			@min_ids[$count]=$topo_min;
                        last;
                }
        }
}
#clean up
$command="rm -f ".$ARGV[0]."temp*";
system($command);

#sort accoring to numerical rank
for($i=0;$i<=$count;$i++){

	for($j=$i+1;$j<=$count;$j++){

		if($ranks[$i] > $ranks[$j]){

			$temp = $ranks[$i];
			$ranks[$i] = $ranks[$j];
			$ranks[$j] = $temp;
			$temp = $min_ids[$i];
			$min_ids[$i] = $min_ids[$j];
			$min_ids[$j] = $temp;
		}
	}
}

#writing the output file
$output_file=$ARGV[1]."Nupack_Rank_Topo.txt";
open(OUTPUT,">$output_file");
for($i=0;$i<=$count;$i++){
	print OUTPUT $ranks[$i]."\t".$min_ids[$i]."\n";
}
close(OUTPUT);
