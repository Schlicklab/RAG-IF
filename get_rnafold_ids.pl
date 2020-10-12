#!/usr/local/bin/perl
#ARGV[0] = diretcory where the rnafold run files are, ARGV[1] directory where RNAfold_Rank_Topo.txt should be written

# get the list of RNAfold output files
$command="ls ".$ARGV[0]."*.rnafold > ".$ARGV[0]."temp_list.txt";
system($command);

@ranks=();
@min_ids=();
@cen_ids=();
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

	#creating bpseq files from RNAfold output
        open(OUT,"<$line");
        $garbage=<OUT>;
	$sequen=<OUT>;
	$sequen =~ s/\n//g;
	$struct=<OUT>;
	@cols_s=split(/\s+/,$struct);
        $min_struct=$cols_s[0]; #dot bracket of the minimum energy structure
        $garbage=<OUT>;
	$struct=<OUT>;
        @cols_s=split(/\s+/,$struct);
        $cen_struct=$cols_s[0]; #dot bracket of the centroid energy structure
        close(OUT);
	$min_seq=$ARGV[0]."temp_min.fa";
	open(FA,">$min_seq");
	print FA ">Temp mfe\n";
	print FA $sequen."\n";
	close(FA);
        $min_txtfile=$ARGV[0]."temp_min.dotbracket";
        open(TXT,">$min_txtfile");
        print TXT $min_struct;
        close(TXT);
        $command="python dotfa2bpseq.py ".$ARGV[0]."temp_min.fa ".$ARGV[0]."temp_min.dotbracket"; #bpseq file for the minimum energy structure
        system($command);
	$min_seq=$ARGV[0]."temp_cen.fa";
        open(FA,">$min_seq");
        print FA ">Temp cen\n";
        print FA $sequen."\n";
        close(FA);
        $cen_txtfile=$ARGV[0]."temp_cen.dotbracket";
        open(TXT,">$cen_txtfile");
        print TXT $cen_struct;
	close(TXT);
        $command="python dotfa2bpseq.py ".$ARGV[0]."temp_cen.fa ".$ARGV[0]."temp_cen.dotbracket"; #bpseq file for the minimum energy structure
        system($command);
	
	#runnning treeGraphs to get the RAG IDs
	$command="python2.7 modified-treeGraph/treeGraphs.py ".$ARGV[0]."temp_min.bpseq";
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
	$command="python2.7 modified-treeGraph/treeGraphs.py ".$ARGV[0]."temp_cen.bpseq";
        $output_text=`$command`;
        @outlines=split(/\n/,$output_text);
        $topo_cen="";
        for(my $i=0; $i<scalar(@outlines); $i++){
                @cols5=split(/:/,$outlines[$i]);
                if($cols5[0] eq "Graph ID"){ # this is the graph id of the bpseq file

                        $topo_cen = $cols5[1];
                        $topo_cen =~ s/\s+//g;
                        $topo_cen =~ s/\n//g;
			@cen_ids[$count]=$topo_cen;
                        last;
                }
        }	
	$command="rm -f ".$ARGV[0]."temp_min*";
	system($command);
	$command="rm -f ".$ARGV[0]."temp_cen*";
	system($command);
}
close(LIST);

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
			$temp = $cen_ids[$i];
			$cen_ids[$i] = $cen_ids[$j];
			$cen_ids[$j] = $temp;
		}
	}
}

#writing the output file
$output_file=$ARGV[1]."RNAfold_Rank_Topo.txt";
open(OUTPUT,">$output_file");
for($i=0;$i<=$count;$i++){
	print OUTPUT $ranks[$i]."\t".$min_ids[$i]."\t".$cen_ids[$i]."\n";
}
close(OUTPUT);
