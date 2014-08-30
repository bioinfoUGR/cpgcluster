#!/usr/bin/perl -w 


###########################################################################
###########################################################################

#              ********************************
#              **  CpGcluster (version 2.0)  **
#			   **		 stand-alone		 **
#              ******************************** 

#   Laboratorio de Genómica Evolutiva y BioInformática
#    Universidad de Granada, Departamento de Genética

#              Web: http://bioinfo2.ugr.es
#              CGI: http://bioinfo2.ugr.es/CpGcluster


#  For questions, feedback, etc please contact to: José L. Oliver (oliver@ugr.es)
#                                                  Michael Hackenberg (mlhack@gmail.com)
#                                                  Francisco Dios (frankdiox@gmail.com)

# To see the options of the algorithm, please launch CpGcluster without any command line arguments
# If you use CpGcluster, please cite...

############################################################################
############################################################################
use strict;


my ($seqst, $ID, $seqlen, $seqlenbruto, $cpgnr, $prob, $output, $assembly_dir, $d, $plimit, $getd, @dd, @getd_hash ,@dist_n,$chrom_intersec,$genome_intersec);

&GetDefault();


################################# 
### Parameters ##################
#################################
my $maxN = 0; ## maximal number of Ns
#################################


my @f_assembly = split(/\//,$assembly_dir);
my @f1_assembly = split(/\./,$f_assembly[$#f_assembly]);


my @dist_all; # holds distances for the genome
my $cpg_genome; # number of CpGs in genome
my $length_genome; # number of dinucleotides in genome
if(opendir (my $DIR,$assembly_dir)){
	$cpg_genome = 0; 
	$length_genome = 0;
	
	my $output_aux = $output;
	print "\n      --- Single chromosome:\n";
	my %chrom_cods = ();
	while(my $file = readdir($DIR)){
		if($file eq "." || $file eq ".."){
			next;
		}else{
			my @ext = split (/\./,$file);
			next if($ext[$#ext] ne "fa");
		}

		print "\n      ".$file."\n";

		my $chrom = ((split(/\./,$file))[0]);
		
		$output_aux = $output;
		
		
		@dist_n = ();
		print "\n      ***                  Reading   Sequence                ***\n";
		($seqst, $ID) = &GetFA($assembly_dir."/".$file);
		print "      ***             Getting   CpG   Coordinates            ***\n";
		my @cod = &GetCoords($seqst,"CG");
		$chrom_cods{$chrom} = \@cod;
		
		print "      ***             Calculating   Seq   Features           ***\n";
		(my $num_dinuc, $cpgnr, $seqlen) = &GetSeqFeat($seqst);
		$seqlenbruto = length($seqst);
		
		if($genome_intersec){ #genome
			$cpg_genome += $cpgnr; 
			$length_genome += $num_dinuc;
		}
		

		## Prob CpG
		my $Ndach = $num_dinuc - $cpgnr; 
		$prob = $cpgnr/$Ndach;
				
		#** Begin: single chrom INTERSEC
		if($chrom_intersec){
			print "      ***          Calculating  Chrom  Intersection          ***\n";
			my ($max, $min) = getMinMaxDistance(\@dist_n, $prob);						
			$d = $max;
			$getd = "Chromosome Intersection";
			### get protoislands
			print "      ***      Detecting CpG clusters  (Chrom Intersec)      ***\n";
			my @protoislas = &GetProtoIslas(\@cod);
			print "      ***       Calculating P-values  (Chrom Intersec)       ***\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob);
			## Get Features like the obs/esp, clustering etc....
			@protoislas = &GetCGI_features(\@protoislas);

			## Writing output
			$output .= $chrom."_chromIntersec_CpGcluster.txt";
			&OUT_f(\@protoislas);
			$output = $output_aux;
		}
		
		#** End: single chrom INTERSEC
	
	
		
		
		#** Begin: single chrom PERCENTILE
		if(@getd_hash > 0){
			@dd = sort {$a <=> $b} @dist_n;
		
			foreach(@getd_hash){
				$getd = $_;
				print "      ***               Calculating Percentile $getd            ***\n";
				$d = &GetPerc(\@dd,$_);
				
				### get protoislands
				print "      ***             Detecting CpG clusters (p$getd)           ***\n";
				my @protoislas = &GetProtoIslas(\@cod);

				print "      ***              Calculating P-values  (p$getd)           ***\n";
				@protoislas = &CalcPvalNB(\@protoislas,$prob);

				## Get Features like the obs/esp, clustering etc....
				@protoislas = &GetCGI_features(\@protoislas);

				## Writing output
				$output .= $chrom."_p".$getd."_CpGcluster.txt";
				&OUT_f(\@protoislas);
				$output = $output_aux;
				
			}
		}
		
		#** End: single chrom PERCENTILE
	
		push @dist_all,@dist_n if($genome_intersec); #genome
	
	}
	


	
	#** Begin: genome
	if($genome_intersec){
		print "\n\n      --- Genome:\n\n";
		my $Ndach_genome = $length_genome - $cpg_genome; 
		my $prob_genome = $cpg_genome/$Ndach_genome;
		my ($max, $min) = getMinMaxDistance(\@dist_all, $prob_genome);

		foreach (keys %chrom_cods){
			my $chrom = $_;
			print "\n      ".$_."\n";
			$d = $max;
			$getd = "Genome Intersection";
			### get protoislands
			print "      ***      Detecting CpG clusters (Genome Intersec)      ***\n";
			($seqst, $ID) = &GetFA($assembly_dir."/".$chrom.".fa");
			my @protoislas = &GetProtoIslas($chrom_cods{$chrom});
			print "      ***       Calculating P-values (Genome Intersec)       ***\n\n\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob_genome);
			## Get Features like the obs/esp, clustering etc....
			@protoislas = &GetCGI_features(\@protoislas);
			
			## Writing output
			$output .= $chrom."_genomeIntersec_CpGcluster.txt";
			&OUT_f(\@protoislas);
			$output = $output_aux;		
		}
	}
	
	closedir($DIR);
}else{
	die "Cannot open $assembly_dir\n"
}

#** End: genome



####################################################################
#######   SUBFUNCTIONS   ###########################################
###################################################################

sub GetDefault{
  print "\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------                                                         ---------\n";
  print "---------   Laboratorio de Genomica Evolutiva y BioInformatica    ---------\n";
  print "---------    Universidad de Granada, Departamento de Genetica     ---------\n";
  print "---------                                                         ---------\n";
  print "---------               Web: http://bioinfo2.ugr.es               ---------\n";
  print "---------        CGI: http://bioinfo2.ugr.es/CpGcluster           ---------\n";
  print "---------                                                         ---------\n";
  print "---------              CpGcluster (2.0) 11/30/11                  ---------\n";
  print "---------                                                         ---------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "\n";

  
	if($#ARGV < 2){
	    print "Example for the usage of CpGcluster:\n\n";
	    print "perl CpGcluster.pl <assembly>  <d>  <P-value>\n\n";
	    
	    print "\nassembly:   Directory containing sequence files in FASTA format\n";
	
	    print "\nd:          The threshold distance on basis of a given percentile.\n";
	    print "            For example: d=25 calculates the percentile 25 of the genomic\n";
	    print "            CpG distance distribution and takes this value as the threshold\n";
	    print "            distance\n";
	    print "            The recommended value is 50 (median distance)\n";
		print "            You can add multiple comma-separated percentile values, \"ci\"\n";
		print "            (chromosome intersection) or \"gi\" (genome intersection)\n";
		print "            Example: gi,25,60,ci,50\n";
	    print "\nP-value:    The maximal P-value under which a CpG cluster is considered as a\n";
	    print "            CpG island\n";
	    print "            The recommended limit is 1E-5\n\n";
	
	    die "\n";
	}

	#ARGV0 - assembly/directory
    $assembly_dir = $ARGV[0];
    if(!(-e $assembly_dir)){
		die "Cannot find the input directory: $assembly_dir\n";
    }

    #ARGV1 - Percentile/ci/gi
	$getd = $ARGV[1];
	foreach(split(',',$ARGV[1])){
		if(/\D/){
			if(lc eq 'ci'){
				$chrom_intersec = 1;
			}elsif(lc eq 'gi'){
				$genome_intersec = 1;
			}
		}elsif($_ < 0 or $_ > 100){
			die "The Percentile must be between 0 and 100\n";
		}else{
			push @getd_hash,$_;
		}
	}

  	#ARGV2 - pLimit
    $plimit = $ARGV[2];
	die "The maximal P-value you have choosen is higher than 1!\nPlease revise the order of the input parameters\n" if($plimit > 1);


	my @f = split(/\//,$assembly_dir);
	my @f1 = split(/\./,$f[$#f]);
	$f[$#f] = "$f1[0]";
	$output = join('/',@f);
	$output .= "/";

 
}

# Read Sequence
sub GetFA{

  my $seqst_temp = "";
  open (I,$_[0]) or die "Can't open $_[0]";
  my $z = <I>;
  my $tes = substr($z,0,1);
  if($tes ne ">"){
    die "Sequence seems not to be in fasta format !!!!";
  }
  my @z = split(/\s+/," $z");
  $z = $z[1];
  $z=~ s/>//;
  $z =~ s/[\n\t\f\r\s]//g;
  my $ID_temp = $z;
  while($z = <I>){
    $z =~ s/[\n\t\f\r_0-9\s]//g;
    $seqst_temp .= $z;
  }
  return ($seqst_temp,$ID_temp);
}

sub OUT_f {
  my $c=0;

  open (OO,">$output") or die "could not open $output";
  print OO "CGI\tFrom\tTo\tLength\tCount\tOEratio\t%G+C\tPatDen\tPValue\tlogPValue\n";
 
  
  while($_[0]->[$c]){
	my $log_pvalue = ($_[0]->[$c]->[8] == 0 ? 0 : (log($_[0]->[$c]->[8])/log(10)));
	my $patden = ($_[0]->[$c]->[3]/$_[0]->[$c]->[2]);
	my $gc = ($_[0]->[$c]->[7]/$_[0]->[$c]->[2]);
	printf OO "%i\t%i\t%i\t%i\t%i\t%.3f\t%.2f\t%.3f\t%.2e\t%.2f\n",$c+1, $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[0]->[$c]->[2], $_[0]->[$c]->[3], $_[0]->[$c]->[4], $gc*100, $patden, $_[0]->[$c]->[8],$log_pvalue;
    $c++;
  }
  close(OO);
  open(O,">$output-log.txt") or die "can't open $output-log.txt";
  print O "Basic statistics of the input sequence: $ID\n";
  printf O "Length: %d\n",$seqlen;
  printf O "Length without Ns: %d\n",$seqlenbruto;
  my $fg = $seqst =~ s/g/g/ig;
  my $fc = $seqst =~ s/c/c/ig;
  my $fa = $seqst =~ s/a/a/ig;
  my $ft = $seqst =~ s/t/t/ig;
  printf O "GC content: %0.3f\n",100*($fg+$fc)/$seqlen;
  printf O "Number of CpGs in sequence: %d\n",$cpgnr;
  printf O "Probability to find a CpG: %.4f\n\n",$prob;
  print O "Parameters used:\n";

  printf O "p-value threshold: $plimit\n";
  if($getd){
    print O "Distance threshold method: ";
    $_ = $getd;
    print O "percentile " if(/\d/);
    print O "$getd\n";
  }

  close(O);
}


## Get CpG cluster
sub GetProtoIslas{

  my @coord = @{$_[0]};
  my @t;
  my ($start, $end);
  my $des = "no";
  for(my $i = 0; $i <= $#coord - 1; $i++){
    
    my $dist = $coord[$i+1]  - ($coord[$i] + 1);

    if($dist <= $d){
      if($des eq "no"){
	$start = $coord[$i];
      }
      $end = $coord[$i+1]+1;
      $des = "yes";
    }
    elsif($dist > $d and $des eq "yes"){
      $des = "no";
      my @f = ($start, $end);
      push @t,\@f;
    }
  }
  if($des eq "yes"){
    my @f = ($start, $end);
    push @t,\@f;
  }
  return @t;
}

sub GetCGI_features{

  my @temp;
  my $c=0;
  while(defined($_[0]->[$c])){
    if($_[0]->[$c]->[4] < $plimit){

      my $len = $_[0]->[$c]->[1] - $_[0]->[$c]->[0] +1;
      (my $oe, my $cpgseq, my $gccont)= &CalcObsEsp($_[0]->[$c]->[0] -1,$len);
      my $coord1 = &GetCoord($cpgseq);
      (my $clust, my $meandist) = &GetClust($coord1);
      
      my $pval = $_[0]->[$c]->[4];
      $_[0]->[$c]->[4] = $oe;
      $_[0]->[$c]->[5] = $meandist;
      $_[0]->[$c]->[6] = $clust;
      $_[0]->[$c]->[7] = $gccont;
      $_[0]->[$c]->[8] = $pval;
      my @t = @{$_[0]->[$c]};
      push @temp,\@t;
    }
    $c++;
  }
  return @temp;
}

sub CalcObsEsp{

  my $cpgseq = substr ($seqst,$_[0],$_[1]);
  my $fc = $cpgseq =~ s/c/c/ig;
  my $fg = $cpgseq =~ s/g/g/ig;
  my $CGICpG = $cpgseq =~ s/CG/CG/ig;
  my $e = $fc*$fg;
  my $oet = $CGICpG*length($cpgseq)/$e;
  return ($oet,$cpgseq,$fc+$fg);
}

sub GetSeqFeat{

  my $n = $_[0];
  my $CpGnr = $n =~ s/CG/Cg/ig;
  my $NN = $n =~ s/N/N/ig;  
  my $seqlen = length($n)-$NN;
  my $num_dinuc = 0;


  my $elem;
  my $limit = length($n) - 1;
  for(my $i = 0; $i < $limit; $i++){
  	$elem = substr($n,$i,2);
  	if($elem !~ m/N/i){
  		$num_dinuc++;
  	}
  }
  
  return ($num_dinuc,$CpGnr,$seqlen);
}


sub GetCoords{

  my $n = $_[0];
  
  $n.="j";
  my @f =split(/$_[1]/i,$n);
  my @t;

  my $lencount = 0;
  for(my $i = 0; $i < $#f; $i++){
    $lencount += length($f[$i]);
    $t[$i] = $lencount + 1;
    $lencount+=2;
    my $nnr = $f[$i] =~ s/n/n/ig;
    if($nnr <= $maxN){
      push @dist_n,length($f[$i])+1;
    }
  }
  return @t;
}

sub GetPerc{

  my @t = @{$_[0]};
  my $totnr = @t;
  for(my $i = 0; $i <= $#t; $i++){
    if(100*$i/$totnr >= $_[1]){
      return $t[$i];
    }
  }
  return $t[$#t];
}


##########################################################################
############## SUBFUNCTIONS for P-value calculations  ####################
##########################################################################

######################################################################
## *** Calculates the negative binomial distribution

sub CalcPvalNB{

  my ($pval, @temp);
  my $c=0;
  my @islas = @{$_[0]};
  for(my $i = 0; $i <= $#islas; $i++){

    my $l = $islas[$i]->[1] - $islas[$i]->[0] + 1;
    my $str = substr ($seqst,$islas[$i]->[0]-1,$l);
    my $cpg = $str =~ s/cg/cg/ig;
    my $pval = &GetNB($l-(2*$cpg),$cpg-1,$_[1]);
    $pval = sprintf("%.5e",$pval);
    my @t = ($islas[$i]->[0],$islas[$i]->[1],$l,$cpg,$pval);
    push @temp, \@t;;
    $c++;
  }
  return @temp;
}
sub GetNB{
  
  my $pval = 0;
  for(my $j = 0; $j <= $_[0]; $j++){
    my $ptemp = &FactorialNB($j,$_[1]) + $_[1]*log($_[2]) + $j*log(1.0-$_[2]);
    $ptemp = exp($ptemp);
    $pval += $ptemp;
    }
  return $pval;
}


sub FactorialNB{
  my $stop = $_[0]+$_[1]-1;
  my $l1 = 0;
  my $l2 = 0;
  for(my $i = $_[0]+1;$i <= $stop; $i++){
    $l1 +=log($i);
  }
  for(my $i = 1;$i <= $_[1]-1; $i++){
    $l2 +=log($i);
  }
  return $l1-$l2;
}
##########################################################################
############## SUBFUNCTIONS for Clustering Calculation  ####################
##########################################################################

sub GetClust{

  my @t = @{$_[0]};
  my $dist = &GetDist($_[0]);
  my $mean = &Normalize($dist);
  return (1, $mean);
}

sub Normalize{

  my @d = @{$_[0]};
  my $tot;
  for(my $i = 0; $i <= $#d; $i++){
    $tot+=$d[$i];
  }
  my $mean = $tot/@d;
  return $mean;
}

sub GetDist{

  my @dist;
  my @d = @{$_[0]};
  for(my $i = 0; $i < $#d; $i++){
    my @f = split (/\-/,$d[$i]);
    my @s = split (/\-/,$d[$i+1]);
    my $dist = $s[0]-$f[1];
    push @dist,$dist;
  }
  return \@dist;
}

sub GetCoord{

  my @coord;
  my @c = split (//,$_[0]);
  for(my $i = 0; $i < $#c; $i++){
    if($c[$i] =~ /[cC]/ and $c[$i+1] =~ /[gG]/){
      my $str = $i.'-'.eval($i+1);
      push @coord,$str;
    }
  }
  return \@coord;
}


sub getMinMaxDistance{
	my @distances = @{$_[0]};
	my $prob = $_[1];
	
	my @distCount = ();
	my $maxDist = 0;
	my $nrDistances = 0;
	my $stop = @distances;
	
	for(my $i = 0; $i < $stop; $i++){
		$nrDistances++;
		my $dist = $distances[$i];
	
		if(defined $distCount[$dist]){
			$distCount[$dist]++;
		}else{
			$distCount[$dist] = 1;
		}
		
		if($dist > $maxDist){
			$maxDist = $dist;
		}
	}
	
	my $obsCum = 0;
	my $teoCum = 0;
	my @dif = ();
	
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $distCount[$i]){
			$obsCum += $distCount[$i]/$nrDistances;
		}
		$teoCum += $prob * (1 - $prob)**($i-1);
		$dif[$i] = ($obsCum - $teoCum);
	}
	
	my $max = -1;
	my $min = 1;
	my $maxD = 0;
	my $minD = 0;
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $dif[$i]){
			my $difT = $dif[$i];
			if($difT > $max){
				$max = $difT;
				$maxD = $i;
			}
			elsif($difT < $min){
				$min = $difT;
				$minD = $i;
			}
		}
	}
	
	return ($maxD, $minD);
	
}