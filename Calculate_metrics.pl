#!/usr/local/bin/perl

#Read in the results of the phasing software
open(INFILE1,"@ARGV[0]")||die;
$headerline1=<INFILE1>;
chomp $headerline1;
@headertoks1=split(",",$headerline1);
$j=0;
foreach $headertok1 (@headertoks1){
	if($j>2){
		@toks1=split("_",$headertok1);
		$name=@toks1[0];
		@{$names{$name}}[0]++;
	}
	$j++;
}
$j=0;
while($line=<INFILE1>){
	chomp $line;
	@toks2=split(",",$line);
	for($i=0;$i<scalar(@headertoks1);$i++){
		$curcol=@headertoks1[$i];
		@{$phasedhash{$curcol}}[$j]=@toks2[$i];
	}
	$j++;
}
$k=$j;
close INFILE1;

#Read in the answers known from simulation
open(INFILE2,"@ARGV[1]")||die;
$headerline2=<INFILE2>;
chomp $headerline2;
@headertoks2=split(",",$headerline2);
$j=0;
foreach $headertok2 (@headertoks2){
	if($j>2){
		@toks3=split("_",$headertok2);
		$name=@toks3[0];
		@{$names2{$name}}[0]++;
	}
	$j++;
}
$j=0;
while($line2=<INFILE2>){
	chomp $line2;
	@toks4=split(",",$line2);
	for($i=0;$i<scalar(@headertoks2);$i++){
		$curcol=@headertoks2[$i];
		@{$knownhash{$curcol}}[$j]=@toks4[$i];
	}
	$j++;
}
close INFILE2;

#Read in the genotype data given to the software
open(INFILE3, "@ARGV[2]")||die;
$headerline3=<INFILE3>;
chomp $headerline3;
@headertoks3=split(",",$headerline3);
$j=0;
while($line3=<INFILE3>){
	chomp $line3;
	@toks5=split(",",$line3);
	for($i=0;$i<scalar(@headertoks3);$i++){
		$curcol=@headertoks3[$i];
		@{$genotypehash{$curcol}}[$j]=@toks5[$i];
	}
	$j++
}

#Find overlap in samples between known and phased
@uniqueids1= keys %names;
@uniqueids2= keys %names2;

foreach $id1 (@uniqueids1){
	foreach $id2 (@uniqueids2){
		if($id1 eq $id2){
			push @idoverlap, $id1;
		}
	}
}

open(OUTFILE,">Switch_error_rate.csv")||die;
print OUTFILE "ID,Correct_switch(n),Incorrect_switch(n),SER,Switch_Accuracy,Percentage_phased,Hetero_percentage_phased,Homo_percentage_phased,Total_phased(n),Hetero_phased(n),Homo_phased(n),Hetero_phased_correctly(n),Homo_phased_correctly(n),Hetero_phased_incorrectly(n),Homo_phased_incorrectly(n),Fixed_genotypes(n),Hetero_genotype_fixed(n),Homo_genotype_fixed(n),Not_fixed_genotypes(n),Hetero_not_fixed(n),Homo_not_fixed(n),Total_not_phased(n),Hetero_not_phased(n),Homo_not_phased(n),Hetero_incorrect_genotype_not_phased(n),Homo_incorrect_genotype_not_phased(n),Simulated_genotype_accuracy,Correct_genotypes_given(n),Incorrect_genotypes_given(n)\n";
#print OUTFILE "$id,$correct,$incorrect,$switcherrorrate,$switchaccuracy,$totalpercentagephased,$heterozygotepercentgephased,$homozygotepercentagephased,$totalphased,$heterophased,$homophased,$heterozygotephasedcorrectly,$homozygotephasedcorrectly,$heterozygotephasedincorrectly,$homozygotephasedincorrectly,$fixedconstruction,$heterozygotephasedandfixedcorrectly,$homozygotephasedandfixedcorrectly,$didnotfixconstruction,$heterozygotephasedandfixedincorrectly,$homozygotephasedandfixedincorrectly,$totalnotphased,$heteronotphased,$homonotphased,$heteronotfixed,$homonotfixed,$correctgenogivenpercentage,$correctgenogiven,$incorrectgenogiven\n";

foreach $id (@idoverlap){
	$correct=0;
	$incorrect=0;
	$first=0;
	$chr1="$id"."_"."1";
	$chr2="$id"."_"."2";
	@known1=@{$knownhash{$chr1}};
	@known2=@{$knownhash{$chr2}};
	@phased1=@{$phasedhash{$chr1}};
	@phased2=@{$phasedhash{$chr2}};
	@genotypegiven=@{$genotypehash{$id}};
	for($i=0;$i<scalar(@known1);$i++){
		if(@known1[$i]+@known2[$i]==1){
			if(@phased1[$i]+@phased2[$i]==1){
				if(@known1[$i] eq '0' || @known1[$i] eq '1'){
					if(@known2[$i] eq '0' || @known2[$i] eq '1'){
						if(@phased1[$i] eq '0' || @phased1[$i] eq '1'){
							if(@phased2[$i] eq '0' || @phased2[$i] eq '1'){
								if($first==0){
									if("@known1[$i]" eq "@phased1[$i]"){
										$firstobs=1;
									}else{
										$firstobs=0;
									}
									$first=1;
									$last=$i;
								}else{
									if("@known1[$i]" eq "@phased1[$i]"){
										$secondobs=1;
									}else{
										$secondobs=0;
									}
									if($firstobs==$secondobs){
										$correct++;
										$last=$i;
									}else{
										$incorrect++;
										#print STDOUT "Firstobs:$firstobs\tSecondobs:$secondobs\t@known1[$last]\t@phased1[$last]\t:\t@known1[$i]\t@phased1[$i]\n";
										$firstobs=$secondobs;
										$last=$i;
									}
								}
							}
						}
					}				
				}
			}
		}
	}
	$switcherrorrate=();
	if(($correct+$incorrect)>0){
		$switcherrorrate=$incorrect/($correct+$incorrect);
	}
	$switchaccuracy=();
	$switchaccuracy=1-$switcherrorrate;
	$totalphased=();
	$heterophased=();
	$homophased=();
	$heterozygotephasedcorrectly=();
	$homozygotephasedcorrectly=();
	$heterozygotephasedincorrectly=();
	$homozygotephasedincorrectly=();
	$fixedconstruction=();
	$heterozygotephasedandfixedcorrectly=();
	$homozygotephasedandfixedcorrectly=();
	$didnotfixconstruction=();
	$heterozygotephasedandfixedincorrectly=();
	$homozygotephasedandfixedincorrectly=();
	$totalnotphased=();
	$heteronotphased=();
	$homonotphased=();
	$heteronotfixed=();
	$homonotfixed=();
	$correctgenogiven=();
	$incorrectgenogiven=();
	for($i=0;$i<scalar(@known1);$i++){
		$switch=0;
		$realgenotype=();
		$givengenotype=();
		$constructedgenotype=();
		$givengenotype=@genotypegiven[$i];
		$realgenotype=@known1[$i]+@known2[$i];
		if(@phased1[$i] eq '0' || @phased1[$i] eq '1'){
			if(@phased2[$i] eq '0' || @phased2[$i] eq '1'){
				$constructedgenotype=@phased1[$i]+@phased2[$i];
				$switch=1;
			}
		}
		if($switch==1){
			$totalphased++;
			if("$realgenotype" eq "1"){
				$heterophased++;
			}else{
				$homophased++;
			}
			if("$realgenotype" eq "$givengenotype"){
				$correctgenogiven++;
				if("$constructedgenotype" eq "$realgenotype"){
					if("$realgenotype" eq "1"){
						$heterozygotephasedcorrectly++;
					}else{
						$homozygotephasedcorrectly++;
					}
				}else{
					if("$realgenotype" eq "1"){
						$heterozygotephasedincorrectly++;
					}else{
						$homozygotephasedincorrectly++;
					}
				}
			}else{
			$incorrectgenogiven++;
				if("$constructedgenotype" eq "$realgenotype"){
					if("$realgenotype" eq "1"){
						$heterozygotephasedandfixedcorrectly++;
					}else{
						$homozygotephasedandfixedcorrectly++;
					}
				}else{
					if("$realgenotype" eq "1"){
						$heterozygotephasedandfixedincorrectly++;
					}else{
						$homozygotephasedandfixedincorrectly++;
					}
				}
			}
		}else{
			$totalnotphased++;
			if("$realgenotype" eq "$givengenotype"){
				$correctgenogiven++;
				if("$realgenotype" eq "1"){
					$heteronotphased++;
				}else{
					$homonotphased++;
				}
			}else{
				$incorrectgenogiven++;
				if("$realgenotype" eq "1"){
					$heteronotfixed++;
				}else{
					$homonotfixed++;
				}
			}	
		}
	}
	$totalpercentagephased=();
	if(($totalphased+$totalnotphased)>0){
		$totalpercentagephased=$totalphased/($totalphased+$totalnotphased);
	}
	$heterozygotepercentgephased=();
	if(($heterophased+$heteronotphased+$heteronotfixed)>0){
		$heterozygotepercentgephased=$heterophased/($heterophased+$heteronotphased+$heteronotfixed);
	}
	$homozygotepercentagephased=();
	if(($homophased+$homonotphased+$homonotfixed)>0){
		$homozygotepercentagephased=$homophased/($homophased+$homonotphased+$homonotfixed);
	}
	$correctgenogivenpercentage=();
	if(($correctgenogiven+$incorrectgenogiven)>0){
		$correctgenogivenpercentage=$correctgenogiven/($correctgenogiven+$incorrectgenogiven);
	}
	
	print OUTFILE "$id,$correct,$incorrect,$switcherrorrate,$switchaccuracy,$totalpercentagephased,$heterozygotepercentgephased,$homozygotepercentagephased,$totalphased,$heterophased,$homophased,$heterozygotephasedcorrectly,$homozygotephasedcorrectly,$heterozygotephasedincorrectly,$homozygotephasedincorrectly,$fixedconstruction,$heterozygotephasedandfixedcorrectly,$homozygotephasedandfixedcorrectly,$didnotfixconstruction,$heterozygotephasedandfixedincorrectly,$homozygotephasedandfixedincorrectly,$totalnotphased,$heteronotphased,$homonotphased,$heteronotfixed,$homonotfixed,$correctgenogivenpercentage,$correctgenogiven,$incorrectgenogiven\n";
}

