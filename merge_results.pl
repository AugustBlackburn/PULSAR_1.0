#!/usr/local/bin/perl

open(INFILE1,"@ARGV[0]")||die;
$header1=<INFILE1>;
chomp $header1;
@headertoks1=split(",",$header1);
$i=0;
while($line1=<INFILE1>){
	chomp $line1;
	@toks1=split(",",$line1);
	for($j=0;$j<scalar(@headertoks1);$j++){
		$columnname=@headertoks1[$j];
		@{$PBRVhaps{$columnname}}[$i]=@toks1[$j];
	}
	$i++;
}
close INFILE1;

open(INFILE2,"@ARGV[1]")||die;
$header2=<INFILE2>;
chomp $header2;
@headertoks2=split(",",$header2);
$i=0;
while($line2=<INFILE2>){
	chomp $line2;
	@toks2=split(",",$line2);
	for($j=0;$j<scalar(@headertoks2);$j++){
		$columnname=@headertoks2[$j];
		@{$SHAPEIThaps{$columnname}}[$i]=@toks2[$j];
	}
	$i++;
}
close INFILE2;

for($i=3;$i<scalar(@headertoks1);$i++){
	@toks3=split("_",@headertoks1[$i]);
	$id=@toks3[0];
	@{$idhash{$id}}[0]++;
}
for($i=3;$i<scalar(@headertoks2);$i++){
	@toks3=split("_",@headertoks2[$i]);
	$id=@toks3[0];
	@{$idhash{$id}}[0]++;
}
@uniqueids=keys %idhash;

#for each unique id
foreach $id (@uniqueids){
	#print STDOUT "$id\n";
	#establish the two chromosome names for this person
	$id1="$id"."_1";
	$id2="$id"."_2";
	#set up an initialization variable that can be switched once the alleles for the initial heterozygous genotype are assigned to haplotypes
	$initialized=0;
	#traverse through the markers along the chromosome
	for($i=0;$i<scalar(@{$PBRVhaps{Position}});$i++){
		#if the PBRV results say that this person is homozygous '0' assign '0' to each haplotype
		if(@{$PBRVhaps{$id1}}[$i] eq '0' & @{$PBRVhaps{$id2}}[$i] eq '0'){
			@{$COMBINEDhaps{$id1}}[$i]=0;
			@{$COMBINEDhaps{$id2}}[$i]=0;
		#if the PBRV results say that this person is homozygous '1' assign '1' to each haplotype
		}elsif(@{$PBRVhaps{$id1}}[$i] eq '1' & @{$PBRVhaps{$id2}}[$i] eq '1'){
			@{$COMBINEDhaps{$id1}}[$i]=1;
			@{$COMBINEDhaps{$id2}}[$i]=1;
		#if the PBRV results say this person is heterozygous
		}elsif(@{$PBRVhaps{$id1}}[$i]+@{$PBRVhaps{$id2}}[$i]==1){
			#print STDOUT "Used PBRV.\n";
			#if the phase of the haplotypes has already been initialized
			if($initialized eq '1'){
				#if the previous heterozygous genotype was also phased by PBRV
				#print STDOUT "$TRUSTPBRV\n";
				if($TRUSTPBRV eq '1'){
					#print STDOUT "Used PBRV.\n";
					#if the '1' allele was on the "_1" haplotype for the last heterozygous marker in the PBRV results
					if($PBRVlast eq '1'){
						#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the PBRV results they are on the same haplotype otherwise they aren't set $same to reflect accordingly
						if(@{$PBRVhaps{$id1}}[$i] eq '1'){
							$same=1;
						}else{
							$same=0;
						}
					#the '1' allele was on the "_2" haplotype
					}else{
						#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the PBRV results they are on opposite haplotypes otherwise they are are on the same haplotype, set $same to reflect accordingly
						if(@{$PBRVhaps{$id1}}[$i] eq '1'){
							$same=0;
						}else{
							$same=1;
						}
					}
					#We have established if the "1" alleles in adjacent heterozygous genotypes are on the same haplotype, assign the new haplotypes accordingly
					if($same eq '1'){
						#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_1" haplotype for the current marker, else the opposite is true
						if($COMBINEDlast eq '1'){
							@{$COMBINEDhaps{$id1}}[$i]=1;
							@{$COMBINEDhaps{$id2}}[$i]=0;
							$COMBINEDlast=1;
						}else{
							@{$COMBINEDhaps{$id1}}[$i]=0;
							@{$COMBINEDhaps{$id2}}[$i]=1;
							$COMBINEDlast=2;
						}
					#We have established if the "1" alleles in adjacent heterozygous genotypes are on opposite haplotypes, assign the new haplotypes accordingly
					}else{
						#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_2" haplotype for the current marker, else the opposite is true
						if($COMBINEDlast eq '1'){
							@{$COMBINEDhaps{$id1}}[$i]=0;
							@{$COMBINEDhaps{$id2}}[$i]=1;
							$COMBINEDlast=2;
						}else{
							@{$COMBINEDhaps{$id1}}[$i]=1;
							@{$COMBINEDhaps{$id2}}[$i]=0;
							$COMBINEDlast=1;
						}
					}
				#The last heterozygous genotype for PBRV was not phased, if the last was phased by SHAPEIT
				}elsif($TRUSTSHAPEIT eq '1'){
					#print STDOUT "Used Shapeit.\n";
					#Make sure that shapeit thinks that the current genotype is heterozygous
					if(@{$SHAPEIThaps{$id1}}[$i]+@{$SHAPEIThaps{$id2}}[$i]==1){
						#if the '1' allele was on the "_1" haplotype for the last heterozygous marker in the shapeit results
						if($SHAPEITlast eq '1'){
							#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the shapeit results they are on the same haplotype otherwise they aren't set $same to reflect accordingly
							if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
								$same=1;
							}else{
								$same=0;
							}
						#the '1' allele was on the "_2" haplotype
						}else{
							#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the shapeit results they are on opposite haplotypes otherwise they are are on the same haplotype, set $same to reflect accordingly
							if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
								$same=0;
							}else{
								$same=1;
							}
						}
						#We have established if the "1" alleles in adjacent heterozygous genotypes are on the same haplotype, assign the new haplotypes accordingly
						if($same eq '1'){
							#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_1" haplotype for the current marker, else the opposite is true
							if($COMBINEDlast eq '1'){
								@{$COMBINEDhaps{$id1}}[$i]=1;
								@{$COMBINEDhaps{$id2}}[$i]=0;
								$COMBINEDlast=1;
							}else{
								@{$COMBINEDhaps{$id1}}[$i]=0;
								@{$COMBINEDhaps{$id2}}[$i]=1;
								$COMBINEDlast=2;
							}
						#We have established if the "1" alleles in adjacent heterozygous genotypes are on opposite haplotypes, assign the new haplotypes accordingly
						}else{
							#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_2" haplotype for the current marker, else the opposite is true
							if($COMBINEDlast eq '1'){
								@{$COMBINEDhaps{$id1}}[$i]=0;
								@{$COMBINEDhaps{$id2}}[$i]=1;
								$COMBINEDlast=2;
							}else{
								@{$COMBINEDhaps{$id1}}[$i]=1;
								@{$COMBINEDhaps{$id2}}[$i]=0;
								$COMBINEDlast=1;
							}
						}
					#neither PBRV or shapeit provide a result, assign randomly
					}elsif(@{$SHAPEIThaps{$id1}}[$i] eq '0' & @{$SHAPEIThaps{$id2}}[$i] eq '0'){
						@{$COMBINEDhaps{$id1}}[$i]=0;
						@{$COMBINEDhaps{$id2}}[$i]=0;
					#PBRV didn't phase this genotype, if SHAPEIT called it homozygous '1', fill in both alleles with '1'
					}elsif(@{$SHAPEIThaps{$id1}}[$i] eq '1' & @{$SHAPEIThaps{$id2}}[$i] eq '1'){
						@{$COMBINEDhaps{$id1}}[$i]=1;
						@{$COMBINEDhaps{$id2}}[$i]=1;
					#PBRV didn't phase this genotype, if SHAPEIT thinks this genotype is heterozygous
					}else{
						print STDOUT "This shouldn't happen, but it did. #1\n";
						if(@{$PBRVhaps{$id1}}[$i] eq '1'){
							@{$COMBINEDhaps{$id1}}[$i]=1;
							@{$COMBINEDhaps{$id2}}[$i]=0;
							$COMBINEDlast=1;
						}else{
							@{$COMBINEDhaps{$id1}}[$i]=0;
							@{$COMBINEDhaps{$id2}}[$i]=1;
							$COMBINEDlast=2;
						}
					}
				#neither PBRV or shapeit provide a result, assign randomly
				}else{
					print STDOUT "This shouldn't happen, but it did. #2\n";
					if(@{$PBRVhaps{$id1}}[$i] eq '1'){
						@{$COMBINEDhaps{$id1}}[$i]=1;
						@{$COMBINEDhaps{$id2}}[$i]=0;
						$COMBINEDlast=1;
					}else{
						@{$COMBINEDhaps{$id1}}[$i]=0;
						@{$COMBINEDhaps{$id2}}[$i]=1;
						$COMBINEDlast=2;
					}
				}	
			#Phase for the haplotypes has not been initialized
			}else{
				if(@{$PBRVhaps{$id1}}[$i] eq '1'){
					@{$COMBINEDhaps{$id1}}[$i]=1;
					@{$COMBINEDhaps{$id2}}[$i]=0;
					$COMBINEDlast=1;
				}elsif(@{$PBRVhaps{$id2}}[$i] eq '1'){
					$COMBINEDlast=2;
					@{$COMBINEDhaps{$id1}}[$i]=0;
					@{$COMBINEDhaps{$id2}}[$i]=1;
				}else{
					print STDOUT "Houston we have a problem.\n";
				}
				$initialized=1;
			}
			$TRUSTPBRV=1;
			if(@{$PBRVhaps{$id1}}[$i] eq '1'){
				$PBRVlast=1;
			}else{
				$PBRVlast=2;
			}
			if(@{$SHAPEIThaps{$id1}}[$i]+@{$SHAPEIThaps{$id2}}[$i]==1){
				$TRUSTSHAPEIT=1;
				if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
					$SHAPEITlast=1;
				}else{
					$SHAPEITlast=2;
				}
			}else{
				$TRUSTSHAPEIT=0;
			}
		#PBRV didn't phase this genotype, if SHAPEIT called it homozygous '0', fill in both alleles with '0'
		}elsif(@{$SHAPEIThaps{$id1}}[$i] eq '0' & @{$SHAPEIThaps{$id2}}[$i] eq '0'){
			@{$COMBINEDhaps{$id1}}[$i]=0;
			@{$COMBINEDhaps{$id2}}[$i]=0;
		#PBRV didn't phase this genotype, if SHAPEIT called it homozygous '1', fill in both alleles with '1'
		}elsif(@{$SHAPEIThaps{$id1}}[$i] eq '1' & @{$SHAPEIThaps{$id2}}[$i] eq '1'){
			@{$COMBINEDhaps{$id1}}[$i]=1;
			@{$COMBINEDhaps{$id2}}[$i]=1;
		#PBRV didn't phase this genotype, if SHAPEIT thinks this genotype is heterozygous
		}elsif(@{$SHAPEIThaps{$id1}}[$i]+@{$SHAPEIThaps{$id2}}[$i]==1){
			#print STDOUT "Used Shapeit.\n";
			#if the phase of the haplotypes has already been initialized
			if($initialized eq '1'){
				if($TRUSTSHAPEIT eq '1'){
					#if the '1' allele was on the "_1" haplotype for the last heterozygous marker in the shapeit results
					if($SHAPEITlast eq '1'){
						#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the shapeit results they are on the same haplotype otherwise they aren't set $same to reflect accordingly
						if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
							$same=1;
						}else{
							$same=0;
						}
					#the '1' allele was on the "_2" haplotype
					}else{
						#if the '1' allele is on the "_1" haplotype for the current heterozygous marker in the shapeit results they are on opposite haplotypes otherwise they are are on the same haplotype, set $same to reflect accordingly
						if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
							$same=0;
						}else{
							$same=1;
						}
					}
					#We have established if the "1" alleles in adjacent heterozygous genotypes are on the same haplotype, assign the new haplotypes accordingly
					if($same eq '1'){
						#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_1" haplotype for the current marker, else the opposite is true
						if($COMBINEDlast eq '1'){
							@{$COMBINEDhaps{$id1}}[$i]=1;
							@{$COMBINEDhaps{$id2}}[$i]=0;
							$COMBINEDlast=1;
						}else{
							@{$COMBINEDhaps{$id1}}[$i]=0;
							@{$COMBINEDhaps{$id2}}[$i]=1;
							$COMBINEDlast=2;
						}
					#We have established if the "1" alleles in adjacent heterozygous genotypes are on opposite haplotypes, assign the new haplotypes accordingly
					}else{
						#if the '1' allele was on the "_1" haplotype previously, put the '1' allele on the "_2" haplotype for the current marker, else the opposite is true
						if($COMBINEDlast eq '1'){
							@{$COMBINEDhaps{$id1}}[$i]=0;
							@{$COMBINEDhaps{$id2}}[$i]=1;
							$COMBINEDlast=2;
						}else{
							@{$COMBINEDhaps{$id1}}[$i]=1;
							@{$COMBINEDhaps{$id2}}[$i]=0;
							$COMBINEDlast=1;
						}
					}
				#neither PBRV or shapeit provide a result, assign randomly
				}else{
					print STDOUT "This shouldn't happen, but it did. #3\n";
					if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
						@{$COMBINEDhaps{$id1}}[$i]=1;
						@{$COMBINEDhaps{$id2}}[$i]=0;
						$COMBINEDlast=1;
					}else{
						@{$COMBINEDhaps{$id1}}[$i]=0;
						@{$COMBINEDhaps{$id2}}[$i]=1;
						$COMBINEDlast=2;
					}
				}
			
			#phase for the haplotypes has not been initialized, PBRV didn't phase this genotype, but shapeit did
			}else{
				if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
					@{$COMBINEDhaps{$id1}}[$i]=1;
					@{$COMBINEDhaps{$id2}}[$i]=0;
					$COMBINEDlast=1;
				}elsif(@{$SHAPEIThaps{$id2}}[$i] eq '1'){
					@{$COMBINEDhaps{$id1}}[$i]=0;
					@{$COMBINEDhaps{$id2}}[$i]=1;
					$COMBINEDlast=2;
				}else{
					print STDOUT "Houston we have a problem.\n";
				}
				$initialized=1;
			}
			$TRUSTPBRV=0;
			$TRUSTSHAPEIT=1;
			if(@{$SHAPEIThaps{$id1}}[$i] eq '1'){
				$SHAPEITlast=1;
			}else{
				$SHAPEITlast=2;
			}
		}
	}
}

#Print out data
open(OUTFILE,">Combined_phased_chromosomes.csv")||die;
print OUTFILE "@headertoks2[0],@headertoks2[1],@headertoks2[2]";
foreach $id (@uniqueids){
	#establish the two chromosome names for this person
	$id1="$id"."_1";
	$id2="$id"."_2";
	print OUTFILE ",$id1,$id2";
}
print OUTFILE "\n";
for($i=0;$i<scalar(@{$PBRVhaps{Position}});$i++){
	print OUTFILE "@{$PBRVhaps{Marker}}[$i],@{$PBRVhaps{Chr}}[$i],@{$PBRVhaps{Position}}[$i]";
	foreach $id (@uniqueids){
	#establish the two chromosome names for this person
	$id1="$id"."_1";
	$id2="$id"."_2";
	print OUTFILE ",@{$COMBINEDhaps{$id1}}[$i],@{$COMBINEDhaps{$id2}}[$i]";
	}
print OUTFILE "\n";
}
close OUTFILE;