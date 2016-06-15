#!/usr/local/bin/perl

###This is the first section of PBRV. The purpose of this section is to parse out the lineages in the pedigree. $ARGV[0] is pedigree file. 
###This program outputs 1 file: 
### All_paths.csv is a list of all of the paths an allele can take through a pedigree from each founder to all of their offspring.

#Open input and output files
#$ARGV[0] is pedigree file
open (INFILE_PED, $ARGV[0]) || die;
#open (OUTFILE, ">All_paths.csv") || die;

#Read in header line in order to discard
$header = <INFILE_PED>;
#Read in the file one line at a time, split each line, and store information in hashes and arrays
while ($trioline = <INFILE_PED>){
	#Remove end of line
	chomp $trioline;
	#Split line into an array
	@triotoks = split(',', $trioline);
	#Set the variable $currid to hold the ID for the current line
	$currid=@triotoks[0];
	#Make an array of each mom
	push @maternal, @triotoks[2];
	#Make a hash where the keys are the current id and the first array position is their mom
	push @{$matid{$currid}}, @triotoks[2];
	#Make an array of each dad
	push @paternal, @triotoks[1];
	#Make a hash where the keys are the current id and the first array position is their dad
	push @{$patid{$currid}}, @triotoks[1];
	#Make an array of the current ids, this array and the previous two arrays holding moms and dads are in the same order
	push @offspr, @triotoks[0];
	#Make an array of the sex for each person, this array is in the same order as the previous three arrays
	push @sex, @triotoks[3];
	#If the person is a founder add them to an array of the founders
	if(@triotoks[2] eq '0'){
		push @founders, @triotoks[0];
	}
}
close INFILE_PED;

#create a variable ($numberofids) that holds the number of individuals in the arrays created above
$numberofids=scalar(@offspr);
#for each founder make a list of offspring
foreach $founder (@founders){
	#blank the array @currfoundoffspring
	@currfoundoffpring=();
	#add the founder to the array
	push @currfoundoffpring, $founder;
	#set a count of the founders offspring
	$thisfoundercount = 0;
	#Call a recursion algorithm that traverses the pedigree and finds each individuals offspring, starting with the founder. For each kid it adds to the array and count directly above.
	&find_paths($founder,$founder);
	#Get the length of the array just created, it will change from founder to founder
	$toks1length = scalar(@currfoundoffpring);
	#Add each offspring to the hash %founder, with the current founder as the key.
	for ($i=1;$i<$toks1length;$i++){
		push @{$founder{$currfoundoffpring[0]}}, @currfoundoffpring[$i];
	}
	
	#insert @currfounderoffpring
}
#close OUTFILE 
#close OUTFILE;

###This is the second section of PBRV. The purpose of this section is to identify a subset of variants in which one allele is lineage specific (ie. introduced into the pedigree by a single founder).
### This program outputs 1 file:
###	Potential_private_variants.rv

#$ARGV[1] is a genotype file in HIPster format
open (INFILE2, $ARGV[1])|| die;
#Open the outfile that will have a list of the potentially lineage specific variants
open (OUTFILE1, ">Potential_private_variants.rv") || die;
#for each lineage
foreach $lineage (@lineages){
	#Split the path into an array containing the individuals in this path
	@toks = split('_',$lineage);
	#Get the length of the array
	$lineagelength=scalar(@toks);
	#If there are individuals between the founder and the last person on this path, store these individuals in a hash with a compression of the founder and the individual as the key.
	if($lineagelength>2){
		#Compress the founder and the individual at the end of the path into a string to use as a key in the hash
		$firstlastcombo=@toks[0].@toks[$lineagelength-1];
		#Store the individuals between the founder and the last person in a hash with $firstlastcombo as the key to the hash.
		@{$pathholder{$firstlastcombo}}=@toks[1..($lineagelength-2)];
	}
	#The founder is the first spot in the array
	$founder=@toks[0];
	#The last person in the lineage
	$laster=@toks[scalar(@toks)-1];
	#Add the current line to the end of an array in the hash %last where the key is the last person in the lineage
	push @{$last{$laster}}, $lineage;
	#Add the founder in the current line to the end of an array in the hash %founderlist where the key is the last person in the lineage
	push @{$founderlist{$laster}},$founder;
	#Create a sting that joins the founder and the last person in the lineage.
	$combo=$founder.$laster;
	#For the combo of founder and lineage store the current lineage line in the hash %combohash where the key is the string created in the last line.
	@{$combohash{$combo}}[0] = $lineage;
}

#Read in the header line of the file with the genotypes in HIPster format
$header2 = <INFILE2>;
#Print the headerline to the new output file. The header is the same in both files.
print OUTFILE1 "$header2";
#Remove the end of line
chomp $header2;
#Split the line and store in the array @toks2
@toks2 = split(',',$header2);
#Find out the length of @toks2
$toks2length=scalar(@toks2);
#Create a hash called %genotypes that has the ids as the keys
#Read in the file and store the genotypes in a hash called %genotypes
while($line2=<INFILE2>){
	#Remove the end of line
	chomp $line2;
	#Split the line
	@toks3=split(",",$line2);
	#For each column, store the appropriate data from the current line in a hash with the column name as the key
	for($i=0;$i<$toks2length;$i++){
		$id=@toks2[$i];
		@{$genotypes{$id}}[$j]=@toks3[$i];
	}
	$j++;
}
$linecount=$j;
close INFILE2;

#Iterate through each variant in the %genotypes hash
for($i=0;$i<$linecount;$i++){
	#Establish counts that will be used to count the states of the genotypes
	$count0=0;
	$count1=0;
	$count2=0;
	#establish a switch that will be used to determine if the genotype line has already been printed, this will be reset each line and insures that each line is only printed once
	$switch=0;
	$switch2=0;
	#Iterate through the genotypes 
	for($j=3;$j<$toks2length;$j++){
		$id=@toks2[$j];
		if (@{$genotypes{$id}}[$i] == 0){
			$count0++;
			#print STDOUT "Increased 0 count\n";
		}elsif(@{$genotypes{$id}}[$i] == 1){
			#print STDOUT "Increased 1 count\n";
			$count1++;
		}elsif(@{$genotypes{$id}}[$i] == 2){
			#print STDOUT "Increased 2 count\n";
			$count2++;
		#if 0,1,or 2 are not observed there are missing genotypes, skip this line, proceed to the next line
		}else{
			###Change the switch2 to 1, exit this for loop
			$switch2=1;
			last;
		}
	}
	#if there was a genotype missing, $switch2 will be set to 1, move forward to the next line in the file
	if($switch2==1){
		#print STDOUT "Switch 2\n";
		next;
	}
	#If both homozygous states are present then this is a common variant, is not lineage specific. If both are present just continue to next line.
	if($count0==0||$count2==0){
		#cycle through the list of founders
		foreach $founderid (@founders){
			#Establish a switch that indicates if their is a genotype issue with specific lineages
			$lineageproblemswitch=0;
			#Establish a count for the number of heterozygotes in this founders lineage
			$foundercount=0;
			#If the founder is genotyped, do the following
			if($genotypes{$founderid}){
				#If the founder is genotyped, they must be a heterozygote
				if($genotypes{$founderid}[$i] == 1){
					#Add one to the founder count
					$foundercount++;
				#Not a heterozygote, continue to next founder
				}else{
					next;
				}
			}
			#Count how many of their offspring are heterozygotes
			foreach $offspring (@{$founder{$founderid}}){
				if($genotypes{$offspring}[$i] == 1){
					$foundercount++;
					#Check that all genotyped individuals on the lineage path between the current founder and current offspring are also heterozygotes.
					$firstlastcombo=$founderid.$offspring;
					#If there are individuals along the lineage path between the current founder and the current offspring
					if(@{$pathholder{$firstlastcombo}}){
						#Set a temporary array equal to the lineage stored in the hash %pathholder
						@currentpatharray=@{$pathholder{$firstlastcombo}};
						#Each individual that is genotyped along this path should be heterozygous. If one isn't then flip the switch to indicate that this variant is not potentially lineage specific
						foreach $lineageid (@currentpatharray){
							if(@{$genotypes{$lineageid}}){
								if($genotypes{$lineageid}[$i] == 0 || $genotypes{$lineageid}[$i] == 2){
									$lineageproblemswitch=1;
								}
							}
						}
					}
				}
			}
			#If the switch $lineageproblemswitch was set to 1, move on to the next founder
			if($lineageproblemswitch == 1){
				#print STDOUT "Lineage Switch problem\n";
				next;
			}			
			#If their lineage contains all of the heterozygotes print the genotype line to the outfile. This line only needs to be printed once, so change the switch so that if additional founders meet the requirements this genotype line isn't printed again.
			if($foundercount == $count1 & $count1 > 0){
				#the genotype line hasn't been printed yet
				if($switch == 0){
					#Print the line to OUTFILE1 and store it in an array of potentially private variants
					@temparray=();
					for($j=0;$j<$toks2length;$j++){
						$id=@toks2[$j];
						push @temparray, @{$genotypes{$id}}[$i];
						push @{$geno{$id}}, @{$genotypes{$id}}[$i];
					}
					print OUTFILE1 join(',', @temparray);
					print OUTFILE1 "\n";
					$prvline=join(',', @temparray);
					push @potentialrarevariants, $prvline;
					$switch = 1;
				#the genotype line has already been printed, don't print again
				}	
			}
		}
	}
}

#This is the third part of PBRV. 
###Lineage specific alleles will be shared on the same haplotype for long stretches of the chromosome. 
###The pattern of individuals carrying the allele is sufficient to identify the haplotype. 
###The purpose of this script is to apply this concept by giving each haplotype a name for each founder specific allele carried in each individual.

#make a list of the ids
@ids=@toks2[3..(scalar(@toks2)-1)];
#traverse through the list as you would going through a file line by line
for($i=0;$i<scalar(@{$geno{Marker}});$i++){
	#create an array that holds a list of the individuals that carry the lineage specific allele
	@currentgroup=();
	#add each individual carrying the allele to an array @currentgroup
	foreach $id (@ids){
		if(@{$geno{$id}}[$i] == 1){
			push @currentgroup, $id;
		}
	}
	#Sort the array @currentgroup so that the name of the haplotype will be consistent 
	@currentgroupsort = sort { lc($a) cmp lc($b) } @currentgroup;
	#join the array of names into a string that will be used as a unique identifier for this haplotype
	$token=join('v',@currentgroupsort);
	#For each person carrying the allele store the current haplotype (ie. $token) in the hash %state
	foreach $id (@ids){
		if(@{$geno{$id}}[$i] == 1){
			@{$state{$id}}[$i]=$token;
		}
	}
}


#This is the fourth part of PBRV.
#This section of code takes the list of haplotypes (from the previous step) each person carries and strings them together into segments


#Output the header line of the out file
push @segmentmatrixarray, "ID,PseudoFounder,Start,End,ChrHap\n";

foreach $id (@ids){
	#Blank a bunch of variables. This doesn't matter for the first file but it matters for subsequent files.
	@toks4=();
	@marker=();
	@Chr=();
	@location=();
	@state=();
	$state1=();
	$state1start=();
	$state1last=();
	$state2=();
	$state2start=();
	$state2last=();
	#Make arrays to hold the marker name, chr, position, and the current haplotype (ie state) for that marker, these arrays will be in the same order
	#In theory we don't have to do this, but its easier to follow what's going on this way
	for($i=0;$i<scalar(@{$geno{Marker}});$i++){
		if(@{$state{$id}}[$i]){
			push @marker, @{$geno{Marker}}[$i];
			push @Chr, @{$geno{Chr}}[$i];
			push @location, @{$geno{Position}}[$i];
			push @state, @{$state{$id}}[$i];
		}
	}
	#establish first two states
	#Set a switch for when both states have been established
	$established=0;
	#Set a counter that will be used to keep track of progression along the chromosome
	$i=0;
	#State 1 is the first state observed, easy peasy
	$state1=@state[$i];
	#And it is also the start of the first segment
	$state1start=@location[$i];
	#Currently it is the end of the segment because we haven't observed any other positions along the chromosome
	$state1last=@location[$i];
	#For as long as both states have not been established, try to find the second state
	while($established==0){
		$i++;
		#The second state will be different from the first state, so this is an easy check
		if(@state[$i] ne $state1){
			#We found the second state. Set the second state
			$state2=@state[$i];
			#Set the second state start point
			$state2start=@location[$i];
			#Set the second state last to the current position
			$state2last=@location[$i];
			#Both states are established, so change the switch so the program can move forward
			$established=1;
		}
	}

	#print STDOUT "State1=$state1\tState2=$state2\n";

	#the tricky part, identify the switches
	$i++;
	#Do this until you reach the end of the array of variants/positions
	while($i<scalar(@state)){
		#If a new state is observed it will not match one of the two current states
		if(@state[$i] ne $state1 && @state[$i] ne $state2){
			##count forward to double check. Call the 'testrealswitch function, which looks forward to determine if we continue to see this new state
			$realswitch=&testrealswitch($i);
			#If the 'testrealswitch' function determined that this is a real switch
			if($realswitch==1){
				#Call the 'findnext' function to determine which state changed
				$newresult=&findnext($i);
				#print STDOUT "$newresult\n";
				#If state 1 changed print the old segment to the outfile, then change the current state settting. If state 1 didn't change, state 2 changed. Do the same for that state.
				if($newresult == 1){
					push @segmentmatrixarray, "$id,$state1,$state1start,$state1last,1\n"; 
					$state1=@state[$i];
					$state1start=@location[$i];
					$state1last=@location[$i];
				}else{
					push @segmentmatrixarray, "$id,$state2,$state2start,$state2last,2\n"; 	
					$state2=@state[$i];
					$state2start=@location[$i];
					$state2last=@location[$i];
				}
			}
		}elsif(@state[$i] eq $state1){
			#The current position is the same as state one, change the last marker seen for state one to the current marker
			$state1last=@location[$i];
		}elsif(@state[$i] eq $state2){
			#The current position is the same as state two, change the last marker seen for state two to the current marker
			$state2last=@location[$i];
		}else{
			#This should never happen, but if it does print something crazy out so they will know something went wrong.
			print STDOUT "WTF? Contact the mothership!!\n";
		}
		$i++;
	}
	#print STDOUT "@state[$i]\n";
	#You have reached the end of the array of positions/states. Print the last two segments to the outfile.
	push @segmentmatrixarray, "$id,$state1,$state1start,$state1last,1\n";
	push @segmentmatrixarray, "$id,$state2,$state2start,$state2last,2\n"; 
}

#This section checks that the haplotype segments are placed on the correct chromosomes based on if they are shared with pedigree members on the probands maternal or paternal side of the pedigree

#First is makes a list of the founders for each proband's maternal or paternal side
foreach $lineage (@lineages){
	#Split the path into an array containing the individuals in this path
	@toks5 = split('_',$lineage);
	#Get the length of the array
	$currlength=scalar(@toks5);
	#figure out the current id
	$currid=@toks5[$currlength-1];
	#if the length is one that means that the current id is a founder
	if($currlength == 1){
		#If the current id is a founder use themselves as the maternal and paternal founder
		push @{$maternal{$currid}}, @toks5[0];
		push @{$paternal{$currid}}, @toks5[0];
	}else{
		#if the person isn't a founder figure out their parent on this lineage
		$currparent=@toks5[$currlength-2];
		#if their parent on this lineage is their mom add the current founder of this lineage to a hash of their maternal founders, and same is true of paternal side. %matid and %patid hashes were established at the beginning of the program
		if($currparent eq @{$matid{$currid}}[0]){
			push @{$maternal{$currid}}, @toks5[0];
		}elsif($currparent eq @{$patid{$currid}}[0]){
			push @{$paternal{$currid}}, @toks5[0];
		}
	}
}

foreach $id (@offspr){
	#clear the matseen hash
	my %matseen;
	#clear the patseen hash
	my %patseen;
	#for each  founder on each side (maternal and paternal) make sure they are unique, then sort them, and store them in a hash
	@{$matfounders{$id}} = grep {!$matseen{$_}++} @{$maternal{$id}};
	@{$patfounders{$id}} = grep {!$patseen{$_}++} @{$paternal{$id}};
}

#This section figures out if the individuals sharing the haplotype are on the maternal or paternal side of the probands relatives
#open(OUTFILE2,">Segment_matrix_2.csv")||die;
#print OUTFILE2 "@segmentmatrixarray[0]";
push @segmentmatrixarray2, @segmentmatrixarray[0];

for($i=1;$i<scalar(@segmentmatrixarray);$i++){
	$line2=@segmentmatrixarray[$i];
	chomp $line2;
	@toks2=split(",",$line2);
	$id=@toks2[0];
	@currmatlineage=();
	@currpatlineage=();
	@tempmatlineage=();
	@temppatlineage=();
	@currmatlineage=@{$matfounders{$id}};
	@currpatlineage=@{$patfounders{$id}};
	if("$id" eq "@{$matfounders{$id}}[0]"){
		push @segmentmatrixarray2, @segmentmatrixarray[$i];
		#print OUTFILE2 "$line2\n";
	}else{
		#get the list of people also sharing the haplotype
		@shared=split("v",@toks2[1]);
		foreach $id2 (@shared){
			if($id ne $id2){
				foreach $matlineage (@currmatlineage){
					foreach $matlineage2 (@{$matfounders{$id2}}){
						if($matlineage eq $matlineage2){
							push @tempmatlineage, $matlineage2;
						}
					}
					foreach $patlineage2 (@{$patfounders{$id2}}){
						if($matlineage eq $patlineage2){
							push @tempmatlineage, $patlineage2;
						}
					}
				}
				foreach $patlineage (@currpatlineage){
					foreach $matlineage2 (@{$matfounders{$id2}}){
						if($patlineage eq $matlineage2){
							push @temppatlineage, $matlineage2;
						}
					}
					foreach $patlineage2 (@{$patfounders{$id2}}){
						if($patlineage eq $patlineage2){
							push @temppatlineage, $patlineage2;
						}
					}
				}
				@currmatlineage=@tempmatlineage;
				@currpatlineage=@temppatlineage;
				@tempmatlineage=();
				@temppatlineage=();
			}
		}
		if(scalar(@currmatlineage) > 0 & scalar(@currpatlineage) == 0){
			push @segmentmatrixarray2, "@toks2[0],@toks2[1],@toks2[2],@toks2[3],1\n";
		}elsif(scalar(@currpatlineage) > 0 & scalar(@currmatlineage) == 0){
			push @segmentmatrixarray2, "@toks2[0],@toks2[1],@toks2[2],@toks2[3],2\n";
		}
	}
}
#close OUTFILE2;

for($i=1;$i<scalar(@segmentmatrixarray2);$i++){
	$line1=@segmentmatrixarray2[$i];
#foreach $line1 (@segmentmatrixarray2){
	#Remove the end of line
	chomp $line1;
	#Split the line into an array
	@toks1=split(",",$line1);
	$id = @toks1[0];
	#print STDOUT "$id\n";
	@{$idlisthash{$id}}[0]=1;
	$chrhap=@toks1[4];
	if($chrhap==1){
		push @{$pseudofounder1{$id}}, @toks1[1];
		push @{$start1{$id}}, @toks1[2];
		push @{$end1{$id}}, @toks1[3];
	}elsif($chrhap==2){
		push @{$pseudofounder2{$id}}, @toks1[1];
		push @{$start2{$id}}, @toks1[2];
		push @{$end2{$id}}, @toks1[3];
	}
}

for($i=0;$i<scalar(@{$genotypes{Position}});$i++){
	#print STDOUT "I=$i\n";
	$loc=@{$genotypes{Position}}[$i];
	@{$loclookup{$loc}}[0]=$i;
}

@idlist= keys %idlisthash;

foreach $id (@idlist){
	@startsarray=();
	@startsarray= sort{$a <=> $b} @{$start1{$id}};
	for($i=0;$i<scalar(@{$start1{$id}});$i++){
		$currtok=@{$start1{$id}}[$i];
		#print "$id:$currtok\n";
		for($j=0;$j<scalar(@startsarray);$j++){
			if($currtok == @startsarray[$j]){
				@{$arrayorder1{$id}}[$i]=$j;
				@{$arrayorderlookup1{$id}}[$j]=$i;
			}
		}
	}
	@startsarray=();
	@startsarray= sort {$a <=> $b} @{$start2{$id}};
	for($i=0;$i<scalar(@{$start2{$id}});$i++){
		$currtok=@{$start2{$id}}[$i];
		for($j=0;$j<scalar(@startsarray);$j++){
			if($currtok == @startsarray[$j]){
				@{$arrayorder2{$id}}[$i]=$j;
				@{$arrayorderlookup2{$id}}[$j]=$i;
			}
		}
	}
}

foreach $id (@idlist){
	#print STDOUT "$id:\n\t";
	for ($i=0;$i<scalar(@{$arrayorderlookup1{$id}});$i++){
		$currtokspot=@{$arrayorderlookup1{$id}}[$i];
		#Dont't forget to do this with both sides
		#print STDOUT "HAPSIDE 1:\n\t\tPseudofounder=@{$pseudofounder1{$id}}[$currtokspot]\n\t\tStart=@{$start1{$id}}[$currtokspot]\n\t\tEnd=@{$end1{$id}}[$currtokspot]\n";
		@sharedidlist=();
		@sharedidlist=split('v',@{$pseudofounder1{$id}}[$currtokspot]);
		$start=@{$start1{$id}}[$currtokspot];
		$startarrayposition=@{$loclookup{$start}}[0];
		$end=@{$end1{$id}}[$currtokspot];
		$endarrayposition=@{$loclookup{$end}}[0];
		$stop=0;
		$goingback=$startarrayposition;
		while($stop==0){
			#print STDOUT "Going back: $goingback\n";
			$zerosum=0;
			$twosum=0;
			foreach $id2 (@sharedidlist){		
				if(@{$genotypes{$id2}}[$goingback]==0){
					$zerosum++;
				}elsif(@{$genotypes{$id2}}[$goingback]==2){
					$twosum++;
				}
			}
			if($zerosum > 0 & $twosum > 0){
				#print STDOUT "Entered 1\n";
				$stop++;
			}elsif($goingback==0){
				$stop++;
			}else{
				#print STDOUT "Entered 2\n";
				$goingback=$goingback-1;
			}
		}
		$stop=0;
		$goingforward=$endarrayposition;
		while($stop==0){
			#print STDOUT "Going forward: $goingforward\n";
			$zerosum=0;
			$twosum=0;
			foreach $id2 (@sharedidlist){		
				if(@{$genotypes{$id2}}[$goingforward]==0){
					$zerosum++;
				}elsif(@{$genotypes{$id2}}[$goingforward]==2){
					$twosum++;
				}
			}
			if($zerosum > 0 & $twosum > 0){
				#print STDOUT "Entered 1\n";
				$stop++;
			}elsif($goingforward>= (scalar(@{$genotypes{$id}})-1)){
				$stop++;
			}else{
				#print STDOUT "Entered 2\n";
				$goingforward++;
			}
		}
		#print STDOUT "Old start: @{$start1{$id}}[$currtokspot]\nGoes back to: @{$genotypes{Position}}[$goingback]\nOld end: @{$end1{$id}}[$currtokspot]\nGoes forward to: @{$genotypes{Position}}[$goingforward]\n";
		@{$LRPstart1{$id}}[$currtokspot]=@{$genotypes{Position}}[$goingback];
		@{$LRPend1{$id}}[$currtokspot]=@{$genotypes{Position}}[$goingforward];
	}
	for ($i=0;$i<scalar(@{$arrayorderlookup2{$id}});$i++){
		$currtokspot=@{$arrayorderlookup2{$id}}[$i];
		#Dont't forget to do this with both sides
		#print STDOUT "HAPSIDE 1:\n\t\tPseudofounder=@{$pseudofounder1{$id}}[$currtokspot]\n\t\tStart=@{$start1{$id}}[$currtokspot]\n\t\tEnd=@{$end1{$id}}[$currtokspot]\n";
		@sharedidlist=();
		@sharedidlist=split('v',@{$pseudofounder2{$id}}[$currtokspot]);
		$start=@{$start2{$id}}[$currtokspot];
		$startarrayposition=@{$loclookup{$start}}[0];
		$end=@{$end2{$id}}[$currtokspot];
		$endarrayposition=@{$loclookup{$end}}[0];
		$stop=0;
		$goingback=$startarrayposition;
		while($stop==0){
			#print STDOUT "Going back: $goingback\n";
			$zerosum=0;
			$twosum=0;
			foreach $id2 (@sharedidlist){		
				if(@{$genotypes{$id2}}[$goingback]==0){
					$zerosum++;
				}elsif(@{$genotypes{$id2}}[$goingback]==2){
					$twosum++;
				}
			}
			if($zerosum > 0 & $twosum > 0){
				#print STDOUT "Entered 1\n";
				$stop++;
			}elsif($goingback==0){
				$stop++;
			}else{
				#print STDOUT "Entered 2\n";
				$goingback=$goingback-1;
			}
		}
		$stop=0;
		$goingforward=$endarrayposition;
		while($stop==0){
			#print STDOUT "Going forward: $goingforward\n";
			$zerosum=0;
			$twosum=0;
			foreach $id2 (@sharedidlist){		
				if(@{$genotypes{$id2}}[$goingforward]==0){
					$zerosum++;
				}elsif(@{$genotypes{$id2}}[$goingforward]==2){
					$twosum++;
				}
			}
			if($zerosum > 0 & $twosum > 0){
				#print STDOUT "Entered 1\n";
				$stop++;
			}elsif($goingforward>= (scalar(@{$genotypes{$id}})-1)){
				$stop++;
			}else{
				#print STDOUT "Entered 2\n";
				$goingforward++;
			}
		}
		#print STDOUT "Old start: @{$start1{$id}}[$currtokspot]\nGoes back to: @{$genotypes{Position}}[$goingback]\nOld end: @{$end1{$id}}[$currtokspot]\nGoes forward to: @{$genotypes{Position}}[$goingforward]\n";
		@{$LRPstart2{$id}}[$currtokspot]=@{$genotypes{Position}}[$goingback];
		@{$LRPend2{$id}}[$currtokspot]=@{$genotypes{Position}}[$goingforward];
	}
}

foreach $id (@idlist){
	#print STDOUT "$id:\n\t";
	for ($i=0;$i<scalar(@{$arrayorderlookup1{$id}});$i++){
		$currtokspot=@{$arrayorderlookup1{$id}}[$i];
		$nextstart=@{$start1{$id}}[$currtokspot+1];
		$lastend=@{$end1{$id}}[$currtokspot-1];
		$curstart=@{$start1{$id}}[$currtokspot];
		$curend=@{$end1{$id}}[$currtokspot];
		$LRPstart=@{$LRPstart1{$id}}[$currtokspot];
		$LRPend=@{$LRPend1{$id}}[$currtokspot];
		$lastLRPend=@{$LRPend1{$id}}[$currtokspot-1];
		$nextLRPstart=@{$LRPstart1{$id}}[$currtokspot+1];
		###Figure out new start
		if($i==0){
			@{$newstart1{$id}}[$currtokspot]=@{$LRPstart1{$id}}[$currtokspot];	
		}else{
			if($lastend > $LRPstart){
				#This LRPstart is unreliable
				if($lastLRPend > $curstart){
					#The last LRTend is unreliable
					#Both LRP extensions are unreliable, go off of original start
					@{$newstart1{$id}}[$currtokspot]=$curstart;
				}else{
					#Use last LRTend + 1
					@{$newstart1{$id}}[$currtokspot]=$lastLRPend+1;
				}
			}else{
				#This LRPstart is reliable
				if($lastLRPend > $curstart){
					#The last LRTend is unreliable
					#Use this LRTstart
					@{$newstart1{$id}}[$currtokspot]=$LRPstart;
				}else{
					#Both LRT extensions are reliablish
					if($lastLRPend > $LRPstart){
						#The LRP extensions overlap
						#Use last LRTend + 1
						@{$newstart1{$id}}[$currtokspot]=$lastLRPend+1;
					}else{
						#The LRP extensions don't overlap
						#Use this LRTstart
						@{$newstart1{$id}}[$currtokspot]=$LRPstart;
					}
				}
			}
		}
		if($i==(scalar(@{$arrayorderlookup1{$id}})-1)){
			@{$newend1{$id}}[$currtokspot]=@{$LRPend1{$id}}[$currtokspot];
		}else{
			if($nextstart < $LRPend){
				#This LRP end is unreliable
				if($nextLRPstart < $curend){
					#The next LRP start is unreliable
					#Both LRP extensions are unreliable, go off of the original end
					@{$newend1{$id}}[$currtokspot]=$curend;
				}else{
					#Use the next LRP start -1
					@{$newend1{$id}}[$currtokspot]=$nextLRPstart-1;
				}	
			}else{
				#This LRP end is reliable
				if($nextLRPstart < $curend){
					#The next LRP start is unreliable
					#Use the current LRP end
					@{$newend1{$id}}[$currtokspot]=$LRPend;
				}else{
					#Both LRP extensions are reliable
					if($nextLRPstart < $LRPend){
						#The LRP extensions overlap
						#Use the next LRP start -1
						@{$newend1{$id}}[$currtokspot]=$nextLRPstart-1;
					}else{
						#The LRP extensions don't overlap
						#Use this LRP end
						@{$newend1{$id}}[$currtokspot]=$LRPend;
					}
				}
			}
		}
		#Dont't forget to do this with both sides
		#print STDOUT "HAPSIDE 1:\n\t\tPseudofounder=@{$pseudofounder1{$id}}[$currtokspot]\n\t\tStart=@{$start1{$id}}[$currtokspot]\n\t\tEnd=@{$end1{$id}}[$currtokspot]\n";
#		@sharedidlist=();
#		@sharedidlist=split('v',@{$pseudofounder1{$id}}[$currtokspot]);
#		$correctsum=();
#		$incorrectsum=();
#		foreach $id2 (@sharedidlist){
#			for($k=0;$k<scalar(@{$pseudofounder1{$id2}});$k++){
#				if("@{$pseudofounder1{$id2}}[$k]" eq "@{$pseudofounder1{$id}}[$currtokspot]"){
#					if(@{$start1{$id2}}[$k] < @{$end1{$id}}[$currtokspot] & @{$start1{$id}}[$currtokspot] < @{$end1{$id2}}[$k]){
#						print STDOUT "Current ID:$id2 HAPSIDE 1:\n\t\tPseudofounder=@{$pseudofounder1{$id2}}[$k]\n\t\tStart=@{$start1{$id2}}[$k]\n\t\tEnd=@{$end1{$id2}}[$k]\n";
#					}
#				}
#			}
#			for($k=0;$k<scalar(@{$pseudofounder2{$id2}});$k++){
#				if("@{$pseudofounder2{$id2}}[$k]" eq "@{$pseudofounder1{$id}}[$currtokspot]"){
#					if(@{$start2{$id2}}[$k] < @{$end1{$id}}[$currtokspot] & @{$start1{$id}}[$currtokspot] < @{$end2{$id2}}[$k]){
#						print STDOUT "Current ID:$id2 HAPSIDE 2:\n\t\tPseudofounder=@{$pseudofounder2{$id2}}[$k]\n\t\tStart=@{$start2{$id2}}[$k]\n\t\tEnd=@{$end2{$id2}}[$k]\n";
#					}
#				}
#			}
#		}
		
		
			
#		for($j=$startarrayposition;$j<($endarrayposition+1);$j++){
#			$zerosum=0;
#			$twosum=0;
#			foreach $id2 (@sharedidlist){
#				if(@{$genotypes{$id2}}[$j]==0){
#					$zerosum++;
#				}elsif(@{$genotypes{$id2}}[$j]==2){
#					$twosum++;
#				}
#			}
#			if($zerosum > 0 & $twosum > 0){
#				$incorrectsum++;
#			}else{
#				$correctsum++;
#			}
#		}
#		$percentcorrect=$correctsum/($correctsum+$incorrectsum);
#		print STDOUT "Percent matching=$percentcorrect\n\n";
	}
	for ($i=0;$i<scalar(@{$arrayorderlookup2{$id}});$i++){
		$currtokspot=@{$arrayorderlookup2{$id}}[$i];
		$nextstart=@{$start2{$id}}[$currtokspot+1];
		$lastend=@{$end2{$id}}[$currtokspot-1];
		$curstart=@{$start2{$id}}[$currtokspot];
		$curend=@{$end2{$id}}[$currtokspot];
		$LRPstart=@{$LRPstart2{$id}}[$currtokspot];
		$LRPend=@{$LRPend2{$id}}[$currtokspot];
		$lastLRPend=@{$LRPend2{$id}}[$currtokspot-1];
		$nextLRPstart=@{$LRPstart2{$id}}[$currtokspot+1];
		###Figure out new start
		if($i==0){
			@{$newstart2{$id}}[$currtokspot]=@{$LRPstart2{$id}}[$currtokspot];	
		}else{
			if($lastend > $LRPstart){
				#This LRPstart is unreliable
				if($lastLRPend > $curstart){
					#The last LRTend is unreliable
					#Both LRP extensions are unreliable, go off of original start
					@{$newstart2{$id}}[$currtokspot]=$curstart;
				}else{
					#Use last LRTend + 1
					@{$newstart2{$id}}[$currtokspot]=$lastLRPend+1;
				}
			}else{
				#This LRPstart is reliable
				if($lastLRPend > $curstart){
					#The last LRTend is unreliable
					#Use this LRTstart
					@{$newstart2{$id}}[$currtokspot]=$LRPstart;
				}else{
					#Both LRT extensions are reliablish
					if($lastLRPend > $LRPstart){
						#The LRP extensions overlap
						#Use last LRTend + 1
						@{$newstart2{$id}}[$currtokspot]=$lastLRPend+1;
					}else{
						#The LRP extensions don't overlap
						#Use this LRTstart
						@{$newstart2{$id}}[$currtokspot]=$LRPstart;
					}
				}
			}
		}
		if($i==(scalar(@{$arrayorderlookup2{$id}})-1)){
			@{$newend2{$id}}[$currtokspot]=@{$LRPend2{$id}}[$currtokspot];
		}else{
			if($nextstart < $LRPend){
				#This LRP end is unreliable
				if($nextLRPstart < $curend){
					#The next LRP start is unreliable
					#Both LRP extensions are unreliable, go off of the original end
					@{$newend2{$id}}[$currtokspot]=$curend;
				}else{
					#Use the next LRP start -1
					@{$newend2{$id}}[$currtokspot]=$nextLRPstart-1;
				}
			}else{
				#This LRP end is reliable
				if($nextLRPstart < $curend){
					#The next LRP start is unreliable
					#Use the current LRP end
					@{$newend2{$id}}[$currtokspot]=$LRPend;
				}else{
					#Both LRP extensions are reliable
					if($nextLRPstart < $LRPend){
						#The LRP extensions overlap
						#Use the next LRP start -1
						@{$newend2{$id}}[$currtokspot]=$nextLRPstart-1;
					}else{
						#The LRP extensions don't overlap
						#Use this LRP end
						@{$newend2{$id}}[$currtokspot]=$LRPend;
					}
				}
			}		
		}
	}
}

push @segmentmatrixarray3, "$header1";
foreach $id (@idlist){
	#print STDOUT "$id:\n";
	for ($i=0;$i<scalar(@{$arrayorderlookup1{$id}});$i++){
		$currtokspot=@{$arrayorderlookup1{$id}}[$i];
		#print STDOUT "Old: @{$start1{$id}}[$currtokspot] - @{$end1{$id}}[$currtokspot]\nNew: @{$newstart1{$id}}[$currtokspot] - @{$newend1{$id}}[$currtokspot]\n" ;
		push @segmentmatrixarray3, "$id,@{$pseudofounder1{$id}}[$currtokspot],@{$newstart1{$id}}[$currtokspot],@{$newend1{$id}}[$currtokspot],1\n";
	}
	#print STDOUT "\n";	
	for ($i=0;$i<scalar(@{$arrayorderlookup2{$id}});$i++){
		$currtokspot=@{$arrayorderlookup2{$id}}[$i];
		#print STDOUT "Old: @{$start2{$id}}[$currtokspot] - @{$end2{$id}}[$currtokspot]\nNew: @{$newstart2{$id}}[$currtokspot] - @{$newend2{$id}}[$currtokspot]\n" ;
		push @segmentmatrixarray3, "$id,@{$pseudofounder2{$id}}[$currtokspot],@{$newstart2{$id}}[$currtokspot],@{$newend2{$id}}[$currtokspot],2\n";
	}
	#print STDOUT "\n";
}

for($i=1;$i<scalar(@segmentmatrixarray3);$i++){
	$line1=@segmentmatrixarray3[$i];
	#Remove the end of line
	chomp $line1;
	#Split the line into an array
	@toks1=split(",",$line1);
	#Create a string that is composed of the haplotypename, start, and end. This string will be used as a key for the hash %hashing
	$holder=@toks1[1]."@".@toks1[2]."@".@toks1[3];
	#Create a variable that holds the ID of the individual for this line
	$id=@toks1[0];
	#Create a key of the string '$holder' created above. Put a 1 in the 0th position as a holder.
	@{$hashing{$holder}}[0]=1;
	#Add the string 'holder' to the end of an array in the hash %hashing2 in where the current individual is the key.
	push @{$hashing2{$id}}, $holder;
}

#This is another step in the PBRV program. This step consolidates redundant haplotypes identified in multiple people into single records. 
#Establish a switch that is used to reset a count of the current array position
#Open the outfile for the condensed format
$switcheroosky=0;
#Create an array (@haps) of the unique haplotypes observed (in this case unique means the same people and same start and end positions
@haps=keys %hashing;
#Find out how many unique haplotypes there are
$originallength=scalar(@haps);
#Go through the list of unique haplotypes
for($i=0;$i<scalar(@haps);$i++){
	#If the switch is flipped reset the current position to 0 to start over. Also, reset the switch.
	if($switcheroosky == 1){
		$i=0;
		$switcheroosky=0;
	}
	#print STDOUT "I=$i\n";
	#Split the current haplotype string and establish variables to hold the hap name, start, and end positions
	$unique1=@haps[$i];
	@toks1=split("@",$unique1);
	$hap1=@toks1[0];
	$start1=@toks1[1];
	$end1=$toks1[2];
	#Go through the list of unique haplotypes to compare to the outer current haplotype
	for($j=0;$j<scalar(@haps);$j++){
		#If $i == $j this is a special case. We don't need to compare something to itself so just move on.
		if($i==$j){
			next;
		}else{
			#Split the current haplotype string and establish variables to hold the hap name, start, and end positions
			$unique2=@haps[$j];
			@toks2=split("@",$unique2);
			$hap2=@toks2[0];
			$start2=@toks2[1];
			$end2=$toks2[2];
			#If the two haplotype names are the same
			if($hap1 eq $hap2){
				#If the two haplotypes overlap
				if($start1 <= $end2 && $start2 <= $end1){
					#Remove the second haplotype
					splice @haps, $j, 1;
					#Set the new breakpoints to be the outer breakpoints
					if($start1 <= $start2){
						$start3=$start1;
					}else{
						$start3=$start2;
					}
					if($end1 >= $end2){
						$end3=$end1;
					}else{
						$end3=$end2;
					}
					#Create a new string with the haplotype name and haplotype breakpoints
					$holder=$hap1."@".$start3."@".$end3;
					#Insert this new string where the first haplotype was
					splice @haps, $i, 1, $holder;
					#Flip the switch to start over at 0
					$switcheroosky=1;
					#print STDOUT "J=$j\nFound a match\n";
					last;
				}
			}
		}
	}
}

##Now that @haps has been cleaned up, find any individuals that are obligated to carry each haplotype given the individuals known to carry the haplotype
foreach $unique (@haps){
	#Split each case into haplotype name, start, and end position
	@toks2=split("@",$unique);
	#The hap
	$hap=@toks2[0];
	#Get the ids for the people carrying this unique hap
	@uniqueids=split("v",$hap);
	#Clear the hash %hectors, which will be used to store ids for obligatory carriers
	%hectors=();
	#for each person carrying the current haplotype
	for($l=0;$l<(scalar(@uniqueids)-1);$l++){
		#for each person after the outter loop until the end
		for($m=$l+1;$m<scalar(@uniqueids);$m++){
			#clear the scalar variable $obligation
			$obligation = ();
			#clear the array @tooters
			@tooters=();
			#call the subroutine 'findobligatories' 
			$obligation = &findobligatories(@uniqueids[$l],@uniqueids[$m]);
			#Split the output from the subroutine 'findobligatories'
			@tooters=split("_",$obligation);
			#For each person in the array of obligatory carriers, hash them as keys in the hash %hectors
			foreach $snippet (@tooters){
				@{$hectors{$snippet}}[0]++;
			}
		}
	}
	#The people known to carry the current haplotype are also obviously carriers that need to be hashed as keys in the %hectors hash
	foreach $id (@uniqueids){
		@{$hectors{$id}}[0]++;
	}
	#Get the unique ids that are carries including the known and obligatory carriers
	@obligatorilyin=keys %hectors;
	#Create a string that can later be parsed with a "_" delimiter that has a list of all of the carriers of the current haplotype
	$all_ids= join("_",@obligatorilyin);
	#Print the current haplotype, start and end locations, and the list of ids including obligatory carriers
	#print OUTFILE3 "@toks2[0],@toks2[1],@toks2[2],$all_ids\n";
	#hash the current hap for future use in the program, also store the starts and ends in the first locations
	@{$haplist{$hap}}[0]++;
	push @{$starts{$hap}}, @toks2[1];
	push @{$ends{$hap}}, @toks2[2];
}
#Done, close the outfile
#close OUTFILE3;
		
#Get a list of unique haplotypes. In this case, uniqueness is not dependent on start and end location, but rather only on the individuals carrying the haplotype.
@uniquehaplist = keys %haplist;


#Impute genotypes into haplotypes
#Go through each unique haplotype, in this case it is defined just by the individuals carrying the haplotype
foreach $hap (@uniquehaplist){
	#print STDOUT "Working on Hap:$hap\n";
	#Create/blank a variable to hold ids
	@uniqueids=();
	#Get a list of the genotyped people that carry this haplotype
	@uniqueids=split("v",$hap);
	#foreach $person (@uniqueids) {
	#	print STDOUT "$person\n";
	#}
	#Iterate through all of the locations that this hap was observed
	for($i=0;$i<scalar(@{$starts{$hap}});$i++){
		#Get the start and end for where the hap was observed
		$start=@{$starts{$hap}}[$i];
		$end=@{$ends{$hap}}[$i];
		#For each marker/variant
		for($j=0;$j<scalar(@{$genotypes{Position}});$j++){
			#If the variant falls within the current hap 
			if(@{$genotypes{Position}}[$j] >= $start && @{$genotypes{Position}}[$j] <= $end){
				$zerosum=0;
				$onesum=0;
				#For each person carrying the current hap
				foreach $person (@uniqueids) {
					#If the person is homozygous then the allele on the haplotype is unambiguous
					if(@{$genotypes{$person}}[$j] eq '0'){
						$zerosum++;
						#@{$haploidgenotypes{$hap}}[$j]=0;
					#Ambiguous
					}elsif(@{$genotypes{$person}}[$j] eq '1'){
						next;
					#If the person is homozygous then the allele on the haplotype is unambiguous
					}elsif(@{$genotypes{$person}}[$j] eq '2'){
						#@{$haploidgenotypes{$hap}}[$j]=1;
						$onesum++;
					#This shouldn't happen, but just in case
					}else{
						next;
					}
				}
				########The problem is right here!!!!! Found that bitch!!!!
				$allelevalue=();
				if(($onesum+$zerosum)>0){
					$allelevalue=$onesum/($onesum+$zerosum);
					if($allelevalue>0.5){
						@{$haploidgenotypes{$hap}}[$j]=1;
					}elsif($allelevalue<0.5){
						@{$haploidgenotypes{$hap}}[$j]=0;
					}
				}
			}
		}
	}
}

#Iteratively loop around pedigree until no more updating can be accomplished

#For each marker/variant in the file
for($j=0;$j<scalar(@{$genotypes{Position}});$j++){
	#Establish/blank 3 variables that are need to determine if the program is finished updating
	$looptestscalarold=();
	$looptestscalar=();
	@looptest=();
	#Get the allele that each haplotype carries and make an array
	foreach $hap (@uniquehaplist){
		push @looptest, @{$haploidgenotypes{$hap}}[$j];
	}
	#Turn the array into a string. This will be used to check if updating is complete by checking it against $looptestscalarold which contains the same string from the previous iteration
	$looptestscalar=join(",",@looptest);
	#while loop exits when no more updating possible. This statement checks if any genotypes changed during the last iteration. If no change was made, then it is finished.
	while($looptestscalar ne $looptestscalarold){
	#while($looptestscalar ne $looptestscalarold){
		#print STDOUT "New:$looptestscalar,Old:$looptestscalarold\n";
		#Set $looptestscalarold, which is used to check if $looptestscalar is updated before the next iteration.
		$looptestscalarold = $looptestscalar;
		#For each hap
		foreach $hap (@uniquehaplist){
			#print STDOUT "Analyzing $hap.\n";
			#Make a list of ids carrying the hap at this location
			@uniqueids=();
			@uniqueids=split("v",$hap);	
			#For each segment that this hap is observed
			for($i=0;$i<scalar(@{$starts{$hap}});$i++){
				#The start and end of this segment
				$start=@{$starts{$hap}}[$i];
				$end=@{$ends{$hap}}[$i];
				#If the current location falls within the segment
				if(@{$genotypes{Position}}[$j] >= $start && @{$genotypes{Position}}[$j] <= $end){
					#print STDOUT "---@{$haploidgenotypes{$hap}}[$j]---\n";
					#If the allele the hap carries hasn't been established
					if(@{$haploidgenotypes{$hap}}[$j] ne '0' & @{$haploidgenotypes{$hap}}[$j] ne '1'){
						#print STDOUT "Entered\n";
						#Go through each person carrying the hap
						foreach $person (@uniqueids) {
							#print STDOUT "Searching ID:$person\n";
							#If the person is a heterozygote
							if(@{$genotypes{$person}}[$j] eq '1'){
								#Go through the list of the haplotypes the person carries
								foreach $item (@{$hashing2{$person}}){
									#Split this string into the haplotype, start, and end
									@itemtoks=split("@",$item);
									#This is the other haplotype
									$otherhap=@itemtoks[0];
									#print STDOUT "Other haplotype:$otherhap\t@{$haploidgenotypes{$otherhap}}[$j]\n";
									#Set a value that $othergeno can not acutally be. This is just a precaution. It will get blanked in the next line of code if the right data isn't available anyway.
									$othergeno=2;
									#Establish othergeno as the genotype of the complementary haplotype. If the haplotype isn't observed at this location it will just blank the scalar $othergeno.. 
									$othergeno=@{$haploidgenotypes{$otherhap}}[$j];
									#If the two haplotypes are not accidently the same one. Ie. ignore it if they are the same one.
									if($otherhap ne $hap ){
										#If the allele for the other haplotype is known
										if($othergeno eq '0' || $othergeno eq '1'){
											#print STDOUT "Checking other haplotype\n";
											#If it is 0, that means that the current haplotype carries the 1 allele
											if(@{$haploidgenotypes{$otherhap}}[$j] eq '0'){
												@{$haploidgenotypes{$hap}}[$j] = 1;
												#print STDOUT "Made an addition!\n";
											#If it is 1, that means that the current haplotype carries the 0 allele
											}elsif(@{$haploidgenotypes{$otherhap}}[$j] eq '1'){
												@{$haploidgenotypes{$hap}}[$j] = 0;
												#print STDOUT "Made an addition!\n";
											#This should never happen. Print an error statement to the command line.
											}else{
												print STDOUT "WTF BATMAN!!\n";
											}
										}	
									}
								}
							}
						}
					}
				}
			}
		}
		#As many alleles have been added to each haplotype as possible this round. Create a string again, to check to see if anything changed. Reminder: If nothing changed it will exit the loop.
		@looptest=();
		foreach $hap (@uniquehaplist){
			push @looptest, @{$haploidgenotypes{$hap}}[$j];
		}
		$looptestscalar=join(",",@looptest);
	}
	#print STDOUT "\n";
}

#Open the outfile
#open(OUTFILE4,">Haplotypes.csv")||die;
#Print the header line
push @tempholderarray, "Marker,Chr,Position";
#Continuation of printing the header line, print each haplotype as a column header
foreach $hap (@uniquehaplist){
	push @tempholderarray, $hap;
	#print OUTFILE4 ",$hap";	
}

$tempholder=join(",",@tempholderarray);
push @haplotypesmatrixarray, "$tempholder\n";
#Print end of line
#print OUTFILE4 "\n";
#For each marker
for($i=0;$i<scalar(@{$genotypes{Position}});$i++){
	@tempholderarray = ();
	#Print the name of the marker, the chromosome, and the position
	push @tempholderarray, "@{$genotypes{Marker}}[$i],@{$genotypes{Chr}}[$i],@{$genotypes{Position}}[$i]";
	#For each haplotype print the allele that the haplotype carries. This will just be blank where segments of the haplotype are not observed.
	foreach $hap (@uniquehaplist){
		push @tempholderarray, @{$haploidgenotypes{$hap}}[$i];
	}
	#print end of line
	$tempholder=join(",",@tempholderarray);
	push @haplotypesmatrixarray, "$tempholder\n";
	
	
}
#close the output file
#close OUTFILE4;

#open(INFILE1,"Haplotypes.csv")||die;

$headerline=@haplotypesmatrixarray[0];
@headertoks=split(',',$headerline);
$j=0;
#while($line1=<INFILE1>){
for($i=1;$i<scalar(@haplotypesmatrixarray);$i++){
	$line1=@haplotypesmatrixarray[$i];
	chomp $line1;
	my @toks1;
	@toks1=split(",",$line1);
	for($k=0;$k<scalar(@headertoks);$k++){
		$currcol=@headertoks[$k];
		@{$haplotypes{$currcol}}[$j]=@toks1[$k];
	}
	$j++;
}
#close INFILE1;

#open(INFILE2,"Segment_matrix_3.csv")||die;
$headerline2=@segmentmatrixarray3[0];
for($i=1;$i<scalar(@segmentmatrixarray3);$i++){
	$line2=@segmentmatrixarray3[$i];
#while($line2=<INFILE2>){
	chomp $line2;
	@toks2=split(",",$line2);
	$chromohap="@toks2[0]"."_"."@toks2[4]";
	push @{$pseudofounder{$chromohap}}, @toks2[1];
	push @{$start{$chromohap}}, @toks2[2];
	push @{$end{$chromohap}}, @toks2[3];
}
#close INFILE2;

@uniquechromohaps= keys %pseudofounder;
foreach $chromohap (@uniquechromohaps){
	for($i=0;$i<scalar(@{$pseudofounder{$chromohap}});$i++){
		for($j=0;$j<scalar(@{$haplotypes{'Position'}});$j++){
			if(@{$haplotypes{'Position'}}[$j]>=@{$start{$chromohap}}[$i] && @{$haplotypes{'Position'}}[$j]<=@{$end{$chromohap}}[$i]){
				$currentpseudohap=@{$pseudofounder{$chromohap}}[$i];
				@{$alleles{$chromohap}}[$j]=@{$haplotypes{$currentpseudohap}}[$j];
			}elsif(@{$haplotypes{'Position'}}[$j]>@{$end{$chromohap}}[$i]){
				last;
			}
		}
	}
}

open(INFILE3, $ARGV[1])||die;
$headerline=<INFILE3>;
chomp $headerline;
@headertoks=split(",",$headerline);
$j=0;
while($line=<INFILE3>){
	chomp $line;
	@toks=split(",",$line);
	for($i=3;$i<scalar(@toks);$i++){
		$chromohap1="@headertoks[$i]"."_1";
		$chromohap2="@headertoks[$i]"."_2";
		if(@{$alleles{$chromohap1}}[$j] ne "0" & @{$alleles{$chromohap1}}[$j] ne "1"){
			if(@{$alleles{$chromohap2}}[$j] eq "0" || @{$alleles{$chromohap2}}[$j] eq "1"){
				@{$alleles{$chromohap1}}[$j] = @toks[$i] - @{$alleles{$chromohap2}}[$j];
			}else{
				if(@toks[$i] eq '0'){
					@{$alleles{$chromohap1}}[$j]=0;
					@{$alleles{$chromohap2}}[$j]=0;
				}elsif(@toks[$i] eq '2'){
					@{$alleles{$chromohap1}}[$j]=1;
					@{$alleles{$chromohap2}}[$j]=1;
				}
			}
		}elsif(@{$alleles{$chromohap2}}[$j] ne "0" & @{$alleles{$chromohap2}}[$j] ne "1"){
			if(@{$alleles{$chromohap1}}[$j] eq "0" || @{$alleles{$chromohap1}}[$j] eq "1"){
				@{$alleles{$chromohap2}}[$j] = @toks[$i] - @{$alleles{$chromohap1}}[$j];
			}
		}
	}
	$j++;
}

####################
#Problem is right here! Use a different source to come up with the list of haplotypes to print, cool thanks.
####################

@haplisttoprint=keys %alleles;

open(OUTFILE,">Phased_chromosomes.csv")||die;
print OUTFILE "Marker,Chr,Position";
foreach $chromohap (@haplisttoprint){
	print OUTFILE ",$chromohap";
}
print OUTFILE "\n";
for($j=0;$j<scalar(@{$haplotypes{'Position'}});$j++){
	print OUTFILE "@{$haplotypes{'Marker'}}[$j],@{$haplotypes{'Chr'}}[$j],@{$haplotypes{'Position'}}[$j]";
	foreach $chromohap (@haplisttoprint){
		print OUTFILE ",@{$alleles{$chromohap}}[$j]";
	}
	print OUTFILE "\n";
}
close OUTFILE;





######Subroutines######
#Recursion algorithm that traverses the pedigree and finds each individuals offspring (this is called starting with the founder). For each kid it adds to the array (@currfoundoffpring) and count ($thisfoundercount).
sub find_paths{
	$thisfoundercount++;
	#The current person to search for their offspring, must use 'my' to declare that each time this function is called that the variable is unique to make the recursion work 
	my $thisperson = $_[0];
	#print the current ID to outfile
	#print OUTFILE "$_[1]\n";
	push @lineages, $_[1];
	#Search through the arrays of moms and dads
	for(my $i = 0; $i < $numberofids;$i++){
		#If the current person is a mom or a dad in the pedigree
		if (@paternal[$i] eq $thisperson || @maternal[$i] eq $thisperson){
			#Add the person to the string contained in $path 'my' has to be called for recursion to work
			my $path = $_[1]."_".@offspr[$i];
			#Add the person to their offspring to the array of offspring
			push @currfoundoffpring, @offspr[$i];
			#Call the algorithm again using the offspring and sending the new path variable
			&find_paths(@offspr[$i],$path);
		}
	}	
}

sub testrealswitch {

	#print STDOUT "Entered new test\n";
	#Set two counters to the location passed to the subroutine in @_[0]
	$j=@_[0];
	$k=@_[0];
	#set some switches and counts
	#This is the count of the number of times the new state has been observed 
	$idpositive=0;
	#A switch that counts if state 1 has been observed
	$state1positive=0;
	#A switch that counts if state 2 has been observed
	$state2positive=0;
	#A switch that counts if a third state has been observed. This will occasionally happen when two recombination events are near each other in two converging lineages.
	$state3positive=0;
	#Start looking forward through the next positions. Depending on results this can be a variable number of spots so I use that 'last' function to exit when it's appropriate.
	while($j<scalar(@state)){
		#print STDOUT "####\nNEW MARKER\n####\n";
		$j++;
		#The new state is observed again, add to the counter
		if(@state[$j] eq @state[$k]){
			$idpositive++;
		}else{
			#The new state was not observed this time. It is either one of three other options:
			#print STDOUT "@state[$k]\t@state[$j]\t$state1\t$state2\n";
			if(@state[$j] eq $state1){
				#Option 1: It is the same as state 1. Flip the state 1 switch.
				#print STDOUT "State 1 positive\n";
				$state1positive=1;
			}elsif(@state[$j] eq $state2){
				#Option 2: It is the same as state 2. Flip the state 2 switch.
				#print STDOUT "State 2 positive\n";
				$state2positive=1;
			}else{
				#Option 3: It is a third state, and there has potentially been another recombination event affecting the other chromosome. I don't really do anything in this case, but if you are having issues this could be helpful to put a print statement here.
				#$state3positive=1;
				#print STDOUT "Two recombinations observed very close to each other. In some cases this may cause a problem. Carefully inspect the output from this software.\n";
			}
		}
		if($state1positive==1 && $state2positive==1){
			#Both of the old states have been observed before the necessary number of observations of the new state. Return '0', meaning that this is not a real switch.
			return 0;
			#Exit the loop
			last;
		}elsif($idpositive==@ARGV[2]){
			#We have observed the new state a sufficient number of times to be convinced it is real. Return 1 indicating it is real and exit the loop.
			#print "Entered\n";
			return 1;
			last;
		}
	}	
}

sub findnext {
	#@_[0] is the current location
	#Set two counters to the location passed to the subroutine in @_[0]
	$j=@_[0];
	$k=@_[0];
	$altreally=0;
	#Start looking forward through the next positions. Depending on results this can be a variable number of spots so I use that 'last' function to exit when it's appropriate.
	while($j<scalar(@state)){
		$j++;
		#Iterate through until you observe the haplotype that didn't change
		if(@state[$j] ne @state[$k]){
			#print "Entered\n";
			if(@state[$j] eq $state1){
				#If the haplotype that didn't change is 1, return 2 indicating it was the second haplotype that changed
				return 2;
				last;
			}elsif(@state[$j] eq $state2){
				#If the haplotype that didn't change is 2, return 1 indicating it was the second haplotype that changed
				return 1;
				last;
			}else{
				#There is a third state which happens to be real
				if($alt=@state[$k]){
					$altreally++;
				}else{
					$alt=@state[$k];
					$altreally=1;
				}
				if($altreally == @ARGV[2]){
					return 1;
					last;
				}
			}
		}
	}
}

#This is a subroutine that given two ids that are carriers of a founder specific allele identifies individuals that are obligatory carriers of the allele as well
sub findobligatories {
	#Establish the ids
	$id1 = @_[0];
	$id2 = @_[1];
	#Establish a switch to speed up the program. It skips segments if the switch is flipped to 1.
	$done=0;
	#Blank some variables
	@returntoks=();
	$returntoks=();
	#For each lineage that id 1 is the last person in
	foreach $lineage (@{$last{$id1}}){
		@lineagetoks2=();
		#Get a list of ids in this lineage, in order
		@lineagetoks2=split("_",$lineageline);
		#Go through each person, and if they are id2, then the people between them are obligatory carries. Establish the return variables and set the switch $done to 1 to skip forward in the program.
		for($i=0;$i<scalar(@lineagetoks2);$i++){
			#For each person in the lineage, if they are $id2
			if(@lineagetoks2[$i] eq $id2){
				#Return the ids between them including them
				@returntoks=@lineagetoks2[$i..(scalar(@lineagetoks2)-1)];
				#Create a string of the array for a simpler return value
				$returntoks=join("_",@returntoks);
				#Change the switch, since the other two scenarios don't need to be run
				$done=1;
			}
		}
	}
	#This section uses the same rationale as the previous section, but in reverse order
	if($done == 0){
		foreach $lineage (@{$last{$id2}}){
			@lineagetoks2=();
			@lineagetoks2=split("_",$lineageline);
			for($i=0;$i<scalar(@lineagetoks2);$i++){
				if(@lineagetoks2[$i] eq $id1){
					@returntoks=@lineagetoks2[$i..(scalar(@lineagetoks2)-1)];
					$returntoks=join("_",@returntoks);
					$done=1;
				}
			}
		}
	}
	#This section identifies all of the potential founders that the two individuals share. If a person is seen in all of these lineages then they are an obligatory carrier. Ex. If half-siblings with the same mom share a potential lineage specific variant, then their mother is an obligatory carrier.
	if($done == 0){
		@founderstogether=();
		#make a list of all of the founders that are shared between the two individuals
		foreach $founder1 (@{$founderlist{$id1}}){
			foreach $founder2 (@{$founderlist{$id2}}){
				if($founder1 eq $founder2){
					push @founderstogether, $founder1;
				}
			}
		}
		%counthash1=();
		%counthash2=();
		%returnhash=();
		$potentialcount=0;
		#go through the list of shared founders 
		foreach $olddude (@founderstogether){
			#Make the two lineage hashes that store the appropriate lineages for the two id's and the current founder
			$combo1=$olddude.$id1;
			$combo2=$olddude.$id2;
			#Establish the lineages by referencing the hash
			$lineageline1=@{$combohash{$combo1}}[0];
			$lineageline2=@{$combohash{$combo2}}[0];
			#Split the lineages
			@lineagelinetoks1=split("_",$lineageline1);
			@lineagelinetoks2=split("_",$lineageline2);
			#Make a count for each individual in each lineage
			foreach $token1 (@lineagelinetoks1){
				@{$counthash1{$token1}}[0]++;
			}
			foreach $token2 (@lineagelinetoks2){
				@{$counthash2{$token2}}[0]++;
			}
			#Count the number of potential founders
			$potentialcount++;
		}
		#For each person in one of these lineages
		foreach $person1 (keys %counthash1){
			#If they are observed in every lineage from themselves to the founders shared between these two individuals
			if(@{$counthash1{$person1}}[0] == $potentialcount){
				#Flip a switch that this individual should be returned
				@{$returnhash{$person1}}[0]++;
			}
		}
		foreach $person2 (keys %counthash2){
			#If they are observed in every lineage from themselves to the founders shared between these two individuals
			if(@{$counthash2{$person2}}[0] == $potentialcount){
				#Flip a switch that this individual should be returned
				@{$returnhash{$person2}}[0]++;
			}
		}
		#Get the ids from the hash of switches 
		@returntoks=keys %returnhash;
		#Create a string if ids to return
		$returntoks2=join("_",@returntoks);
	}
	return $returntoks2;
}

