#!/usr/local/bin/perl
srand(localtime);

open (INFILE_PED, @ARGV[0]) || die;
$header = <INFILE_PED>;
#Identify founders
while ($trioline = <INFILE_PED>){
	chomp $trioline;
	@triotoks = split(',', $trioline);
	push @maternal, @triotoks[2];
	push @paternal, @triotoks[1];
	push @offspr, @triotoks[0];
	push @sex, @triotoks[3];
	push @sequenced, @triotoks[6];
	if(@triotoks[2] eq '0'){
		push @founders, @triotoks[0];
	}
}
for ($i=0;$i<scalar(@offspr);$i++){
	$id=@offspr[$i];
	push @{$idsexhash{$id}}, @sex[$i];
}
close INFILE_PED;

####Drop data through pedigree####
####Initiate founder phase####
foreach $founder (@founders){
		$nameone="$founder"."_1";
		$nametwo="$founder"."_2";
		@{$phase1{$founder}}[0] = $nameone;
		@{$phase2{$founder}}[0] = $nametwo;
		#print STDOUT "@{$phase1{$founder}}[0],@{$phase2{$founder}}[0]\n";
}

#Initiate dropping
$numberofids=scalar(@offspr);
##@ready holds ids that are ready to be dropped


for($q=0;$q<@ARGV[1];$q++){
#print STDOUT "Count:$q\n";
%ready1=();
%ready2=();
@ready=@founders;
foreach $id (@ready){
	@currfoundoffpring=();
	&find_kids($id);
	for $kid (@currfoundoffpring){
		#pick a random chromosome to start on
		$fiftyfifty=rand();
		if($fiftyfifty<0.5){
			$oneortwo=1;
		}else{
			$oneortwo=2;
		}
		if (@{$idsexhash{$id}}[0] == 1){
			#print STDOUT "Entered1\n";
			if ($oneortwo == 1){
			#	print STDOUT "Entered2\n";
			#	print STDOUT "@{$phase1{$kid}}[0],@{$phase1{$id}}[0]\n";
				@{$phase1{$kid}}[0] = @{$phase1{$id}}[0];
			#	print STDOUT "@{$phase1{$kid}}[0],@{$phase1{$id}}[0]\n";
			}else{
			#	print STDOUT "Entered3\n";
				@{$phase1{$kid}}[0] = @{$phase2{$id}}[0];
			}
			$ready1{$kid}=1;
		}else{
			#print STDOUT "Entered2\n";
			if ($oneortwo == 1){
				@{$phase2{$kid}}[0] = @{$phase1{$id}}[0];
			}else{
				@{$phase2{$kid}}[0] = @{$phase2{$id}}[0];
			}
			$ready2{$kid}=1;
		}	
		if($ready1{$kid}==1 & $ready2{$kid}==1){
			#print STDOUT "Kid ready!\n";
			push @ready, $kid;
		}		
	}
}

for($i=0;$i<scalar(@offspr);$i++){
	if(@sequenced[$i] ne '1'){
		next;
	}else{
		#print STDOUT "Entered1\n";
		$id=@offspr[$i];
		$onecount=0;
		$twocount=0;
		for($j=0;$j<scalar(@offspr);$j++){
			#print STDOUT "Entered2\n";
			if($i==$j){
		#		print STDOUT "Nexted\n";
				next;
			}elsif(@sequenced[$j] ne '1'){
				next;
			}else{
		#		print STDOUT "Entered3\n";
				$id2=@offspr[$j];
				if(@{$phase1{$id}}[0] eq @{$phase1{$id2}}[0]){
		#			print STDOUT "Entered4\n";
		#			print STDOUT "@{$phase1{$id}}[0],@{$phase1{$id2}}[0]\n";
					$onecount++;
				}elsif(@{$phase1{$id}}[0] eq @{$phase2{$id2}}[0]){
		#			print STDOUT "Entered5\n";
					$onecount++;
				}elsif(@{$phase2{$id}}[0] eq @{$phase1{$id2}}[0]){
		#			print STDOUT "Entered6\n";
					$twocount++;
				}elsif(@{$phase2{$id}}[0] eq @{$phase2{$id2}}[0]){
		#			print STDOUT "Entered7\n";
					$twocount++;
				}
			}
		}
		#@{$onecountsum{$id}}[0]=@{$onecountsum{$id}}[0]+$onecount;
		#@{$twocountsum{$id}}[0]=@{$twocountsum{$id}}[0]+$twocount;
		if($onecount eq '0' & $twocount eq '0'){
			
			@{$neithershared{$id}}[0]=@{$neithershared{$id}}[0]+1;
		}
	}
}
}
for($i=0;$i<scalar(@offspr);$i++){
	if(@sequenced[$i] ne '1'){
		next;
	}else{
	$id=@offspr[$i];
	$percentacceptable=1-(@{$neithershared{$id}}[0]/@ARGV[1]);
	print STDOUT "$id,$percentacceptable\n";
	
	}
}




######

sub find_kids{
	my $thisperson = $_[0];
	for(my $i = 0; $i < $numberofids;$i++){
		if (@paternal[$i] eq $thisperson || @maternal[$i] eq $thisperson){
			push @currfoundoffpring, @offspr[$i];			
		}
	}
}




















