#!/usr/local/bin/perl

open(INFILE1,"@ARGV[0]")||die;
open(OUTFILE1,">@ARGV[1]")||die;
$headerline=<INFILE1>;
while($line=<INFILE1>){
	chomp $line;
	@toks1=split(",",$line);
	print OUTFILE1 " 1 @toks1[0] 0 @toks1[2]\n";
}
close INFILE1;
close OUTFILE1;

open(INFILE2,"@ARGV[0]")||die;
open(INFILE3,"@ARGV[2]")||die;
open(OUTFILE2,">@ARGV[3]")||die;
#Read in pedigree file and store the parents for each kid
$pedheader=<INFILE3>;
while($pedline=<INFILE3>){
	chomp $pedline;
	@pedtoks=split(",",$pedline);
	$id=@pedtoks[0];
	@{$father{$id}}[0]=@pedtoks[1];
	@{$mother{$id}}[0]=@pedtoks[2];
	@{$sex{$id}}[0]=@pedtoks[3];
	@{$pednum{$id}}[0]=@pedtoks[5];
}

$j=0;
$headerline=<INFILE2>;
chomp $headerline;
@headertoks=split(",",$headerline);
while($line2=<INFILE2>){
	chomp $line2;
	@toks=split(",",$line2);
	for($i=3;$i<scalar(@toks);$i++){
		$id=@headertoks[$i];
		$R[$id][$j]=@toks[$i];
	}
	$j++;
}
$k=$j;
for($i=3;$i<scalar(@toks);$i++){
	$id=@headertoks[$i];
	print OUTFILE2 " @{$pednum{$id}}[0] $id @{$father{$id}}[0] @{$mother{$id}}[0] @{$sex{$id}}[0] 0";
	for($j=0;$j<$k;$j++){
		$genotype=$R[$id][$j];
		if($genotype eq '0'){
			print OUTFILE2 " 1 1";
		}elsif($genotype eq '1'){
			print OUTFILE2 " 1 2";
		}elsif($genotype eq '2'){
			print OUTFILE2 " 2 2";
		}else{
			print OUTFILE2 " 0 0";
		}
	}
	print OUTFILE2 "\n";
}
close INFILE2;
close INFILE3;
close OUTFILE2;
