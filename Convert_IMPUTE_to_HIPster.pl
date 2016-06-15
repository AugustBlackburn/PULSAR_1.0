#!/usr/local/bin/perl

open(INFILE1,"@ARGV[0]")||die;
open(INFILE2,"@ARGV[1]")||die;
$headerline=<INFILE1>;
$headerline2=<INFILE1>;

while($line=<INFILE1>){
	chomp $line;
	@toks1=split(/\s+/,$line);
	$new1="@toks1[1]"."_1";
	$new2="@toks1[1]"."_2";
	#print STDOUT "$new1 $new2\n";
	push @haplist, $new1;
	push @haplist, $new2;
}	
close INFILE1;

open(OUTFILE,">Phased_chromosomes.csv")||die;
print OUTFILE "Marker,Chr,Position";
foreach $hap (@haplist){
	print OUTFILE ",$hap";
}
print OUTFILE "\n";

while($line = <INFILE2>){
	chomp $line;
	@toks=split(/\s+/,$line);
	print OUTFILE "@toks[1],@toks[0],@toks[2]";
	if(@toks[3] eq '1'){
		for($i=5;$i<scalar(@toks);$i++){
			print OUTFILE ",@toks[$i]";
		}
	}else{
		for($i=5;$i<scalar(@toks);$i++){
			if(@toks[$i] eq '1'){
				print OUTFILE ",0";
			}elsif(@toks[$i] eq '0'){
				print OUTFILE ",1";
			}else{
				print OUTFILE ",";
			}
		}
	}
	print OUTFILE "\n";
}
close INFILE2;
close OUTFILE;
