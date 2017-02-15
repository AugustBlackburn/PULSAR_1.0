#!/usr/local/bin/perl

open (OUTFILE, ">out_012.csv") || die;
$rootname = $ARGV[0];
$POS = "$rootname".".012.pos";
$DATA = "$rootname".".012";
$INDV = "$rootname".".012.indv";

########TRANSPOSE THE 012 FILE###########
open INFILE, $DATA or die "File $DATA not found";
$outfile = "t"."$DATA";
open tINFILE, ">$outfile";
$i=0;
foreach $line (<INFILE>){
  chomp $line;
  @line=split /\t/,$line;
  for ($j=0;$j<=$#line;$j++){
    $R[$i][$j]=$line[$j];
  }
  $i++;
}
$nlin=$i-1;
$ncol=$#line;
for ($i=0;$i<=$ncol;$i++){
  for ($j=0;$j<=$nlin;$j++){
    print tINFILE "$R[$j][$i]";
    if ($j<$nlin){
      print tINFILE ",";
    }
  }
  print tINFILE "\n";
}
print "Transpose of $DATA created\n";
close INFILE;
close tINFILE;

###############Combine the files##################
open (INLOCS, $POS) || die;
open (INGENO, $outfile) || die;
open (INID, $INDV) || die;
print OUTFILE "Marker,Chr,Position,";
while($line = <INID>){
	chomp $line;
	push @headerlist, $line;
}
$headlisttoprint=join(",",@headerlist);
print OUTFILE "$headlisttoprint";
print OUTFILE "\n";
close INID;
#read in header line
$line = <INGENO>;
while($lineloc = <INLOCS>){
	chop $lineloc;
	@toklocs = split('\t', $lineloc);
	$marker = "$toklocs[0]" . "_" . "$toklocs[1]";
	print OUTFILE "$marker,";
	print OUTFILE "$toklocs[0],$toklocs[1],";
	$linegeno = <INGENO>;
	$linegeno =~ s/-1/-/g;
	print OUTFILE $linegeno;
}

close INLOCS;
close INGENO;
close OUTFILE;













