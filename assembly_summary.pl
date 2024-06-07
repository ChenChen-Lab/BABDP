#!/usr/bin/perl -w
use strict;
my $input=shift@ARGV;
my $output=shift@ARGV;
open IN,"$input" or die$!;
open OUT,">$output" or die$!;
$/="#";
<IN>;
print OUT "Sample\tTotalNum\tMinLen\tScfNum\tScfLen\tScfN50\tScfN90\tScfMax\tmScfNum\tmScfLen\tCtgNum\tCtgLen\tCtgN50\tCtgN90\tCtgMax\tGC(%)\tGapNum\tGapLen\tMaxGap\n";
while (<IN>){
	my $line=$_;
	chomp $line;
	my @tmp=split/\n/,$line;	
	my $sample=(split/\//,((split/\t/,$tmp[0])[0]))[-2];
	my $Scaffold=$tmp[3];
	$Scaffold=~s/^\s+//;
	$Scaffold=~s/\s+/\t/g;
	my $Contigs=$tmp[6];
	$Contigs=~s/^\s+//;
 	$Contigs=~s/\s+/\t/g;	
	my $Insert_Gap=$tmp[9];
	$Insert_Gap=~s/^\s+//;
	$Insert_Gap=~s/\s+/\t/g;
	print OUT "$sample\t$Scaffold\t$Contigs\t$Insert_Gap\n";
}
close IN;
close OUT;

