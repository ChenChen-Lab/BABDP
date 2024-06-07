#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw/abs_path/;
use FindBin qw($Bin $Script);

if (@ARGV<2) {
	print "perl $0 <in dir> <out dir> <min length>\n";
	exit;
}

my ($id, $od, $len)=@ARGV;
$len=1000 unless ($len);

unless (-e $od) {
	mkdir $od;
}

$id=abs_path($id);
$od=abs_path($od);

foreach my $contig (glob "$id/*/contigs.fasta") {
	my @tmp=split /\//, $contig;
	my $tag=$tmp[@tmp-2];
	print "perl $Bin/sequence.pl -contig $od/$tag"."_contigs.fas -length $len $contig\n";
}
