#!/usr/bin/perl -w
use File::Basename;

if (@ARGV!=1){
	print "perl $0 <dir>\n";
	exit;
}

my ($dir)=@ARGV;
my (%hash,%file);
print "Filename\tTotal Sequences\tSequence length\tBasic Statistics\tPer base sequence quality\tPer tile sequence quality\tPer sequence quality scores\tPer base sequence content\tPer sequence GC content\tPer base N content\tSequence Length Distribution\tSequence Duplication Levels\tOverrepresented sequences\tAdapter Content\n";
@name=("Filename","Total Sequences","Sequence length","Basic Statistics","Per base sequence quality","Per tile sequence quality","Per sequence quality scores","Per base sequence content","Per sequence GC content","Per base N content","Sequence Length Distribution","Sequence Duplication Levels","Adapter Content");

my @file1=glob "$dir/*/summary.txt";
for my $file (@file1){
	open (IN,$file)||die;
	while (<IN>){
		chomp;
		my ($result,$label,$file)=(split /\t/,$_)[0,1,2];
		$hash{$file}{$label}=$result;
		$file{$file}=1;
	}
	close IN;
}

my @file2=glob "$dir/*/fastqc_data.txt";
for my $file (@file2){
	my $tag=(split /\//,dirname($file))[-1];
	$tag=~s/_fastqc/.fq.gz/g;
	open (IN1,$file)||die;
	while (<IN1>){
		chomp;
		if ($_=~/Filename/){
			my $filename=(split /\t/,$_)[1];
			$tmp_1="Filename";
			$hash{$tag}{$tmp_1}=$filename;
		}
		if ($_=~/Total Sequences/){
			my $total_seq=(split /\t/,$_)[1];
			$tmp_2="Total Sequences";
			$hash{$tag}{$tmp_2}=$total_seq;
		}
		if ($_=~/Sequence length/){
			my $seq_len=(split /\t/,$_)[1];
			$tmp_3="Sequence length";
			$hash{$tag}{$tmp_3}=$seq_len;
		}
	}
	close IN1;
}

for my $a (keys %file){
	#print "$a\t";
	for my $b (@name){
		if (exists $hash{$a}{$b}){
			print "$hash{$a}{$b}\t";
		}else{
			print "NA\t";
		}
	}
	print "\n";
}

