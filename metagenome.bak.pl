#!/usr/bin/perl 
=head1 Description
	metagenome
=head1 Usage
	perl metagenome.pl [options]
	general input and output:
	-sl  <a file>       list of fastq from rawdata(file suffix can be ".fq.gz"), one strain per line, PE read seperated by ",", and different by "\n"
	-host               host[human, mice or none]
	-o   <a directory>  output dir, default current directory [./]
	-annlist	    annotation database	[default vfdb,resfinder]
	-thd <num>          thread for dsub
	-h                  show help
	-notrun		    only write the shell, but not run
=head1 Example
	perl metagenome.pl  -sl sl.list -host human  -o ./ -annlist  annotation.list

=head1 Version
        Author: guchaoyang0826@163.com
        Date: 2022-11-21
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw/abs_path/;

my ($sllist,$host,$outdir,$annotationlist,$thd,$help,$notrun);

GetOptions (
	"sl:s"=>\$sllist,
	"host:s"=>\$host,
	"o:s"=>\$outdir,
	"annlist:s"=>\$annotationlist,
	"h"=>\$help,
	"notrun"=>\$notrun
);
die `pod2text $0` if ($help || !$sllist);

unless ($outdir) {
	$outdir="./";
}

unless (-e $outdir) {
	mkdir $outdir;
}

unless ($thd){
	$thd=8;
}

$outdir=abs_path($outdir);

my $shelldir="$outdir/shell";
my $fastqc_00="$outdir/00.fastqc";
my $qc_01="$outdir/01.qc";
my $assembly_02="$outdir/02.assembly";
my $annotation_03="$outdir/03.annotation";
for($shelldir,$fastqc_00,$qc_01,$assembly_02,$annotation_03){
    (-s $_) || mkdir $_;
    $_ = abs_path($_);
}
#============================================================================================================================================================
if (!$sllist){
	print "No input read files\n";
	exit;
}

if (!$host){
        print "No input host type: 'human' or 'mice' or 'none'\n";
        exit;
}

my $fastqc="/datapool/software/anaconda3/envs/qiime2/bin/fastqc";
my $bowtie2="/datapool/software/anaconda3/bin/bowtie2";
my $metaspades="/datapool/software/anaconda3/bin/metaspades.py";
my $fastqc_summary="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/fastqc_summary.pl";
my $pick_contigs_from_spades_result="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/pick_contigs_from_spades_result.pl";
my $assembly_summary="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/assembly_summary.pl";
my $seqkit="/datapool/software/anaconda3/envs/qiime2/bin/seqkit";
my $run_abricate="/datapool/stu/yuejl/genome-analysis-pipeline/bin/run_abricate.pl";

open (OUT0,">$shelldir/step0_fastqc.sh")||die;
open (OUT1,">$shelldir/step1_qc.sh")||die;
open (OUT2,">$shelldir/step2_assembly.sh")||die;
open (OUT3,">$shelldir/step3_annotation.sh")||die;
open (OUT,">$shelldir/metagenome.sh")||die;

my ($sample,@fq,$annotation_db,@annotation_db);
if ($sllist && $host eq "human"){
	(-d "$fastqc_00/result") || `mkdir "$fastqc_00/result"`;
	(-d "$assembly_02/assembly") || `mkdir "$assembly_02/assembly"`;
	(-d "$assembly_02/final_contig") || `mkdir "$assembly_02/final_contig"`;
	open (IN,$sllist)||die;
	while(my $line=<IN>){
		chomp($line);
		$sample=(split/\t/,$line)[0];
		my $fq=(split/\t/,$line)[1];
		my $fq1=(split/\,/,$fq)[0];
		my $fq2=(split/\,/,$fq)[1];
		push @fq,$fq1;
		push @fq,$fq2;
		for my $i (@fq){
                	my $cmd_0="$fastqc --extract $i -o $fastqc_00/result";
                        print OUT0 $cmd_0."\n";
		}
	my $dir1="$qc_01/$sample";
	my $dir2="$assembly_02/assembly/$sample/";
	(-d $dir1) || `mkdir $dir1`;
	(-d $dir2) || `mkdir $dir2`;

	my $cmd_1="cd $qc_01/$sample\n$bowtie2  --very-sensitive  -p 10 -x /datapool/db/hg38/hg38  -1 $fq1 -2 $fq2 --al-conc-gz $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz";
	print OUT1 $cmd_1."\n";
	my $cmd_2="$metaspades -1 $qc_01/$sample/$sample.unmap.1.fq.gz -2 $qc_01/$sample/$sample.unmap.2.fq.gz  -o $assembly_02/assembly/$sample/metaspades";
	print OUT2 $cmd_2."\n";
	}

	my $cmd_3="perl $fastqc_summary $fastqc_00/result >$fastqc_00/fastqc_summary.txt\nperl $pick_contigs_from_spades_result $assembly_02/assembly  $assembly_02/final_contig/ > $assembly_02/pick_1k_contig.sh\nsh $assembly_02/pick_1k_contig.sh > $assembly_02/ass_stats.txt\nperl $assembly_summary\n$seqkit stats $assembly_02/final_contig/$sample.fa >> $assembly_02/seqkit.txt\nls $assembly_02/final_contig/*.fa >> genomes.list\necho \"all steps finished\" >finished.txt\n";
	open (G,$annotationlist) || die;
        while(my $line2=<G>){
                chomp$line2;
                $annotation_db=$line2;
                push @annotation_db,$annotation_db;
	}
	for my $l(@annotation_db){
		(-s $annotation_03/$l) || `mkdir "$annotation_03/$l"`;
		$cmd_3.="perl $run_abricate -abricate  -sl genomes.list  -db $l  -o $annotation_03/$l";
		print OUT3 $cmd_3."\n";
	}
}
	close IN;
	close OUT0;
	close OUT1;
	close OUT2;
	close OUT3;

	my $cmd="### step0_fastqc\ncd $shelldir\ndate + \"\%D \%T -> Start 0) step0_fastqc\" >>$shelldir/log\n$Bin/dsub_batch.pl -thd 10 -mem 4 $shelldir/step0_fastqc.sh\ndate + \"\%D \%T -> Finish 0) step0_fastqc\" >>$shelldir/log\n";
	$cmd.="### step1_qc\ncd $shelldir\ndate + \"\%D \%T -> Start 1) step1_qc\" >>$shelldir/log\n$Bin/dsub_batch.pl -thd 10 -mem 4 $shelldir/step1_qc.sh\ndate + \"\%D \%T -> Finish 1) step1_qc\" >>$shelldir/log\n";
	$cmd.="### step2_assembly\ncd $shelldir\ndate + \"\%D \%T -> Start 2) step2_assembly\" >>$shelldir/log\n$Bin/dsub_batch.pl -thd 10 -mem 4 $shelldir/step2_assembly.sh\ndate + \"\%D \%T -> Finish 2) step2_assembly\" >>$shelldir/log\n";
	print OUT $cmd. "### step3_annotation\ncd $shelldir\n",
          "date + \"\%D \%T -> Start 3) step3_annotation\" >>$shelldir/log\n",
          "$Bin/dsub_batch.pl  -thd 10 -mem 4 $shelldir/step3_annotation.sh \n",
          "date + \"\%D \%T -> Finish 3) step3_annotation\" >>$shelldir/log\n";

close OUT;
	
$notrun && exit;
