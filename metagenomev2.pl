#!/usr/bin/perl 
=head1 Description
	metagenome
=head1 Usage
	perl metagenome.pl [options]
	general input and output:
	-fql  <a file>       list of fastq from rawdata(file suffix can be ".fq.gz"), one strain per line, PE read seperated by ",", and different by "\n"
	-host               host[human, mice or none]
	-o   <a directory>  output dir, default current directory [./]
	-annlist	    annotation database	[default vfdb,resfinder]
	-thd <num>          thread for dsub
	-h                  show help
	-notrun		    only write the shell, but not run
=head1 Example
	perl metagenome.pl  -fql fq.list -host human  -o ./ -annlist  annotation.list

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

my ($fqlist,$host,$outdir,$annotation,$thd,$help,$notrun);

GetOptions (
	"fql:s"=>\$fqlist,
	"host:s"=>\$host,
	"o:s"=>\$outdir,
	"annlist:s"=>\$annotation,
	"h"=>\$help,
	"notrun"=>\$notrun
);
die `pod2text $0` if ($help || !$fqlist);

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
my $kraken_02="$outdir/02.kraken";
my $assembly_03="$outdir/03.assembly";
my $annotation_04="$outdir/04.annotation";
for($shelldir,$fastqc_00,$qc_01,$kraken_02,$assembly_03,$annotation_04){
    (-s $_) || mkdir $_;
    $_ = abs_path($_);
}
#============================================================================================================================================================
if (!$fqlist){
	print "No input read files\n";
	exit;
}

if (!$host){
        print "No input host type: 'human' or 'mice' or 'none'\n";
        exit;
}

my $fastqc="/datapool/software/anaconda3/envs/qiime2/bin/fastqc";
my $bowtie2="/datapool/software/anaconda3/bin/bowtie2";
my $kraken2="/datapool/software/anaconda3/bin/kraken2";
my $metaspades="/datapool/software/anaconda3/bin/metaspades.py";
my $fastqc_summary="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/fastqc_summary.pl";
my $samtools="/datapool/software/anaconda3/bin/samtools";
#my $pick_contigs_from_spades_result="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/pick_contigs_from_spades_result.pl";
my $assembly_summary="/datapool/bioinfo/guchaoyang/pipeline/metagenome/bin/assembly_summary.pl";
my $seqkit="/datapool/software/anaconda3/envs/qiime2/bin/seqkit";
my $run_abricate="/datapool/stu/yuejl/genome-analysis-pipeline/bin/run_abricate.pl";

open (OUT,">$shelldir/metagenome.sh")||die;

my $sample;
if ($fqlist && $host eq "human" && $annotation eq "vfdb_resfinder"){


open (OUT0,">$shelldir/step0_fastqc.sh")||die;
open (OUT1,">$shelldir/step1_nohost.sh")||die;
open (OUT1_1,">$shelldir/step1_1_nohost.sh")||die;
open (OUT2,">$shelldir/step2_kraken.sh")||die;
open (OUT3,">$shelldir/step3_assembly.sh")||die;
open (OUT3_1,">$shelldir/step3_1_assemblystats.sh")||die;
open (OUT4,">$shelldir/step4_annotation.sh")||die;
open (OUT_T, ">$assembly_03/pick_1k_sequence.sh")||die;

	(-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
	(-d "$qc_01/nohostqc_result") || `mkdir "$qc_01/nohostqc_result"`;
	(-d "$kraken_02/kraken_output") || `mkdir "$kraken_02/kraken_output"`;
	(-d "$assembly_03/assembly") || `mkdir "$assembly_03/assembly"`;
	(-d "$assembly_03/final_contig") || `mkdir "$assembly_03/final_contig"`;
	(-d "$annotation_04/vfdb") ||`mkdir "$annotation_04/vfdb"`;
	(-d "$annotation_04/resfinder") ||`mkdir "$annotation_04/resfinder"`;
	open (LIST,">$annotation_04/genomes.list") ||die;
	open (IN,$fqlist)||die;
	while(my $line=<IN>){
		chomp($line);
		$sample=(split/\t/,$line)[0];
		my $fq=(split/\t/,$line)[1];
		my $fq1=(split/\,/,$fq)[0];
		my $fq2=(split/\,/,$fq)[1];
                
		my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result";
                print OUT0 $cmd_0."\n";
		
		my $dir1="$qc_01/$sample";
		my $dir2="$assembly_03/assembly/$sample/";
		(-d $dir1) || `mkdir $dir1`;
		(-d $dir2) || `mkdir $dir2`;

	my $cmd_1="cd $qc_01/$sample\n$bowtie2  --very-sensitive  -p 10 -x /datapool/db/hg38/hg38  -1 $fq1 -2 $fq2 --al-conc-gz $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -F 4 -Sb $sample.sam  > $sample.bam\n$samtools  sort $sample.bam  -o $sample.sort.bam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz";
	print OUT1 $cmd_1."\n";

	my $cmd_2="$kraken2  -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $qc_01/$sample/$sample.unmap.1.fq.gz  $qc_01/$sample/$sample.unmap.2.fq.gz  --threads 10 --output $kraken_02/kraken_output/$sample.kraken --use-names --report $kraken_02/kraken_output/$sample.label --use-mpa-style";
	print OUT2 $cmd_2."\n";

	my $cmd_1_1="$fastqc --extract $qc_01/$sample/$sample.unmap.1.fq.gz -o $qc_01/nohostqc_result\n$fastqc --extract $qc_01/$sample/$sample.unmap.2.fq.gz -o $qc_01/nohostqc_result";
	print OUT1_1 $cmd_1_1."\n";

	my $cmd_3="$metaspades -1 $qc_01/$sample/$sample.unmap.1.fq.gz -2 $qc_01/$sample/$sample.unmap.2.fq.gz  --only-assembler  -o $assembly_03/assembly/$sample";
	print OUT3 $cmd_3."\n";

	my $cmd_T="$seqkit  seq -m 1000 $assembly_03/assembly/$sample/metaspades/contigs.fasta > $assembly_03/final_contig/$sample\_contigs.fas";
	print OUT_T $cmd_T."\n";

	my $cmd_3_1="$seqkit stats $assembly_03/assembly/$sample/contigs.fasta >> $assembly_03/metaspades.txt";
        print OUT3_1 $cmd_3_1."\n";

	print LIST "$assembly_03/final_contig/$sample\_contigs.fas\n";
	}
	my $cmd_4="perl $fastqc_summary $fastqc_00/qc_result >$fastqc_00/fastqc_summary.txt\nperl $fastqc_summary $qc_01/nohostqc_result >$qc_01/fastqc_summary.txt\nsh $assembly_03/pick_1k_sequence.sh &\necho \"all steps finished\" >finished.txt\n";
	$cmd_4.="perl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db vfdb  -o $annotation_04/vfdb\nperl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db resfinder_new  -o $annotation_04/resfinder";
	print OUT4 $cmd_4."\n";

        close IN;
        close OUT0;
        close OUT1;
        close OUT1_1;
        close OUT2;
        close OUT3;
	close OUT_T;        
	close OUT3_1;
        close OUT4;

	 my $cmd="### step0_fastqc\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_fastqc\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step0_fastqc.sh \ndate +\"\%D \%T -> Finish 0) step0_fastqc\" >>$shelldir/log\n";
        $cmd.="### step1_nohost\ncd $shelldir\ndate +\"\%D \%T -> Start 1) step1_nohost\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step1_nohost.sh \ndate +\"\%D \%T -> Finish 1) step1_nohost\" >>$shelldir/log\n";
        $cmd.="### step2_kraken\ncd $shelldir\ndate +\"\%D \%T -> Start 2) step2_kraken\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step2_kraken.sh \ndate +\"\%D \%T -> Finish 2) step2_kraken\" >>$shelldir/log\n";
        $cmd.="### step1_1_nohostqc\ncd $shelldir\ndate +\"\%D \%T -> Start 3) step1_1_nohostqc\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step1_1_nohostqc.sh \ndate +\"\%D \%T -> Finish 3) step1_1_nohostqc\" >>$shelldir/log\n";
        $cmd.="### step3_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_assembly\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_assembly.sh \ndate +\"\%D \%T -> Finish 4) step3_assembly\" >>$shelldir/log\n";
	$cmd.="### step3_1_assemblystats\ncd $shelldir\ndate +\"\%D \%T -> Start 5) step3_1_assemblystats\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_1_assemblystats.sh \ndate +\"\%D \%T -> Finish 5) step3_1_assemblystats\" >>$shelldir/log\n";
 
       print OUT $cmd. "### step4_annotation\ncd $shelldir\n",
          "date +\"\%D \%T -> Start 6) step4_annotation\" >>$shelldir/log\n",
          "nohup sh $shelldir/step4_annotation.sh >$shelldir/step4_annotation.log \n",
          "date +\"\%D \%T -> Finish 6) step4_annotation\" >>$shelldir/log\n";


}


if ($fqlist && $host eq "mice" && $annotation eq "vfdb_resfinder"){
open (OUT0,">$shelldir/step0_fastqc.sh")||die;
open (OUT1,">$shelldir/step1_nohost.sh")||die;
open (OUT1_1,">$shelldir/step1_1_nohost.sh")||die;
open (OUT2,">$shelldir/step2_kraken.sh")||die;
open (OUT3,">$shelldir/step3_assembly.sh")||die;
open (OUT3_1,">$shelldir/step3_1_assemblystats.sh")||die;
open (OUT4,">$shelldir/step4_annotation.sh")||die;
open (OUT_T, ">$assembly_03/pick_1k_sequence.sh")||die;

        (-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
 	(-d "$qc_01/nohostqc_result") || `mkdir "$qc_01/nohostqc_result"`;
	(-d "$kraken_02/kraken_output") || `mkdir "$kraken_02/kraken_output"`;
	(-d "$assembly_03/assembly") || `mkdir "$assembly_03/assembly"`;
        (-d "$assembly_03/final_contig") || `mkdir "$assembly_03/final_contig"`;
        (-d "$annotation_04/vfdb") ||`mkdir "$annotation_04/vfdb"`;
        (-d "$annotation_04/resfinder") ||`mkdir "$annotation_04/resfinder"`;
        open (LIST,">$annotation_04/genomes.list") ||die;
        open (IN,$fqlist)||die;
        while(my $line=<IN>){
                chomp($line);
                $sample=(split/\t/,$line)[0];
                my $fq=(split/\t/,$line)[1];
                my $fq1=(split/\,/,$fq)[0];
                my $fq2=(split/\,/,$fq)[1];

                my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result";
                print OUT0 $cmd_0."\n";

                my $dir1="$qc_01/$sample";
                my $dir2="$assembly_03/assembly/$sample/";
                (-d $dir1) || `mkdir $dir1`;
                (-d $dir2) || `mkdir $dir2`;
        my $cmd_1="cd $qc_01/$sample\n$bowtie2  --very-sensitive  -p 10 -x /datapool/db/mm39  -1 $fq1 -2 $fq2 --al-conc-gz $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -F 4 -Sb $sample.sam  > $sample.bam\n$samtools  sort $sample.bam  -o $sample.sort.bam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz";
        print OUT1 $cmd_1."\n";
	my $cmd_2="$kraken2  -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $qc_01/$sample/$sample.unmap.1.fq.gz  $qc_01/$sample/$sample.unmap.2.fq.gz  --threads 10 --output $kraken_02/kraken_output/$sample.kraken --use-names --report $kraken_02/kraken_output/$sample.label --use-mpa-style";
        print OUT2 $cmd_2."\n";
        my $cmd_1_1="$fastqc --extract $qc_01/$sample/$sample.unmap.1.fq.gz -o $qc_01/nohostqc_result\n$fastqc --extract $qc_01/$sample/$sample.unmap.2.fq.gz -o $qc_01/nohostqc_result";
        print OUT1_1 $cmd_1_1."\n";	

	my $cmd_3="$metaspades -1 $qc_01/$sample/$sample.unmap.1.fq.gz -2 $qc_01/$sample/$sample.unmap.2.fq.gz  --only-assembler  -o $assembly_03/assembly/$sample";
        print OUT3 $cmd_3."\n";

        my $cmd_T="$seqkit  seq -m 1000 $assembly_03/assembly/$sample/metaspades/contigs.fasta > $assembly_03/final_contig/$sample\_contigs.fas";
        print OUT_T $cmd_T."\n";

	my $cmd_3_1="$seqkit stats $assembly_03/assembly/$sample/contigs.fasta >> $assembly_03/metaspades.txt";
        print OUT3_1 $cmd_3_1."\n"; 

       print LIST "$assembly_03/final_contig/$sample\_contigs.fas\n";
        }
        my $cmd_4="perl $fastqc_summary $fastqc_00/qc_result >$fastqc_00/fastqc_summary.txt\nperl $fastqc_summary $qc_01/nohostqc_result >$qc_01/fastqc_summary.txt\nsh $assembly_03/pick_1k_sequence.sh &\necho \"all steps finished\" >finished.txt\n";
        $cmd_4.="perl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db vfdb  -o $annotation_04/vfdb\nperl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db resfinder_new  -o $annotation_04/resfinder";
        print OUT4 $cmd_4."\n";

        close IN;
        close OUT0;
        close OUT1;
	close OUT1_1;
	close OUT2;
        close OUT3;
        close OUT_T;
	close OUT3_1;
        close OUT4;

	 my $cmd="### step0_fastqc\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_fastqc\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step0_fastqc.sh \ndate +\"\%D \%T -> Finish 0) step0_fastqc\" >>$shelldir/log\n";
        $cmd.="### step1_nohost\ncd $shelldir\ndate +\"\%D \%T -> Start 1) step1_nohost\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step1_nohost.sh \ndate +\"\%D \%T -> Finish 1) step1_nohost\" >>$shelldir/log\n";
        $cmd.="### step2_kraken\ncd $shelldir\ndate +\"\%D \%T -> Start 2) step2_kraken\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step2_kraken.sh \ndate +\"\%D \%T -> Finish 2) step2_kraken\" >>$shelldir/log\n";
        $cmd.="### step1_1_nohostqc\ncd $shelldir\ndate +\"\%D \%T -> Start 3) step1_1_nohostqc\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step1_1_nohostqc.sh \ndate +\"\%D \%T -> Finish 3) step1_1_nohostqc\" >>$shelldir/log\n";
        $cmd.="### step3_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_assembly\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_assembly.sh \ndate +\"\%D \%T -> Finish 4) step3_assembly\" >>$shelldir/log\n";
	$cmd.="### step3_1_assemblystats\ncd $shelldir\ndate +\"\%D \%T -> Start 5) step3_1_assemblystats\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_1_assemblystats.sh \ndate +\"\%D \%T -> Finish 5) step3_1_assemblystats\" >>$shelldir/log\n";

        print OUT $cmd. "### step4_annotation\ncd $shelldir\n",
          "date +\"\%D \%T -> Start 6) step4_annotation\" >>$shelldir/log\n",
          "nohup sh $shelldir/step4_annotation.sh >$shelldir/step4_annotation.log \n",
          "date +\"\%D \%T -> Finish 6) step4_annotation\" >>$shelldir/log\n";

}


if ($fqlist && $host eq "none" && $annotation eq "vfdb_resfinder"){
open (OUT0,">$shelldir/step0_fastqc.sh")||die;
open (OUT2,">$shelldir/step2_kraken.sh")||die;
open (OUT3,">$shelldir/step3_assembly.sh")||die;
open (OUT3_1,">$shelldir/step3_1_assemblystats.sh")||die;
open (OUT4,">$shelldir/step4_annotation.sh")||die;
open (OUT_T, ">$assembly_03/pick_1k_sequence.sh")||die;

        (-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
        (-d "$kraken_02/kraken_output") || `mkdir "$kraken_02/kraken_output"`;
        (-d "$assembly_03/assembly") || `mkdir "$assembly_03/assembly"`;
        (-d "$assembly_03/final_contig") || `mkdir "$assembly_03/final_contig"`;
        (-d "$annotation_04/vfdb") ||`mkdir "$annotation_04/vfdb"`;
        (-d "$annotation_04/resfinder") ||`mkdir "$annotation_04/resfinder"`;
        open (LIST,">$annotation_04/genomes.list") ||die;
        open (IN,$fqlist)||die;
        while(my $line=<IN>){
                chomp($line);
                $sample=(split/\t/,$line)[0];
                my $fq=(split/\t/,$line)[1];
                my $fq1=(split/\,/,$fq)[0];
                my $fq2=(split/\,/,$fq)[1];

                my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result";
                print OUT0 $cmd_0."\n";

                my $dir2="$assembly_03/assembly/$sample/";
                (-d $dir2) || `mkdir $dir2`;
        my $cmd_2="$kraken2  -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $fq1  $fq2  --threads 10 --output $kraken_02/kraken_output/$sample.kraken --use-names --report $kraken_02/kraken_output/$sample.label --use-mpa-style";
        print OUT2 $cmd_2."\n";

        my $cmd_3="$metaspades -1 $fq1 -2 $fq2  --only-assembler  -o $assembly_03/assembly/$sample";
        print OUT3 $cmd_3."\n";

        my $cmd_T="$seqkit  seq -m 1000 $assembly_03/assembly/$sample/metaspades/contigs.fasta > $assembly_03/final_contig/$sample\_contigs.fas";
        print OUT_T $cmd_T."\n";

        
	my $cmd_3_1="$seqkit stats $assembly_03/assembly/$sample/contigs.fasta >> $assembly_03/metaspades.txt";
        print OUT3_1 $cmd_3_1."\n";

	print LIST "$assembly_03/final_contig/$sample\_contigs.fas\n";
        }
        my $cmd_4="perl $fastqc_summary $fastqc_00/qc_result >$fastqc_00/fastqc_summary.txt\nsh $assembly_03/pick_1k_sequence.sh &\necho \"all steps finished\" >finished.txt\n";
        $cmd_4.="perl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db vfdb  -o $annotation_04/vfdb\nperl $run_abricate -abricate  -sl $annotation_04/genomes.list  -db resfinder_new  -o $annotation_04/resfinder";
        print OUT4 $cmd_4."\n";


	close IN;
	close OUT0;
	close OUT2;
	close OUT3;
	close OUT_T;
	close OUT3_1;
	close OUT4;

	my $cmd="### step0_fastqc\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_fastqc\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step0_fastqc.sh \ndate +\"\%D \%T -> Finish 0) step0_fastqc\" >>$shelldir/log\n";
	$cmd.="### step2_kraken\ncd $shelldir\ndate +\"\%D \%T -> Start 2) step2_kraken\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step2_kraken.sh \ndate +\"\%D \%T -> Finish 2) step2_kraken\" >>$shelldir/log\n";
	$cmd.="### step3_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 3) step3_assembly\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_assembly.sh \ndate +\"\%D \%T -> Finish 3) step3_assembly\" >>$shelldir/log\n";
	$cmd.="### step3_1_assemblystats\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_1_assemblystats\" >>$shelldir/log\nperl $Bin/dsub_batch2.pl -thd 10 -mem 4 $shelldir/step3_1_assemblystats.sh \ndate +\"\%D \%T -> Finish 4) step3_1_assemblystats\" >>$shelldir/log\n";
	print OUT $cmd. "### step4_annotation\ncd $shelldir\n",
          "date +\"\%D \%T -> Start 5) step4_annotation\" >>$shelldir/log\n",
          "nohup sh $shelldir/step4_annotation.sh >$shelldir/step4_annotation.log \n",
          "date +\"\%D \%T -> Finish 5) step4_annotation\" >>$shelldir/log\n";

}
close OUT;
	
$notrun && exit;
