package SEQUENCE;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(OpenFile ReadFa ReadFq WriteFq);

# return file head
sub OpenFile{
        my ($file,$property) = @_;
        my $filehand;
        if ($property eq "r"){
                open($filehand,"<$file") or die "OpenError: Read \'$file\',$!";
        }elsif ($property eq "w"){
                open($filehand,">$file") or die "OpenError: Write \'$file\',$!";
        }else{
                die "PropertyError: Unknow property \'$property\',$!";
        }
        return(\*$filehand);
}

# SEQUENCE::Readfa(\%hash_seq,$file) or SEQUENCE::Readfa(\@list_seq,$file)
sub ReadFa{
    my ($refer,$file) = @_;
    open(FA,"<$file") or die "OpenError: Read \'$file\',$!";
    $/ = ">";
    if(ref($refer) eq "ARRAY"){
	while(<FA>){
	    chomp;
	    next if ($_ eq "");
	    my $name = (split(/\s+/,$_))[0];
	    my $seq = (split(/\n/,$_,2))[1];
	    $name =~ s/\r//g;
	    $seq = uc($seq);
	    $seq =~ s/[\n\r]//g;
	    push(@$refer,[$name,$seq]);
	}
	$/ = "\n";
	close FA;
	return();
    }elsif(ref($refer) eq "HASH"){
	while(<FA>){
	    chomp;
	    next if ($_ eq "");
	    my $name = (split(/\s+/,$_))[0];
	    my $seq = (split(/\n/,$_,2))[1];
	    $name =~ s/\r//g;
	    $seq = uc($seq);
	    $seq =~ s/[\n\r]//g;
	    $refer->{$name} = $seq;
	}
	$/ = "\n";
	close FA;
	return();
    }else{
        die "ParameterError: Use SEQUENCE::Readfa(\\\%hash_seq,\$file) or SEQUENCE::Readfa(\\\@list_seq,\$file)";
    }
}


# SEQUENCE::Readfq(\%hash_seq,$file) or SEQUENCE::Readfq(\@list_seq,$file)
sub ReadFq{
    my ($refer,$file) = @_;
    open(FQ,"<$file") or die "OpenError: Read \'$file\',$!";
    if(ref($refer) eq "ARRAY"){
	my $flag = 0;
	my ($name,$seq);
	while(<FQ>){
	    chomp;
	    if($flag == 0 and /^\@/){
		$name = $_;
		$name =~ s/\@//g;
		$flag ++;
	    }elsif($flag == 1){
		$seq = $_;
		$flag ++;
	    }elsif($flag == 2 and /^\+/){
		$flag ++;
	    }elsif($flag == 3){
		push(@$refer,[$name,$seq,$_]);
		$flag = 0;
	    }else{
		next;
	    }
	}
	close FQ;
	return();
    }elsif(ref($refer) eq "HASH"){
	my $flag = 0;
	my ($name,$seq);
	while(<FQ>){
	    chomp;
	    if($flag == 0 and /^\@/){
		$name = $_;
		$name =~ s/\@//g;
		$flag ++;
	    }elsif($flag == 1){
		$seq = $_;
		$flag ++;
	    }elsif($flag == 2 and /^\+/){
		$flag ++;
	    }elsif($flag == 3){
		push(@{$refer->{$name}},[$seq,$_]);
		$flag = 0;
	    }else{
		next;
	    }
	}
	close FQ;
	return();
    }else{
        die "ParameterError: Use SEQUENCE::Readfq(\\\%hash_seq,\$file) or SEQUENCE::Readfq(\\\@list_seq,\$file)";
    }
}

sub WriteFq{
	my ($fh,$name,$seq,$qual) = @_;
	print $fh "\@$name\n$seq\n\+\n$qual\n";
}



1;
