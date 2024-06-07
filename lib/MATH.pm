package MATH;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Sum Avg NumN50 NumN90 Max Min FormatInt CharNum);

sub Sum{
    my $sum = 0;
    foreach my $var (@_) {
        $sum += $var;
    }
    return $sum;
}

sub Avg{
    return(sum(@_)/@_);
}

sub NumN50{
    my @len = @_;
    my ($total,$tmp) = (Sum(@len),0);
    foreach (sort {$b<=>$a} @len){
		#print "$_\n";
        $tmp += $_;
        if($tmp >= $total/2){
            return($_);
        }
    }
}

sub NumN90{
    my @len = @_;
    my ($total,$tmp) = (Sum(@len),0);
    foreach (sort {$b<=>$a} @len){
        $tmp += $_;
        if($tmp >= $total*0.9){
            return($_);
        }
    }
}

sub SeqN50{
    return;
}

sub Max{
    my $max = shift;
    foreach my $var (@_) {
        $max = $var if ($var > $max);
    }
    return $max;
}

sub Min{
    my $min = shift;
    foreach my $var (@_) {
        $min = $var if ($var < $min);
    }
    return $min;
}

sub FormatInt{
    my $num = shift;
    if($num >= 1000){
        1 while $num =~ s/^(-?\d+)(\d{3})/$1,$2/;
    }
    return $num;
}

sub CharNum{
    my ($string,$char) = @_;
    $string =~ s/~$char//g;
    return length($string);
}


1;
