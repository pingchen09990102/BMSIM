#!/usr/bin/perl
use strict;
use warnings;

die "#usage:perl $0 <input.fa> <output.fa> \n" unless @ARGV==2; 

my $FA=$ARGV[0];
open (FA, "<$FA")or die "can't open $FA !";
my @G = <FA>;
close FA;
my $head=shift @G;#去掉表头，将其存入
chomp $head;
my @Gs;
foreach my $g (@G) 
	{
		chomp $g;
		if ($g=~/>.*/) 
			{
				push (@Gs,"X");#替换fa的表头为X
			}
		else
			{
				push (@Gs,$g);
			}
	}
my $G=join ("",@Gs);
#将大于100Ns的gap改为X
#$G =~ s/N+N{99}/X/ig;
#my $count=0;
#$count=$G =~ s/N+N{99}/X/ig;#替换>=100N的gap
#print "$count\n";
my $outfile=$ARGV[1];
open (OUT, ">$outfile")or die "can't open $outfile !";
print OUT "$head\n$G\n";
close OUT;
