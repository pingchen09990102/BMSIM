#!/usr/bin/env perl
package ErrorFP;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);

sub FP
{
	my ($sublamdaFP,$mol_FP) = @_; 
	my $FpApart=0;
	my $FpPos1=0;
	my $FpPos2=0;
	my @sitePos=();
	my @insert=();
	$FpApart=random_exponential(1,$sublamdaFP);#the initial interval of FP, the unit kb should transfer to bp
	$FpPos2=int($FpPos1+$FpApart*1000);#the loci of FP is the initial position plus the interval
	my $mol_len_FP=$$mol_FP[-1];#molecule length
	my $Numsites_FP=scalar(@$mol_FP)-1;#number of labels
	while ($FpPos2 < $mol_len_FP) #if the FP loci shorter than the molecule length
	{
		while (grep{$_ == $FpPos2} @$mol_FP) #Get rid of overlapping with the original point
		{
			$FpApart=random_exponential(1,$sublamdaFP);#produce a new FP site
			$FpPos2=int($FpPos1+$FpApart*1000);
		}
		$Numsites_FP++;
		push (@$mol_FP,$FpPos2);#insert the FP
		#search for the loci of insert site
		#@insert=grep{$sitePos[$_]>$FpPos2} 0..$#sitePos;
		#splice(@$mol_FP,$insert[0],0,$FpPos2);
		@sitePos=();
		$FpPos1=$FpPos2;		
		$FpApart=random_exponential(1,$sublamdaFP);
		$FpPos2=int($FpPos1+$FpApart*1000);
	}
	my @sort_Pos=(sort {$a<=>$b} @$mol_FP);
	return @sort_Pos;
}

1;
