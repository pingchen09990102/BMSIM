#!/usr/bin/env perl
package Fragile;

use strict;
use warnings;
use Math::Random qw(:all);

sub FragileSite 
{
	my ($a,$b,$PosSet,$DistanceSet) = @_; 
	my $PosSet_length=@$PosSet;#position of fragile site
	my @FragilePos=(); 
	my ($p,$rbinom,$x)=0;
	for (my $l=0;$l<$PosSet_length;$l++) 
	{
		#if ($$OriSet[$l] eq "+") #type I fragile
		#{
			#completely break
			#$p=0;

			#don't break
			#$p=1;

			#cumulative normal model
			#$Nlabel_size=($$DistanceSet[$l]-$mean1)/$sd1;
			#$p=1-&uprob ($Nlabel_size);

			#weibull model 
			#$x=$Nlabel_size/100;
			#$p=$a*$b*$x**($b-1)*exp(-$a*$x**$b)
			#my $p=$a*$b*$x**($b-1)*exp(-$a*$x**$b);

			#Exponential model
			$x=$$DistanceSet[$l];
			$p=1-$a*exp($x*$b);
			print "x=$x;p=$p\n";

			#$rbinom=1;
			$rbinom=random_binomial(1,1,$p);
		#}
		#else #type II fragile
		#{
			#cumulative normal model
			#$Nlabel_size=($$DistanceSet[$l]-$mean2)/$sd2;
			#$p=1-&uprob ($Nlabel_size);
				
			#completely Break
			#$p=0;

			#$rbinom=1;
			#$rbinom=random_binomial(1,1,$p);，
		#}
		if ($rbinom==1)#don't break when it is 1
		{
			next;
		}
		else #break when it is not 1, sign the corresponding position with X
		{
			push (@FragilePos,$$PosSet[$l]);
		}
		
	}
	return (@FragilePos);
}

1;
