#!/usr/bin/env perl
package ErrorStretch;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);


###Stretch model, rnorm(n,m=1.0,sd=0.02), each label site is multiplied by the stretch ratio#####
sub Stretch
{
	my ($sub_m_str,$sub_sd_str,$mol_Str) = @_; 
	my $mol_Str_len=scalar(@$mol_Str);
	my @sub_pos_str=();
	push (@sub_pos_str,$$mol_Str[0]);
	my $posPre=$$mol_Str[0];
	for (my $s=1;$s<$mol_Str_len;$s++) 
	{
		my $rnorm=random_normal(1, $sub_m_str, $sub_sd_str);
		$posPre=$posPre+($$mol_Str[$s]-$$mol_Str[$s-1])*$rnorm;
		push (@sub_pos_str,$posPre);
	}

	return @sub_pos_str;

}

1;

