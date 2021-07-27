#!/usr/bin/env perl
package ErrorFN;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);

########Add FNs and FPs first, then add stretch bias########
########Add FNS, Binomial(p), rbinom(n,1,0.8), n is total labels, delate labels with value 0######## 	
sub FN
{
	my ($p_FNsub,$mol_pre_FN) = @_; 
	my @mol_FN=();
	my $n_FN=scalar(@$mol_pre_FN)-1;#number of labels
	for (my $t=0;$t<$n_FN;$t++) 
	{
		my $rbinom=random_binomial(1,1,$p_FNsub);
		if ($rbinom == 1) 
		{
			push (@mol_FN,$$mol_pre_FN[$t]);
		}

	}
	push (@mol_FN,$$mol_pre_FN[-1]);#add the last label
	my $mol_pre_FN_len=@$mol_pre_FN;
	my $mol_FN=@mol_FN;
		
	return @mol_FN;
}

1;
