#!/usr/bin/perl
package RES;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned

#########keep some labels with intervals under 1kb or between 1kb to 1.5kb with user defined propotions, the positon of other site with interval <1.5kb take the mean  averaged under 1500 bp###############
sub RES1 
{
	
	my ($subunder1000_che,$subunder1500_che,$subpos_str) = @_; 
	my $dataset_length=@$subpos_str;
	my $confused_num=0;
	my $confused_length=0;
	my @Numsites_new=();
	my @label_size=();

	for (my $u=0;$u<$dataset_length-1;$u++) 
	{
		my $label_size=$$subpos_str[$u+1]-$$subpos_str[$u]+1;
		push (@label_size,$label_size);
	}

	my @Size1000_pre=grep{$label_size[$_]<=1000} 0..$#label_size;
	my @Size1500_pre=grep{$label_size[$_]>1000 && $label_size[$_]<=1500} 0..$#label_size;
	my $S1000_len=scalar(@Size1000_pre);
	my $S1500_len=scalar(@Size1500_pre);
	my $S1000_remain=int($S1000_len*$subunder1000_che);
	my $S1500_remain=int($S1500_len*$subunder1500_che);
	my $high1000=0;
	my $high1500=0;
	my @S1000_remain=();
	my @S1500_remain=();
	if ($S1000_len==0) 
	{
			@S1000_remain=();
		}
	else
		{
			$high1000=$S1000_len-1;
			@S1000_remain=random_uniform_integer($S1000_remain, 0, $high1000);

		}

	if ($S1500_len==0) 
		{
			@S1500_remain=();
		}
	else
		{
			$high1500=$S1500_len-1;
			@S1500_remain=random_uniform_integer($S1500_remain, 0, $high1500);
		}
	
	my @merge=(); 

	#Looking for the label intervals <1kb, add to @merge
	for (my $j=0;$j<$S1000_len;$j++) 
	{
		if (grep{$j eq $_} @S1000_remain) 
		{
			next;
		}
		else 
		{
			push (@merge,$Size1000_pre[$j]);
		}
	}

	#Looking for the label intervals 1~1.5kb, add to @merge
	for (my $k=0;$k<$S1500_len;$k++) 
	{
		if (grep{$k eq $_} @S1500_remain) 
		{
			next;
		}
		else 
		{
			push (@merge,$Size1500_pre[$k]);
		}
	}
	
	my @PosNew=();
	my @MerPos0=();
	my ($MerStr,$MerEnd,$MerLen,$MerPos,$v1,$v2,$MerPos0);
	for (my $v=0;$v<$dataset_length;$v++) 
	{
		if (grep{$_ eq $v} @merge) 
		{
			$v1=$v+1;
			while (grep{$_ eq $v1} @merge) 
			{
				$v2=$v1+1;
				$v1=$v2;
			}
			$MerStr=$v;
			$MerEnd=$v1;
			$MerLen=$MerEnd-$MerStr+1;
			@MerPos0=@$subpos_str[eval($MerStr)..eval($MerEnd)];
			$MerPos0=sum(@MerPos0);
			$MerPos=$MerPos0/$MerLen;
			push (@PosNew,$MerPos);
	
			$v=$v1;
		}
		else 
		{
			$MerPos=$$subpos_str[$v];
			push (@PosNew,$MerPos);
		}
	}
	
	return (@PosNew);
}

########Resolution model is Gaussian. The likelihood of resolving in general is modeled by the cumulative Gaussian##########
sub RES2 
{
	my ($mean,$sd,$subpos_str) = @_; 
	my $dataset_length=@$subpos_str;#total number of labels=$dataset_length-1
	my $confused_num=0;
	my $confused_length=0;
	my @Numsites_new=();
	my @label_size=();
	my @merge=(); 

	for (my $u=0;$u<$dataset_length-1;$u++) 
	{
		my $label_size=$$subpos_str[$u+1]-$$subpos_str[$u]+1;
		push (@label_size,$label_size);
	}

	my $label_size_len=@label_size;
	for (my $l=0;$l<$label_size_len;$l++) 
	{
		my $Nlabel_size=($label_size[$l]/1000-$mean)/$sd;
		my $p=1-&uprob ($Nlabel_size);
		my $rbinom=random_binomial(1,1,$p);
		if ($rbinom==1)#when 1 means to discern
		{
			next;
		}
		else 
		{
			push (@merge,$l);
		}
	}

	my @PosNew=();
	my @MerPos0=();
	my ($MerStr,$MerEnd,$MerLen,$MerPos,$v1,$v2,$MerPos0);
	for (my $v=0;$v<$dataset_length;$v++) 
	{
		if (grep{$_ eq $v} @merge) 
		{
			$v1=$v+1;
			while (grep{$_ eq $v1} @merge) 
			{
				$v2=$v1+1;
				$v1=$v2;
			}
			$MerStr=$v;
			$MerEnd=$v1;
			$MerLen=$MerEnd-$MerStr+1;
			@MerPos0=@$subpos_str[eval($MerStr)..eval($MerEnd)];
			$MerPos0=sum(@MerPos0);
			$MerPos=$MerPos0/$MerLen;
			push (@PosNew,$MerPos);
	
			$v=$v1;
		}
		else 
		{
			$MerPos=$$subpos_str[$v];
			push (@PosNew,$MerPos);
		}
	}

	return (@PosNew);
}

sub uprob 
{ # Upper probability   N(0,1^2)
	my ($x) = @_;
	return &precision_string(&subuprob($x));
}

sub subuprob 
{
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) 
	{
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} 
	elsif ($absx <= 100) 
	{
		for (my $i = 18; $i >= 1; $i--) 
		{
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}

sub log10 
{
	my $n = shift;
	return log($n) / log(10);
}

sub precision 
{
	my ($x) = @_;
	return abs int(&log10(abs $x) - SIGNIFICANT);
}

sub precision_string 
{
	my ($x) = @_;
	if ($x) 
	{
		return sprintf "%." . &precision($x) . "f", $x;
	} 
	else 
	{
		return "0";
	}
}

1;
