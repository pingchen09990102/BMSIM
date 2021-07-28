#!/usr/bin/env perl
package Digest;
use strict;
use warnings;

#Single enzyme solution
sub Single 
{
	my ($plus_strand1_single,$minus_strand1_single,$nickpos1_single,$str1_single) = @_; 
	my $i_single=0;
	my $j_single=0;
	my $position_01_single=0;
	my $position_02_single=0;
	my @Pos_01_single=();
	my $str_length_single=length ($str1_single);
	my $EPL=length($plus_strand1_single);

	while ($i_single<$str_length_single) 
	{
		$position_01_single= index($str1_single, $plus_strand1_single,$i_single)+1;
		push (@Pos_01_single,$position_01_single+$nickpos1_single);
		if ($position_01_single==0) 
		{
			$position_01_single=$str_length_single;
			pop @Pos_01_single;
		}
		$i_single=$position_01_single+$EPL-1;
	}
	while ($j_single<$str_length_single) 
	{		
		$position_02_single= index($str1_single, $minus_strand1_single,$j_single)+1;
		push (@Pos_01_single,$position_02_single+$EPL-2-$nickpos1_single);
		if ($position_02_single==0) 
		{
			$position_02_single=$str_length_single;
			pop @Pos_01_single;
		}
		$j_single=$position_02_single+$EPL-1;
	}
	return @Pos_01_single;
}


#Dual enzyme solution
sub Dual 
{
	my ($plus_strand1,$minus_strand1,$nickpos1,$plus_strand2,$minus_strand2,$nickpos2,$str2) = @_; 
		
	my $str_length=length ($str2);
	my @Pos_01=();
	my @Pos_02=();

	#argements of enzyme I
	my $i_1=0;
	my $j_1=0;
	my $position_1_01=0;
	my $position_1_02=0;
	my $e1PL=length($plus_strand1);
	print "e1PL=$e1PL\n";
	
	while ($i_1<$str_length) 
	{
		$position_1_01= index($str2, $plus_strand1,$i_1)+1;
		push (@Pos_01,$position_1_01+$nickpos1);
		if ($position_1_01==0) 
		{
			$position_1_01=$str_length;
			pop @Pos_01;
		}
		$i_1=$position_1_01+$e1PL-1;
	}
	while ($j_1<$str_length) 
	{
		$position_1_02= index($str2, $minus_strand1,$j_1)+1;
		push (@Pos_01,$position_1_02-2-$nickpos1);
		if ($position_1_02==0) 
		{
			$position_1_02=$str_length;
			pop @Pos_01;
		}
		$j_1=$position_1_02+$e1PL-1;
	}

	#argements of enzyme II
	my $i_2=0;
	my $j_2=0;
	my $position_2_01=0;
	my $position_2_02=0;
	my $e2PL=length($plus_strand2);
	print "e2PL=$e2PL\n";
	while ($i_2<$str_length) 
	{
		$position_2_01=index($str2,$plus_strand2,$i_2)+1;
		push (@Pos_02,$position_2_01+$nickpos2);
		if ($position_2_01==0) 
		{
			$position_2_01=$str_length;
			pop @Pos_02;
		}
		$i_2=$position_2_01+$e2PL-1;
	}
	while ($j_2<$str_length) 
	{
		$position_2_02=index($str2,$minus_strand2,$j_2)+1;
		push (@Pos_02,$position_2_02-2-$nickpos2);
		if ($position_2_02==0) 
		{
			$position_2_02=$str_length;
			pop @Pos_02;
		}
		$j_2=$position_2_02+$e2PL-1;
	}


	#Combined total enzyme loci
	push (@Pos_01,@Pos_02);
	return @Pos_01;
}

1;
