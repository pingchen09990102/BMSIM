#!/usr/bin/perl
package Digest;

use strict;
use warnings;

#Single enzyme solution
sub Single 
{
	#@Pos_all=Digest::Single($plus1,$minus1,$nickpos1,$string);
	my ($plus_strand1_single,$minus_strand1_single,$nickpos1_single,$str1_single) = @_; 

		#初始化Nt.BspQI的参数
		my $i_single=0;
		my $j_single=0;
		my $position_01_single=0;
		my $position_02_single=0;
		my @Pos_01_single=();
		my $str_length_single=length ($str1_single);#以bp为单位
		my $EPL=length($plus_strand1_single);#计算酶切pattern里的碱基数

		#单酶切
		while ($i_single<$str_length_single) 
		{
			$position_01_single= index($str1_single, $plus_strand1_single,$i_single)+1;
			push (@Pos_01_single,$position_01_single+$nickpos1_single);#正链切割要加上一些碱基,如NtBspQI为7，$nickpos1_single
			if ($position_01_single==0) 
			{
				$position_01_single=$str_length_single;
				pop @Pos_01_single;
			}
			$i_single=$position_01_single+$EPL-1;
			#my $plus_strand1_single_len=length($plus_strand1_single);#计算酶切pattern里的碱基数
			#$i_single=$position_01_single+$plus_strand1_single_len-1;
			#my $plus_strand1_single_len_1=$plus_strand1_single_len-1;
			#print "plus_strand1_single_len_1=$plus_strand1_single_len_1\n";
		}
		while ($j_single<$str_length_single) 
		{
			
			$position_02_single= index($str1_single, $minus_strand1_single,$j_single)+1;
			push (@Pos_01_single,$position_02_single+$EPL-2-$nickpos1_single);#反链切割要减去一些碱基
			if ($position_02_single==0) 
			{
				$position_02_single=$str_length_single;
				pop @Pos_01_single;
			}
			$j_single=$position_02_single+$EPL-1;
			#my $minus_strand1_single_len=length($minus_strand1_single);#计算酶切pattern里的碱基数
			#$j_single=$position_02_single+$minus_strand1_single_len-1;
			#my $minus_strand1_single_len_1=$minus_strand1_single_len-1;
			#print "minus_strand1_single_len_1=$minus_strand1_single_len_1\n";
		}
	return @Pos_01_single;

}


#Dual enzyme solution
sub Dual 
{
	#@Pos_all=Digest::Dual($plus1,$minus1,$np1,$plus2,$minus2,$np2,$string);
	my ($plus_strand1,$minus_strand1,$nickpos1,$plus_strand2,$minus_strand2,$nickpos2,$str2) = @_; 
		
		my $str_length=length ($str2);#以bp为单位
		my @Pos_01=();
		my @Pos_02=();

		#酶1的参数
		my $i_1=0;
		my $j_1=0;
		my $position_1_01=0;
		my $position_1_02=0;
		my $e1PL=length($plus_strand1);#计算酶切pattern里的碱基数
		print "e1PL=$e1PL\n";
	
		#Nt.BspQI切
		while ($i_1<$str_length) 
		{
			$position_1_01= index($str2, $plus_strand1,$i_1)+1;
			push (@Pos_01,$position_1_01+$nickpos1);#正链切割要加上一些碱基,如NtBspQI为7，$nickpos1_single
			if ($position_1_01==0) 
				{
					$position_1_01=$str_length;
					pop @Pos_01;
				}
			#my $plus_strand1_len=length($plus_strand1);
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
			#my $minus_strand1_len=length($minus_strand1);
			$j_1=$position_1_02+$e1PL-1;
		}

	#酶2的参数
	my $i_2=0;
	my $j_2=0;
	my $position_2_01=0;
	my $position_2_02=0;
	my $e2PL=length($plus_strand2);#计算酶切pattern里的碱基数
	print "$e2PL=$e2PL\n";

	#Nb.BtvCI切
	while ($i_2<$str_length) 
	{
		$position_2_01=index($str2,$plus_strand2,$i_2)+1;
		push (@Pos_02,$position_2_01+$nickpos2);
		if ($position_2_01==0) 
			{
				$position_2_01=$str_length;
				pop @Pos_02;
			}
		#my $plus_strand2_len=length($plus_strand2);
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
		#my $minus_strand2_len=length($minus_strand2);
		$j_2=$position_2_02+$e2PL-1;
	}


	#合计总的酶切位点
	push (@Pos_01,@Pos_02);#将两种酶切位点合并到@Pos_01
	
	return @Pos_01;

}

1;