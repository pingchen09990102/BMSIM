#!/usr/bin/perl
package ErrorMOD;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);#链接所有随机数模块

sub FN
	{
		my ($p_FNsub,$mol_pre_FN) = @_; 
		my @mol_FN=();

		###由于是先标记酶切位点的，所以先加入FN、FP再考虑拉伸
		####制造FN假阴性位点，Binomial(p)，利用rbinom(n,1,0.8)产生n个随机数（n为每条分子的总酶切位点),将值为0处的点miss掉
		#输入长度分布参数，个数，均值，方差，产生长度随机数
		my $n_FN=scalar(@$mol_pre_FN)-1;#酶切位点个数
		#将位点miss掉
		#my $mol_len=$n_FN+1;
		for (my $t=0;$t<$n_FN;$t++) 
		{
			my $rbinom=random_binomial(1,1,$p_FNsub);
			if ($rbinom == 1) 
				{
					push (@mol_FN,$$mol_pre_FN[$t]);
				}

		}
		push (@mol_FN,$$mol_pre_FN[-1]);#把最后一个酶切位点放入数组
		my $mol_pre_FN_len=@$mol_pre_FN;
		my $mol_FN=@mol_FN;
		#print "FN0=$mol_pre_FN_len\tFN1=$mol_FN\n";
		#foreach my $mol_pre_FN (@$mol_pre_FN) 
		#	{
		#		print "mol_pre\t$mol_pre_FN\t";
		#	}
		#	print "\n";
		#foreach my $mol_FN (@mol_FN) 
		#	{
		#		print "mol\t$mol_FN\t";
		#	}
		#	print "\n";
		
		return @mol_FN;#注意返回的数组中酶切位点总个数和id号是不对的
	}

sub FP
	{
		my ($sublamdaFP,$mol_FP) = @_; 

		#产生插入位点间距
		my $FpApart=0;#第一个FP位点
		my $FpPos1=0;
		my $FpPos2=0;
		my @sitePos=();
		my @insert=();

		#产生插入位点间距
		$FpApart=random_exponential(1,$sublamdaFP);#第一个FP位点间距，单位为kb，后面要转化为bp
		$FpPos2=int($FpPos1+$FpApart*1000);#此时的插入位点为起始点加上间距(单位)
		my $mol_len_FP=$$mol_FP[-1];#分子长度
		my $Numsites_FP=scalar(@$mol_FP)-1;#酶切位点个数
		#print "FP_before\t";
		#foreach my $FP_before (@$mol_FP) 
		#	{
				#print "$FP_before\t";
		#	}
		#print "\n";
		while ($FpPos2 < $mol_len_FP) #如果插入位置小于总长
			{
				while (grep{$_ == $FpPos2} @$mol_FP) #去掉与原始点重叠的点
					{
						$FpApart=random_exponential(1,$sublamdaFP);#重新产生一个FP位点
						$FpPos2=int($FpPos1+$FpApart*1000);
					}
				$Numsites_FP++;#酶切位点个数加1
				#插入位点FP
				push (@$mol_FP,$FpPos2);
				#寻找插入位置下标
				#@insert=grep{$sitePos[$_]>$FpPos2} 0..$#sitePos;########修改此处减小内存和时间###############
				#splice(@$mol_FP,$insert[0],0,$FpPos2);#插入FP位点
				@sitePos=();
				$FpPos1=$FpPos2;
		
				#再产生FP位点
				$FpApart=random_exponential(1,$sublamdaFP);#第一个FP位点间距
				$FpPos2=int($FpPos1+$FpApart*1000);#此时的插入位点为起始点加上间距
		
			}

		#将位点按大小重新排序
		my @sort_Pos=(sort {$a<=>$b} @$mol_FP);#对酶切位点进行排序
		#print "FP_after\t";
		foreach my $FP_after (@sort_Pos) 
			{
				#print "$FP_after\t";
			}
		#print "\n";
		return @sort_Pos;
	}


sub Stretch
	{
		my ($sub_m_str,$sub_sd_str,$mol_Str) = @_; 
		###分子拉伸，根据rnorm(n,m=1.0,sd=0.02)产生n个分子拉伸百分比随机数，将每条分子的酶切位点位置都乘以拉伸比例
		#my @SitePos=();
		#分子总长有拉伸，因为荧光标记位点不均
		#my @rnorm=random_normal(1, $sub_m_str, $sub_sd_str);#产生一个拉伸百分数
		my $mol_Str_len=scalar(@$mol_Str);
		#每个酶切位点间距有拉伸，因为布朗运动
		my @sub_pos_str=();
		push (@sub_pos_str,$$mol_Str[0]);
		my $posPre=$$mol_Str[0];
		for (my $s=1;$s<$mol_Str_len;$s++) 
			{
				my $rnorm=random_normal(1, $sub_m_str, $sub_sd_str);#产生一个拉伸百分数
				$posPre=$posPre+($$mol_Str[$s]-$$mol_Str[$s-1])*$rnorm;
				push (@sub_pos_str,$posPre);
			}

		#my @sub_pos_str = map {$_ * $rnorm[0]} @$mol_Str;#对每个位点进行拉伸
		#print "Stretch_before\t@$mol_Str\n";
		#print "Stretch_after\t@sub_pos_str\n";
		return @sub_pos_str;

	}

1;
