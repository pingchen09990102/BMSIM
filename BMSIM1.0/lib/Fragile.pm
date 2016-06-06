#!/usr/bin/perl
package Fragile;

use strict;
use warnings;
use Math::Random qw(:all);#链接所有随机数模块
#use constant PI => 3.1415926536;
#use constant SIGNIFICANT => 5; # number of significant digits to be returned

sub FragileSite 
	{
		#传入<250bp易碎点均值和方差
		#my ($mean1,$sd1,$mean2,$sd2,$PosSet,$DistanceSet,$OriSet) = @_; #分辨率均值和方差,FragileSite位置，间距
		my ($a,$b,$PosSet,$DistanceSet) = @_; #指数分布两个参数,FragileSite位置，间距

		my $PosSet_length=@$PosSet;#易碎点位置
		#my $DistanceSet_length=@$DistanceSet;#易碎点间距
		my @FragilePos=(); 
		my ($p,$rbinom,$x)=0;
		for (my $l=0;$l<$PosSet_length;$l++) 
		{
			#先运行易碎点计算perl
			#判断易碎点类型
			#if ($$OriSet[$l] eq "+") #type I fragile
			#{
				#完全断裂
				#$p=0;

				#完全不断裂
				#$p=1;

				#累积正态断裂法
				#$Nlabel_size=($$DistanceSet[$l]-$mean1)/$sd1;#先将易碎点间距标准化
				#$p=1-&uprob ($Nlabel_size);#根据距离计算易碎点不断的概率，长度越长不断即成功概率越大

				#韦伯分布断裂法
				#$x=$Nlabel_size/100;#将长度除以100使其小于10
				#$p=$a*$b*$x**($b-1)*exp(-$a*$x**$b)#根据weibull公式计算断裂概率，a,b分别为位置参数和尺度参数
				#my $p=$a*$b*$x**($b-1)*exp(-$a*$x**$b);#根据weibull公式计算断裂概率，a,b分别为位置参数和尺度参数

				#指数分布断裂法得到成功概率
				$x=$$DistanceSet[$l];
				$p=1-$a*exp($x*$b);#根据指数分布公式计算断裂概率，1-断裂概率得成功概率用于白努力实验
				print "x=$x;p=$p\n";

				#进行白努力实验
				#$rbinom=1;
				$rbinom=random_binomial(1,1,$p);#根据前面算的概率进行白努力实验，
			#}
			#else #type II fragile
			#{
				#累积正态断裂法
				#$Nlabel_size=($$DistanceSet[$l]-$mean2)/$sd2;#先将易碎点间距标准化
				#$p=1-&uprob ($Nlabel_size);#根据距离计算易碎点不断的概率，长度越长不断即成功概率越大
				
				#完全断裂
				#$p=0;

				#$rbinom=1;
				#$rbinom=random_binomial(1,1,$p);#根据前面算的概率进行白努力实验，
			#}
			if ($rbinom==1)#为1的时候表示不断，此处不加X
			{
				next;
			}
			else #不为1的时候表示断了，调用替换脚本对应位点处碱基序列变为X
			{
				push (@FragilePos,$$PosSet[$l]);#把断裂位置存入数组
			}
		
		}
		return (@FragilePos);#返回断裂位置数组

	}#对应sub



1
