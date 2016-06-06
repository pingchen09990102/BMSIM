#!/usr/bin/perl
package RES;


use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);#链接所有随机数模块
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned


sub RES1 
	{
	
		my ($subunder1000_che,$subunder1500_che,$subpos_str) = @_; 

		#########分辨率：保留10%的1kb以下的酶切位点，80%的1kb-1.5kb的位点，其他1500bp以下的位点位置取均值（相加/2）
		my $dataset_length=@$subpos_str;#总酶切位点数=$dataset_length-1
		my $confused_num=0;
		my $confused_length=0;
		my @Numsites_new=();
		my @label_size=();

		#将所有酶切位点间距存入数组
		for (my $u=0;$u<$dataset_length-1;$u++) 
		{
			my $label_size=$$subpos_str[$u+1]-$$subpos_str[$u]+1;
			push (@label_size,$label_size);
		}
		#print "RES_before\t@$subpos_str\n";


		#统计不同间距的酶切位点
		my @Size1000_pre=grep{$label_size[$_]<=1000} 0..$#label_size;#将小于1k的间距下标存入数组
		my @Size1500_pre=grep{$label_size[$_]>1000 && $label_size[$_]<=1500} 0..$#label_size;#将大于1k小于1.5k的间距下标存入数组
		my $S1000_len=scalar(@Size1000_pre);
		my $S1500_len=scalar(@Size1500_pre);
		my $S1000_remain=int($S1000_len*$subunder1000_che);#$under1000由参数文件获得
		my $S1500_remain=int($S1500_len*$subunder1500_che);#$under1500由参数文件获得
		my $high1000=0;
		my $high1500=0;
		my @S1000_remain=();
		my @S1500_remain=();
		#产生保留的酶切位点随机数
		if ($S1000_len==0) #如果没有小于1k的点
			{
				@S1000_remain=();
			}
		else#有小于1k的点
			{
				$high1000=$S1000_len-1;
				@S1000_remain=random_uniform_integer($S1000_remain, 0, $high1000);

			}

		if ($S1500_len==0) #如果没有大于1k小于1.5k的点
			{
				@S1500_remain=();
			}
		else#如果有大于1k小于1.5k的点
			{
				$high1500=$S1500_len-1;
				@S1500_remain=random_uniform_integer($S1500_remain, 0, $high1500);
			}
		
		my @merge=(); 
		#寻找要保留的近邻点<1kb的id，存入@merge
		#my $Size1000_pre_len=@Size1000_pre;
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

		#寻找要保留的近邻点大于1kb小于1.5kb的id，存入@merge
		#my $Size1500_pre_len=@Size1500_pre;
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
		
		#合并位点，取均值
		my @PosNew=();
		my @MerPos0=();
		my ($MerStr,$MerEnd,$MerLen,$MerPos,$v1,$v2,$MerPos0);
		for (my $v=0;$v<$dataset_length;$v++) 
		{
			if (grep{$_ eq $v} @merge) 
			{
				$v1=$v+1;#检测下一个点是否在merge中
				while (grep{$_ eq $v1} @merge) #merge里面存的是数组坐标
				{
					$v2=$v1+1;
					$v1=$v2;#再检测下一个位点是否在merge中
				}
				$MerStr=$v;
				$MerEnd=$v1;
				$MerLen=$MerEnd-$MerStr+1;
				@MerPos0=@$subpos_str[eval($MerStr)..eval($MerEnd)];
				$MerPos0=sum(@MerPos0);
				$MerPos=$MerPos0/$MerLen;
				push (@PosNew,$MerPos);
		
				$v=$v1;
				#if ($v1==$dataset_length-1) #如果为最后一个点，则跳出循环
				#	{
						#push (@PosNew,$$subpos_str[$v1]);#把最后一个点加入进数组
				#		last;
				#	}
			}
			else 
			{
				$MerPos=$$subpos_str[$v];#不用合并，
				push (@PosNew,$MerPos);
			}
		}
		#print "RES_after\t@PosNew\n";

		#return ($return,$subM_ID);
		return (@PosNew);
	}#对应sub


sub RES2 
	{
		#修改传入参数，部分参数可去除
		my ($mean,$sd,$subpos_str) = @_; #分辨率均值和方差

		#########分辨率：保留10%的1kb以下的酶切位点，80%的1kb-1.5kb的位点，其他1500bp以下的位点位置取均值（相加/2）
		my $dataset_length=@$subpos_str;#总酶切位点数=$dataset_length-1
		my $confused_num=0;
		my $confused_length=0;
		my @Numsites_new=();
		my @label_size=();
		my @merge=(); 

		#将所有酶切位点间距存入数组
		for (my $u=0;$u<$dataset_length-1;$u++) 
		{
			my $label_size=$$subpos_str[$u+1]-$$subpos_str[$u]+1;
			push (@label_size,$label_size);
		}
		#print "RES_before\t@$subpos_str\n";

		#寻找要保留的近邻点大于1kb小于1.5kb的id，存入@merge
		my $label_size_len=@label_size;
		for (my $l=0;$l<$label_size_len;$l++) 
		{
			#数据先标准化
			my $Nlabel_size=($label_size[$l]/1000-$mean)/$sd;#先将长度转化为kb
			my $p=1-&uprob ($Nlabel_size);#根据正态分布概率密度函数计算概率密度,长度越长成功概率越大
			#print "label_size=$label_size[$l]\tp=$p\n";
			my $rbinom=random_binomial(1,1,$p);#根据前面算的概率进行白努力实验，
			if ($rbinom==1)#为1的时候表示能分辨 
			{
				next;
			}
			else 
			{
				push (@merge,$l);#把下标存入数组
			}
		}
		
		#合并位点，取均值
		my @PosNew=();
		my @MerPos0=();
		my ($MerStr,$MerEnd,$MerLen,$MerPos,$v1,$v2,$MerPos0);
		for (my $v=0;$v<$dataset_length;$v++) 
		{
			if (grep{$_ eq $v} @merge) 
			{
				$v1=$v+1;#检测下一个点是否在merge中
				while (grep{$_ eq $v1} @merge) #merge里面存的是数组坐标
				{
					$v2=$v1+1;
					$v1=$v2;#再检测下一个位点是否在merge中
				}
				$MerStr=$v;
				$MerEnd=$v1;
				$MerLen=$MerEnd-$MerStr+1;
				@MerPos0=@$subpos_str[eval($MerStr)..eval($MerEnd)];
				$MerPos0=sum(@MerPos0);
				$MerPos=$MerPos0/$MerLen;
				push (@PosNew,$MerPos);
		
				$v=$v1;
				#if ($v1==$dataset_length-1) #如果为最后一个点，则跳出循环
				#	{
						#push (@PosNew,$$subpos_str[$v1]);#把最后一个点加入进数组
				#		last;
				#	}
			}
			else 
			{
				$MerPos=$$subpos_str[$v];#不用合并，
				push (@PosNew,$MerPos);
			}
		}
		#print "RES_after\t@PosNew\n";

		return (@PosNew);

	}#对应sub

sub uprob { # Upper probability   N(0,1^2)
	my ($x) = @_;
	return &precision_string(&subuprob($x));
}

sub subuprob {
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
		} elsif ($absx <= 100) 
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

sub log10 {
	my $n = shift;
	return log($n) / log(10);
}

sub precision {
	my ($x) = @_;
	return abs int(&log10(abs $x) - SIGNIFICANT);
}

sub precision_string {
	my ($x) = @_;
	if ($x) {
		return sprintf "%." . &precision($x) . "f", $x;
	} else {
		return "0";
	}
}

1
