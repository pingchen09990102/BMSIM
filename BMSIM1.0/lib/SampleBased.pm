#!/usr/bin/perl
package SampleBased;

use strict;
use warnings;


###lengths (up 100kb) and labeling signal score (up 阈值) 
##从用户提供的数据库中提取大于基因组长度的分子条数
	while($total<$chr_len)#检验所有分子总长是否大于染色体总长
	{
		$n_len=$n_len+100;#条数不够，参加条数
		@MolLen=random_exponential($n_len,$EXPav);#产生指数分布
		#@MolLenUp100=grep{$_>=100} @MolLen;#截取大于100kb的分子
		@MolLenUp100=grep{$_>0} @MolLen;#截取大于0kb的分子
		$sum=sum @MolLenUp100;
		my $num=@MolLenUp100;
		$total=$sum*1000;
	}

#从用户给的SNR和intensity库中抽取相应个数值，加入的RES模块中
