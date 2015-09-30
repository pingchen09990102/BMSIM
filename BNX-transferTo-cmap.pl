#!/usr/bin/perl
use strict;
use warnings;

die "#usage:perl $0 <in.bnx> <out.cmap>\n" unless @ARGV==2;  

my $infile_map=$ARGV[0];
open (FILE1,"< $infile_map") or die "Can't open file:$infile_map!";
my @BNXdataset_pre = <FILE1>;
chop @BNXdataset_pre;
#cmap数量等于BNX减去表头行数除以2
my @BNXdataset=grep {/^\d.+/} @BNXdataset_pre;#只抓取前面两行
my $Num_cmap=(scalar(@BNXdataset))/2;
close FILE1;

my $outfile_map=$ARGV[1];
open FILE2,"> $outfile_map";

#打印表头
print FILE2 "# hostname=bionano\n";
print FILE2 "# \$ cd /home/bionano; /home/bionano/tools/RefAligner -maxthreads 16 -o /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal1/exp_refineFinal1 -reff /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group6 -T 1e-7 -usecolor 1 -A 5 -extend 1 -MaxCov 100 -MultiMode -nosplit 2 -EndTrim 6.99 -biaswt 0 -LRbias 1e2 -deltaX 4 -deltaY 6 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -outlier 1e-5 -endoutlier 0 -Mprobeval -f -refine 3 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 1.2 -insertThreads 4 -maxmem 4 -usecolor 1 -FP 1.5 -FN 0.15 -sd 0.2 -sf 0.2 -sr 0.03 -res 3.3 -id 6 -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group1_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group2_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group3_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group4_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group5_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group6_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group7_mapped_group6.bnx -i /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/contigs/exp_refineFinal0/exp_refineFinal0_group8_mapped_group6.bnx -stdout -stderr -XmapStatRead /home/bionano/Ecoli_test/A7_11409_6244_s4_2012-08-22_17_26_2/output/molecule_stats.txt\n";
print FILE2 "# CompileDir= /home/users/vdergachev/tools/tools.3147.3265/3147 SVNversion=3265 \$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/3147/RefAligner.cpp 3243 2014-09-16 22:36:17Z tanantharaman \$\n";
print FILE2 "# FLAGS: USE_SSE=0 USE_AVX=0 USE_MIC=0 USE_PFLOAT=1 USE_RFLOAT=1 DEBUG=1 VERB=1\n";
print FILE2 "# CMAP File Version:	0.1\n";
print FILE2 "# Label Channels:	1\n";
print FILE2 "# Nickase Recognition Site 1:	unknown\n";
print FILE2 "# Number of Consensus Nanomaps:	$Num_cmap\n";
print FILE2 "#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n";
print FILE2 "#f int	float	int	int	int	float	float	float	float\n";


my $len_BNX=@BNXdataset;
for (my $i=0;$i<$len_BNX;$i=$i+2) 
	{
		my @data_01=split(/\t/,$BNXdataset[$i]);
		my @data_02=split(/\t/,$BNXdataset[$i+1]);
	
		for (my $k=1;$k<$data_01[5]+1;$k++) 
			{
				if ($data_02[$k]==0.0) 
					{
						print FILE2 "$data_01[1]\t$data_01[2]\t$data_01[5]\t$k\t1\t20.0\t0\t1\t1\n";
						next;
					}
				print FILE2 "$data_01[1]\t$data_01[2]\t$data_01[5]\t$k\t1\t$data_02[$k]\t0\t1\t1\n";
			}
			my $lastID=$data_01[5]+1;
			print FILE2 "$data_01[1]\t$data_01[2]\t$data_01[5]\t$lastID\t0\t$data_01[2]\t0\t1\t1\n";
	}
close FILE2;