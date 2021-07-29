#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);
use Math::Random::MT qw(srand rand);
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Chimera;
use ErrorFN;
use ErrorFP;
use ErrorStretch;
use Digest;
use RES;
use Fragile;

my $start_time=time();

my ($out1,$out2,$Cov,$chr_fa,$bnx,$FragileIN,$str,$RtI,$RtII,$snr,$Ints,$out3,$out4,$out5,$out6,$out7,$e,$np1,$np2,$FragileArg);
my $man = 0;
my $help = 0;
my $version = 0;
my $EXPav = 70;
my $p_FN = 0.85;
my $lamdaFP = 100;
my $res=2;
my $Chi_perc=0.1;
my $Chi_Bi=0.85;
my $Chi_Tri=0.13;
my $Chi_Qua=0.02;
my $ProName="out";

GetOptions (
		'h|help' => \$help,
		'v|version' => \$version,
		'man' => \$man,
		'cov|Coverage:i' => \$Cov,
		'ca|chr_fa:s' => \$chr_fa,#fa input
		'bnx|bnx_outfile:s' => \$bnx,#bnx output
		'fragile|Fragile_infile:s' => \$FragileIN,#Fragile Input
		'lm|length_mean:f' => \$EXPav,#possion
		'FNp|FN_probability:f' => \$p_FN,#add FN, enzyme nick effiency
		'FPi|FP_intensity:f' => \$lamdaFP,#add FP, possion
		'str|stretch_normal:s' => \$str,#stretch~N(ave,std),comma-separated lists of values
		'Rt|Res_type:i' => \$res,
		'RtI|Res_typeI:s' => \$RtI,#type I, under1000,under1500
		'RtII|Res_typeII:s' => \$RtII,#type II, resolution~N(ave,std)
		'snr|SNR:s' => \$snr,#type II, SNR~N(ave,std)
		'Ints|Intensity:s' => \$Ints,#type II, Intensity~N(ave,std)
		'Chi_perc|Chimeric_proportion:f' => \$Chi_perc,
		'Chi_Bi|Bimera_proportion:f' => \$Chi_Bi,
		'Chi_Tri|Trimera_proportion:f' => \$Chi_Tri,
		'Chi_Qua|Quadramera_proportion:f' => \$Chi_Qua,
		'e|enzymes_pattern:s' => \$e,#enzyme recognized pattern split with comma
		'np1|nick_position1:i' => \$np1,#5' to 3', nicking position of the enzyme pattern for enzyme 1
		'np2|nick_position2:i' => \$np2,#5' to 3', nicking position of the enzyme pattern for enzyme 2
		'p|proj:s' => \$ProName,#project name for the simulation
		'f|FragileArg:s' => \$FragileArg
              )  
	      
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
if ($version)
{
    print "BMSIM Version 1.0.0\n";
    exit;
}
#if no new argument input(split with comma), keep the argument default
my ($m_str,$sd_str,$under1000,$under1500,$ResMean,$ResSdt,$QX11av,$QX11sd,$QX12av,$QX12sd,$enzyme_len,$plus1,$minus1,$plus2,$minus2,$FragilePa,$FragilePb);

#default stretch  arguments
my @Str_arguments=(0.98,0.08);
if ($str)
{
    @Str_arguments = split(/,/,$str);
	$m_str=$Str_arguments[0];
	$sd_str=$Str_arguments[1];
}
else
{
	$m_str=$Str_arguments[0];
	$sd_str=$Str_arguments[1];
}

#default Resolution model type I
my @RtI_arguments=(0,0.9);
if ($RtI)
{
    @RtI_arguments = split(/,/,$RtI);
	$under1000=$RtI_arguments[0];
	$under1500=$RtI_arguments[1];
}
else
{
	$under1000=$RtI_arguments[0];
	$under1500=$RtI_arguments[1];
}

#default Resolution model type II:resolution~N(ave,std)
my @RtII_arguments=(1.2,0.9);
if ($RtII)
{
    @RtII_arguments = split(/,/,$RtII);
	$ResMean=$RtII_arguments[0];
	$ResSdt=$RtII_arguments[1];
}
else
{
	$ResMean=$RtII_arguments[0];
	$ResSdt=$RtII_arguments[1];
}

#default SNR
my @snr_arguments=(3,0.66);
if ($snr)
{
    @snr_arguments = split(/,/,$snr);
	$QX11av=$snr_arguments[0];
	$QX11sd=$snr_arguments[1];
}
else
{
	$QX11av=$snr_arguments[0];
	$QX11sd=$snr_arguments[1];
}

#default intensity
my @Ints_arguments=(1,0.2);
if ($Ints)
{
    @Ints_arguments = split(/,/,$Ints);
	$QX12av=$Ints_arguments[0];
	$QX12sd=$Ints_arguments[1];
}
else
{
	$QX12av=$Ints_arguments[0];
	$QX12sd=$Ints_arguments[1];
}

#default enzymes_pattern
my @e_arguments=("GCTCTTC","GAAGAGC");
if ($e)#single 
{
    @e_arguments = split(/,/,$e);
	$enzyme_len=@e_arguments;
	print "$enzyme_len\n";
	if ($enzyme_len==2) #single
		{
			$plus1=$e_arguments[0];
			$minus1=$e_arguments[1];
			print "$plus1\t$minus1\t$np1\n";
		}
	else #dual
		{
			die "Option -np1 or --nick_position1 must be specified with new enzyme1.\n" unless $np1; # report missing required variables
			die "Option -np2 or --nick_position2 must be specified with new enzyme2.\n" unless $np2; # report missing required variables
			$plus1=$e_arguments[0];
			$minus1=$e_arguments[1];
			$plus2=$e_arguments[2];
			$minus2=$e_arguments[3];
			print "$plus1\t$minus1\t$np1\n";
			print "$plus2\t$minus2\t$np2\n";
		}
}
else #default single
{
	$plus1=$e_arguments[0];
	$minus1=$e_arguments[1];
	$np1=7;
	$enzyme_len=2;
	print "$plus1\t$minus1\t$np1\n";
}


#default FragileArg a and b
my @f_arguments=(0.7758,-0.006984);
if ($FragileArg)
{
    @f_arguments = split(/,/,$FragileArg);
	$FragilePa=$f_arguments[0];
	$FragilePb=$f_arguments[1];
}
else
{
	$FragilePa=$f_arguments[0];
	$FragilePb=$f_arguments[1];
}

#check the input files
die "Option -ca or --assembly_dir not specified.\n" unless $chr_fa; # report missing required variables
die "Option -bnx or --proj not specified.\n" unless $bnx; # report missing required variables
die "Option -fragile or --genome not specified.\n" unless $FragileIN; # report missing required variables
die "Option -cov or --ref not specified.\n" unless $Cov; # report missing required variables

####BNX file######
open (BNX, ">$bnx")or die "can't open $bnx !";

###add output bnx file head###
print BNX "# BNX File Version:	1\n";
print BNX "# Label Channels:	1\n";
print BNX "# Nickase Recognition Site 1:\n";
print BNX "# Min Molecule Length (Kb):	0\n";
print BNX "# Label SNR Filter Type:	0\n";
print BNX "# Min Label SNR:	0\n";
print BNX "# rh	SourceFolder	InstrumentSerial	Time	NanoChannelPixelsPerScan	StretchFactor	BasesPerPixel	NumberofScans	ChipId	Flowcell\n";
print BNX "# Run Data	Z:/2014-03/$ProName/Detect Molecules	A001	3/12/2014 9:11:48 PM	68819821	0.85	490.646637	1	20249,11694,10/2/2013,0011926	1\n";					
print BNX "# Quality Score QX01:	SNR\n";	
print BNX "# Quality Score QX02:	Ave Intensity\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
print BNX "#0h	LabelChannel	MoleculeId	Length	AvgIntensity	SNR	NumberofLabels	OriginalMoleculeId	ScanNumber	ScanDirection	ChipId	Flowcell\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																													
print BNX "#0f	int	int	float	float	float	int	int	int	int	string	int\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																													
print BNX "#1h	LabelChannel	LabelPositions[N]\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
print BNX "#1f	int	float\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
print BNX "#2h	LabelChannel	LabelPositions[N]\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
print BNX "#2h	int	float\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
print BNX "#Qh	QualityScoreID	QualityScores[N]\n";																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																					
print BNX "#Qf	str	float\n";

srand(time|$$);

my $Mol_ID=0;
my $MoleculeID;
my @Che_pos=();
my @ids=();
my @BNXprint=();
my @pos_str=();

open (IN1, "<$FragileIN")or die "can't open $FragileIN !";
my @FragileIN = <IN1>;
my @FragileINdata=();
shift @FragileIN;#delate head
foreach my $FragileINdata(@FragileIN) 
	{
		chomp $FragileINdata;
	}
close IN1;

my @PosData=();
my @DistanceData=();
my @genomePre=();
print "enzyme_len total=$enzyme_len\n";
for (my $r=0;$r<$Cov;$r++) #cycle according to coverage
{

	#input genome .fa
	open (GenomeFa,"< $chr_fa") or die "Can't open file:$!";
	@genomePre = <GenomeFa>;
	foreach my $s (@genomePre) 
	{
		chomp $s;
	}
	close GenomeFa;

	#cheak the number of chr and cycle with the counted number
	my @ChrSet=grep {$genomePre[$_]=~/^>.+/} 0..$#genomePre;
	my $ChrNum=@ChrSet;#count the chromosome
	my $chr_len=0;
	my @XChr=();
	for (my $v=1;$v<$ChrNum+1;$v++) 
	{
		my $chrStr=$ChrSet[$v-1]+1;#start from the next line of chr headline
		my $chrEnd;
		if ($v==$ChrNum) 
		{
			$chrEnd=@genomePre-1;
		}
		else
		{
			$chrEnd=$ChrSet[$v]-1;
		}
		my @ChrGrep=@genomePre[eval($chrStr)..eval($chrEnd)];
		my $chrPre = join("",@ChrGrep);
		$chr_len=length ($chrPre);
		my @chrfragile=grep {$_=~/^$v\t.+/} @FragileIN;
		my $ppp=@chrfragile;
		@PosData=();
		@DistanceData=();
		foreach my $chrfragile(@chrfragile) 
		{
			my @ChrFragileSet=split(/\t/,$chrfragile);
			push (@PosData,$ChrFragileSet[5]);
			push (@DistanceData,$ChrFragileSet[6]);
		}

		#add fragile sites
		my @FragileSites=Fragile::FragileSite($FragilePa,$FragilePb,\@PosData,\@DistanceData);

		my $FraStart=0;
		my $FraEnd=1;
		my $FraLen=0;
		my $FraStr=0;
		my @Chr=();
		foreach my $FraSite (@FragileSites) 
		{
			$FraLen=$FraSite-$FraEnd;
			$FraStr=substr ($chrPre, $FraStart, $FraLen);
			push (@Chr,$FraStr);
			my $e_pat_Len=length($plus1);#enzyme pattern length
			print "e_pat_Len=$e_pat_Len\n";
			$FraStart=$FraSite+$e_pat_Len-1;
			$FraEnd=$FraSite+$e_pat_Len;
		}
		$FraStr=substr ($chrPre, $FraStart);
		push (@Chr,$FraStr);
		my $Chrl=@Chr;
		my $XChr = join("X",@Chr);#join fragment produced by fragile sites with X to form the chr
		push (@XChr,$XChr);#Deposited the chr with fragile sites marked with X
	}

	my $genome = join("X",@XChr);#join the chr with fragile sites marked with X, and also the different chr with X
	
	my $phase=rand(10);
	random_set_seed_from_phrase($phase);

	my $str;
	my $dataset_length=0;
	my $rand;
	my $geno_len=length($genome);
	my $Mol_l=0;
	my @Fa=();
	my $reads;
	
	while ($geno_len!=0) 
	{
		my $read_length01_random=random_exponential(1,$EXPav);#Produce a molecular length according to the exponential distribution
		my $read_length01=int($read_length01_random*1000);#kb transfer to bp
		$str=substr ($genome, 0, $read_length01);

		#check whether it is closing to the end
		if ($geno_len<=$read_length01) #approaching to the end
		{
			if ($genome=~m/X/) 
			{
				my @reads=split(/X/,$genome); 
				foreach $reads (@reads) 
				{
					$Mol_l=length($reads);
					push (@Fa,$reads);
				}
			}
			else 
			{
				$Mol_l=length($genome);
				push (@Fa,$genome);
			}
			$geno_len=0;#done with cutting
		}
		else #not the end
		{
			if ($str!~m/X/)#no X, i.e.in the same chr. and exclude fragile sites
			{
				substr($genome, 0, $read_length01)="";#delate the output reads
				$Mol_l=length($str);
				push (@Fa,$str);
			}
			if ($str=~m/X/ && $str!~m/^X/)#with X, but not with X in head, indicating not in the same chr. or contain fragile sites
			{
				substr($genome, 0, $read_length01)="";#delate the output reads
				#break in X: 1) for the breakage of fragile sites; 2) for different chr.
				my @strSplit=split(/X/,$str); 
				foreach my $strSplit (@strSplit) 
				{
					$Mol_l=length($strSplit);
					push (@Fa,$strSplit);
				}
			}
			if ($str=~m/^X/) #with X in head 
			{
				substr($str, 1, 1)="";
				substr($genome, 0, $read_length01)="";
				$Mol_l=length($str);
				push (@Fa,$str);
			}
			$geno_len=length($genome);
		}
	}
	#done with cutting for one 

	#enzyme digest
	my $fasta_length=@Fa;
	my @CMAPdataset=();
	my $CMAP;
	my $CmapID=1;
	my $string;
	my $str_length;
	my $Numsites;
	my @Pos_all=();
	my @sort_Pos=();
	for (my $t=0;$t<$fasta_length;$t++) 
	{
		$string=$Fa[$t];
		chomp $string;
		$string=uc($string);
		$str_length=length ($string);
		@sort_Pos=();
		@Pos_all=();
		$Numsites=0;
		if ($enzyme_len==2) #single
		{
			@Pos_all=Digest::Single($plus1,$minus1,$np1,$string);
		}
		else #dual enzymes digest
		{
			@Pos_all=Digest::Dual($plus1,$minus1,$np1,$plus2,$minus2,$np2,$string);
		}
		my $Pos_length=@Pos_all;
		if ($Pos_length==0) 
		{
			next;
		}
		if ($Pos_length>0) 
		{
			@sort_Pos=(sort {$a<=>$b} @Pos_all);#To sort by enzyme loci
			$Numsites=@sort_Pos;
			my $siteid=0;
			for (my $z=0;$z<=$Numsites;$z++) 
			{
				if ($z==$Numsites) 
				{
					$siteid=$z+1;
					$CMAP=join "\t",$CmapID,$str_length,$Numsites,$siteid,0,$str_length,0,1,1;
					push (@CMAPdataset,$CMAP);
				}
				else
				{
					$siteid=$z+1;
					$CMAP=join "\t",$CmapID,$str_length,$Numsites,$siteid,1,$sort_Pos[$z],0,1,1;
					push (@CMAPdataset,$CMAP);
				}
			}
		}
		$CmapID++;
	}	

	#To produce the needed chimeric molecule ID first
	my @chemericaID=();
	my $mol_Num=$CmapID-1;
	my $cycle=int($mol_Num*$Chi_perc);
	print "$mol_Num\t$cycle\n";
	my %sns;
	my $chemericaID=0;
	my $minimum=1;
	for (my $q=0;$q<$cycle;$q++) 
		{
			$chemericaID=int(rand($mol_Num))+$minimum;
			while(defined $sns{$chemericaID})
				{
					$chemericaID=int(rand($mol_Num))+$minimum;
				}
			$sns{$chemericaID} = 1;
			push (@chemericaID,$chemericaID);
		}


####################Add error models#######################
	for (my $id=1;$id<$CmapID;$id++)
	{
		my @mol=();
		my @mol_1=();
		my @mol_pre=();
		my @pos_str=();
		my @SitePos=();
		@mol_pre=grep /^$id\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+$/, @CMAPdataset;
		foreach my $molecule (@mol_pre) 
		{
			my @data_03=split(/\t/,$molecule);
			push (@SitePos,$data_03[5]);
		}

		#add FNs
		@mol=ErrorFN::FN ($p_FN,\@SitePos);
		
		#add FPs
		@mol_1=ErrorFP::FP ($lamdaFP,\@mol);
		
		#add chimera
		if (grep {$_==$id} @chemericaID) 
			{
				my $Che_pos=join ("t",@mol_1);
				push (@Che_pos,$Che_pos);
				$Mol_ID++;
				push (@ids,$Mol_ID);
				next;
			}

		#add stretch model
		@pos_str=ErrorStretch::Stretch ($m_str,$sd_str,\@mol_1);

		#add resolution models,there is two options
		my @RES=();
		if ($res==1) 
			{
				@RES=RES::RES1($under1000,$under1500,\@pos_str);
			}
		else 
			{
				@RES=RES::RES2($ResMean,$ResSdt,\@pos_str);
			}

		@BNXprint=&PRINT($Mol_ID,$QX11av,$QX11sd,$QX12av,$QX12sd,\@RES);
		$Mol_ID=$BNXprint[1];
		print BNX "$BNXprint[0]\n";
	}
####################Add error models#######################

}

	#get all the molecules that need to be chimera
	my $isd_len=@ids;
	my $Che_pos_len=@Che_pos;
	my @Che=Chimera::Che($res,$under1000,$under1500,$Chi_Bi,$Chi_Tri,$Chi_Qua,\@ids,\@Che_pos);
	foreach my $Che (@Che) 
		{
			my @Che_dataSet=split(/\t/,$Che);
			#add stretch model 
			@pos_str=ErrorStretch::Stretch ($m_str,$sd_str,\@Che_dataSet);
		
			#add resolution model
			my @RE;
			if ($res==1) 
			{
				@RE=RES::RES1($under1000,$under1500,\@pos_str);
			}
			else
			{
				@RE=RES::RES2($ResMean,$ResSdt,\@pos_str);
			}

			@BNXprint=&PRINT($Mol_ID,$QX11av,$QX11sd,$QX12av,$QX12sd,\@RE);
			$Mol_ID=$BNXprint[1];
			print BNX "$BNXprint[0]\n";
		}
close BNX;

my $end_time = time();
my $elapsed_time = $end_time - $start_time;    
print "elapsed_time=$elapsed_time\n";


sub PRINT
{
	my ($subMol_ID,$subQX11av,$subQX11sd,$subQX12av,$subQX12sd,$PosNew) = @_; 

	#print to BNX file
	my $PosNew_len=@$PosNew;
	my ($Numsite_new,$MolSNR_total,$MolSNR,$MolIntensity_total,$MolIntensity);
	my @QX11=();
	my @QX12=();
	my @PosNew_new=();
	if ($PosNew_len<=1) 
	{
		$Numsite_new=0;#with only one nick site, which is not the true nick but the molecule length
		$MolSNR=5.0;#the average snr of molecule
		$MolIntensity=1.0;#the average intensity of molecule
	}
	else 
	{
		$Numsite_new=$PosNew_len-1;
		@QX11=random_normal($Numsite_new, $QX11av, $QX11sd);
		$MolSNR_total=sum @QX11;
		$MolSNR=$MolSNR_total/$Numsite_new;
		@QX12=random_normal($Numsite_new, $subQX12av, $subQX12sd);
		$MolIntensity_total=sum @QX12;
		$MolIntensity=$MolIntensity_total/$Numsite_new;			
	}

	$subMol_ID++;
	my $MoleculeLen=0; 
	if ($PosNew_len==0) 
	{
		$MoleculeLen=1;
	}
	else 
	{
		$MoleculeLen=sprintf "%.1f",$$PosNew[-1];
	}
	my $MolSNR_print=sprintf "%.6f",$MolSNR;
	my $MolIntensity_print=sprintf "%.3f",$MolIntensity;
	my @PosNew_print=();
	foreach my $Pos_New(@$PosNew) 
	{
		my $PosNew_print = sprintf "%.1f",$Pos_New;
		push (@PosNew_print,$PosNew_print);
	}
	my @QX11_print=();
	foreach my $QX11 (@QX11) 
	{
		my $QX11_print = sprintf "%.4f",$QX11;
		push (@QX11_print,$QX11_print);
	}
	my @QX12_print=();
	foreach my $QX12 (@QX12) 
	{
		my $QX12_print = sprintf "%.4f",$QX12;
		push (@QX12_print,$QX12_print);
	}
	@QX11_print=join("\t",@QX11_print);
	@QX12_print=join("\t",@QX12_print);
	@PosNew_print=join("\t",@PosNew_print);
			
	my $str="20249,11694,10\/2\/2013,0011926";
	my $return01=join ("\t","0",$subMol_ID,$MoleculeLen,$MolSNR_print,$MolIntensity_print,$Numsite_new,$subMol_ID,"3","-1",$str,"1");
	my $return02=join ("\t","1",@PosNew_print);
	my $return03=join ("\t","QX11",@QX11_print);
	my $return04=join ("\t","QX12",@QX12_print);
	my $return=join ("\n",$return01,$return02,$return03,$return04);
	return ($return,$subMol_ID);
}

################################################################################
##############                  Documentation                 ##################
################################################################################
##  
__END__

=head1 NAME
	
	BMSIM is a simulation tool for BioNano molecule data of the BioNano optical mapping platform
	
	Version: 1.0.0
	
=head1 USAGE
	
	perl BMSIM.pl [options]

	Example 1：single nicking: NtBspQI: 5' ...GCTCTTCN^... 3' the nick position is 7;
                                            3' ...CGAGAAGN... 5'                            
	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -fragile ./example/MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 -p EcoliSim1 -bnx out1.bnx >log.txt

	Example 2：Dual nicking  NtBspQI: 5' ...GCTCTTCN^... 3' the nick position is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the nick position is -1
	                                  3' ...CGAGAAGN... 5'                                       3' ...CGTTAC^NN... 5'
	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -fragile ./example/MG1655.txt -e GCTCTTC,GAAGAGC,GCAATG,CATTGC  -np1 7 -np2 -1 -p EcoliSim2 -bnx out2.bnx >log.txt

	Example 3：Change parameters (mean and standard deviation) for the stretch model:
	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -fragile ./example/MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 -str 0.95,0.05 -p EcoliSim3 -bnx out3.bnx >log.txt


=head1 OPTIONS

	Documentation options:
		-help	brief help message

	Required options:
		-cov,--Coverage int	Coverage depth for the simulated genome, e.g. 10, represent simulating 10X coverage depth for the simulated genome
		-ca,--chr_fa string	file with genome sequence(.fa), e.g. example/MG1655.fa
		-fragile,--Fragile_infile string	file with fragile site, e.g. example/MG1655.txt
		-bnx,--bnx_outfile string	output file(.bnx), e.g. example/EcoliSim.bnx

	Optional options:
		-lm,--length_mean int	mean length (kb) of BioNano molecule, produce BioNano molecular length according to the exponential distribution (Default:70)
		-FNp,--FN_probability float	parameter for the false negative site model, i.e. Binomil distibution, FNp is the probability of getting a sucess in one trial 
		(Default:0.85)
		-FPi,--FP_intensity float	the initial interval of false positive site (kb) (Default:100)
		-str,--stretch_normal string	parameters (ave,std) for the stretch model,i.e. Normal distribution, Stretch~N(ave,std), ave is mean, std is standard deviation, 
		comma-separated lists of values (Default:0.98,0.08)
		-Rt,--Res_type int	resolution type, "1" for "Res_typeI", "2" for "Res_typeII" (Default:2)
		-RtI,--Res_typeI string	resolution model type I:p1,p2, p1 is the proportion of sites within 1kb distance that are retained, p2 is the proportion of sites within 
		1.5kb distance that are retained, comma-separated lists of the values p1 and p2 (Default:0,0.9)
		-RtII,--Res_typeII string	parameters (ave,std) for the resolution model type II,i.e. Normal distribution, Res_typeII~N(ave,std), ave is mean, std is 
		standard deviation, comma-separated lists of values (Default:1.2,0.9)
		-snr,--SNR string	parameters (ave,std) for signal noise ratio (SNR) model of BioNano molecule,i.e. Normal distribution, SNR~N(ave,std), ave is mean, std is 
		standard deviation, comma-separated lists of values (Default:3,0.66)
		-Ints,--Intensity string	parameters (ave,std) for Intensity model of BioNano molecule, i.e. Normal distribution, Intensity~N(ave,std), ave is mean, std is 
		standard deviation, comma-separated lists of values (Default:1,0.2)
		-Chi_prop,--Chimeric_proportion	float	proportion of the total chimeric molecule in raw BioNano molecule (Default:0.1)
		-Chi_Bi,--Bimera_proportion float	proportion of Bimera (Default:0.85)
		-Chi_Tr,--Trimera_proportion	float	proportion of Trimera (Default:0.13)
		-Chi_Qua,--Quadramera_proportion float	proportion of Quadramera (Default:0.02)
		-e,--enzymes_pattern string	nicking pattern of the enzyme, e.g. GCTCTTC,GAAGAGC for NtBspQI; GCAATG,CATTGC for Nb.BsrDI; GCTCTTC,GAAGAGC,GCAATG,CATTGC for 
		dual nicking with NtBspQI and Nb.BsrDI (Default:GCTCTTC,GAAGAGC)
		-np1,--nick_position1 int	nicking position of the enzyme pattern from 5' to 3' for enzyme 1,
         		e.g. NtBspQI: 5' ...GCTCTTCN^... 3' the np1 is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the np1 is -1.
				      3' ...CGAGAAGN... 5'                             3' ...CGTTAC^NN... 5'	 
		-np2,--nick_position2 int	nicking position of the enzyme pattern from 5' to 3' for enzyme 2,
			e.g. NtBspQI: 5' ...GCTCTTCN^... 3' the np2 is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the np2 is -1.
				      3' ...CGAGAAGN... 5'                             3' ...CGTTAC^NN... 5'

		-p,--proj string	project name for the simulation, e.g. EcoliSim


=cut

