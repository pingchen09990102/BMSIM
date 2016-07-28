#!/usr/bin/perl
package Chimera;

use strict;
use warnings;
use List::Util qw/sum/;
use Math::Random qw(:all);
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned

sub Che 
{
my %Chime;
my ($res_che,$optRes1,$optRes2,$Chi_Bi,$Chi_Tri,$Chi_Qua,$list1,$list2)=@_; 
my $list_len=@$list1;
print "list_len=$list_len\n";
for (my $l=0;$l<$list_len;$l++) 
{
	my $che_key=$$list1[$l];
	$Chime{$che_key}=$$list2[$l];
}

my @keys=keys %Chime;
my $key_len=@keys;
my $key01;
my $key02;
my $key03;
my $key04;
my @mol_01;
my @mol_02;
my @mol_03;
my @mol_04;
my @ID=();
my $ID;
my %sns;
my $r;
my $i;
my $point1=0;
my $point2=0;
my $point3=0;
my @mol_02_new=();
my @mol_03_new=();
my @mol_04_new=();
my @pos_str=();
my @return=();
my $return;
my $key_id;
my @Chi_return;

######Bimeras######
my $Bi_len_all=int($list_len*$Chi_Bi);
print "Bi_len_all=$Bi_len_all\n";
my $Bi_len_int=int($Bi_len_all/2);
print "Bi_len_int=$Bi_len_int\n";
for ($i=0;$i<$Bi_len_int;$i++) 
{
	@ID=();
	$ID=0;
	%sns=();
	@mol_02_new=();
	for ($r=0;$r<2;$r++) 
	{
		$ID=int(rand($key_len));
		while(defined $sns{$ID})#if definde indicating, need new 
		{
			$ID=int(rand($key_len));
		}
			$sns{$ID} = 1;
			push (@ID,$ID);
	}

	$key_id=$ID[0];
	$key01=$keys[$key_id];
	$key_id=$ID[1];
	$key02=$keys[$key_id];
	@mol_01=split(/t/,$Chime{$key01});
	@mol_02=split(/t/,$Chime{$key02});
	$point1=int(rand($mol_01[-1]));
	@mol_02_new= map {$_+ $point1} @mol_02;
	push (@mol_01,@mol_02_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
	print "Chi:2\t$pos_str[-1]\n";
	push (@Chi_return,$pos_str[-1]);
	@keys=keys %Chime;
	$key_len=@keys;
	delete $Chime{$key01}; 
	delete $Chime{$key02}; 
	@keys=keys %Chime;
	$key_len=@keys;
	$return=join ("\t",@pos_str);
	push (@return,$return);
}

######Trimeras######
my $Tri_len_all=int($list_len*$Chi_Tri);
print "Tri_len_all=$Tri_len_all\n";
my $Tri_len_int=int($Tri_len_all/3);
print "Tri_len_int=$Tri_len_int\n";
for ($i=0;$i<$Tri_len_int;$i++) 
{
	@ID=();
	$ID=0;
	%sns=();
	$sns{$ID} = 0;
	@mol_02_new=();
	@mol_03_new=();
	for ($r=0;$r<3;$r++) 
	{
		$ID=int(rand($key_len));
		while(defined $sns{$ID})
		{
			$ID=int(rand($key_len));
		}
		$sns{$ID} = 1;
		push (@ID,$ID);
	}
	$key_id=$ID[0];
	$key01=$keys[$key_id];
	$key_id=$ID[1];
	$key02=$keys[$key_id];
	$key_id=$ID[2];
	$key03=$keys[$key_id];
	@mol_01=split(/t/,$Chime{$key01});
	@mol_02=split(/t/,$Chime{$key02});
	@mol_03=split(/t/,$Chime{$key03});
	$point1=int(rand($mol_01[-1]));
	@mol_02_new= map {$_+ $point1} @mol_02;
	push (@mol_01,@mol_02_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
	$point2=int(rand($mol_01[-1]));
	@mol_03_new= map {$_+ $point2} @mol_03;
	push (@mol_01,@mol_03_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
	print "Chi:3\t$pos_str[-1]\n";
	push (@Chi_return,$pos_str[-1]);
	delete $Chime{$key01}; 
	delete $Chime{$key02}; 
	delete $Chime{$key03}; 
	@keys=keys %Chime;
	$key_len=@keys;
	$return=join ("\t",@pos_str);
	push (@return,$return);
}

######Quadrameras######
my $Qua_len_all=int($list_len*$Chi_Qua);
print "Qua_len_all=$Qua_len_all\n";
my $Qua_len_int=int($Qua_len_all/4);
print "Qua_len_int=$Qua_len_int\n";
for ($i=0;$i<$Qua_len_int;$i++) 
{
	@ID=();
	$ID=0;
	%sns=();
	$sns{$ID} = 0;
	@mol_02_new=();
	@mol_03_new=();
	@mol_04_new=();
	for ($r=0;$r<4;$r++) 
	{
		$ID=int(rand($key_len));
		while(defined $sns{$ID})
		{
			$ID=int(rand($key_len));
		}
		$sns{$ID} = 1;
		push (@ID,$ID);
	}
	$key_id=$ID[0];
	$key01=$keys[$key_id];
	$key_id=$ID[1];
	$key02=$keys[$key_id];
	$key_id=$ID[2];
	$key03=$keys[$key_id];
	$key_id=$ID[3];
	$key04=$keys[$key_id];

	@mol_01=split(/t/,$Chime{$key01});
	@mol_02=split(/t/,$Chime{$key02});
	@mol_03=split(/t/,$Chime{$key03});
	@mol_04=split(/t/,$Chime{$key04});

	$point1=int(rand($mol_01[-1]));
	@mol_02_new= map {$_+ $point1} @mol_02;
	push (@mol_01,@mol_02_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
	
	$point2=int(rand($mol_01[-1]));
	@mol_03_new= map {$_+ $point2} @mol_03;
	push (@mol_01,@mol_03_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
  
	$point3=int(rand($mol_01[-1]));
	@mol_04_new= map {$_+ $point3} @mol_04;
	push (@mol_01,@mol_04_new);
	@pos_str=(sort {$a<=>$b} @mol_01);
	print "Chi:4\t$pos_str[-1]\n";
	push (@Chi_return,$pos_str[-1]);

	delete $Chime{$key01}; 
	delete $Chime{$key02}; 
	delete $Chime{$key03}; 
	delete $Chime{$key04}; 
	@keys=keys %Chime;
	$key_len=@keys;

	$return=join ("\t",@pos_str);
	push (@return,$return);
}

print "Chimera-last=$key_len\n";
#count length greater than or equal to 100 KB
my @data=grep (/$_ >= 100000/,@Chi_return);
my $true_chi=@data;
print "Number_Chimera_true=$true_chi\n";
if ($key_len>0) 
{
	foreach my $last (keys %Chime) 
	{
		@mol_01=split(/t/,$Chime{$last});
		$return=join ("\t",@mol_01);
		push (@return,$return);
	}
}

return @return;

}

1;



