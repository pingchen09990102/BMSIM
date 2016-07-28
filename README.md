# BMSIM
BioNano data simulator
BMSIM requires Perl >= 5.10.0
================

scripts to simulate BioNano molecule data

---------------

SUMMARY

BioNano molecules were simulated from the in silico reference genome according to the BioNano experiments process. We simulated noisy maps from ‘perturbed’ versions of the reference map. As shown in Workflow diagram, BMSIM system included five main steps:

Workflow diagram

 1) Random fragmentation (including a small portion of specific fragmentation in fragile sites loci). Provided the genomic sequence (fasta file) of an organism as input, the initial step is random fragmentation. Since the dsDNA is more likely to break at fragile site block during labeling, we estimated the breakage probability of each candidate fragile sites and introduced markers to corresponding genomic region before random fragmentation. After fragile sites marking, we fragment the reference genome sequence with homogeneous poisson process. Subsequently, molecules containing fragile sites markers are broken up in the marker loci.
 
 2) In silico enzyme digestion. The second step of BioNano experiment is labeling nicking sites on DNA. We simulate this step by In silico digestion of the DNA sequence obtained above with user-defined nick enzyme pattern (one or more nick enzyme pattern are supported). The in silico digestion process computationally mimics how each nicking enzyme would cleave the segment of DNA defined by the molecules, returning “in silico BioNano molecules” that can be used in following steps.
 
 3) Add error models. Errors or bias are introduced by labeling and imaging steps. For missing label due to failed digestion, we treat the digestion of each nicking site as a Bernoulli event with probability of success p. False labels are modeled as Poisson Process with rate  . For sizing bias, we calculated the reported fragment size between each pair of labels as X=RZ, where Z is the bp size between the paired labels, and R is a stretch scale factor, which modeled as  . Chimeras are generated by appending consecutive parent molecules at the junction, where the parent molecules and junctions are randomly selected in frame. The proportion of every chimera type, i.e.‘Bimera’, ‘Trimera’, or ‘Quadramera’, can be defined by users. To model resolution limitation, since the BioNano system resolution is ~1kb, and the likelihood of resolving two sites are different (the far the two neighbor sites apart from, the more possible they can be recognized), we implement resolution model with Gaussian distribution  . The likelihood of resolving in general is modeled by the cumulative Gaussian.
 
 4) Add label signal score SNR Label signal scores are generated based on two log-normal distributions. BMSIM assigns the one (with relative lower mean and std) to introduced errors (i.e. FPs), and the other (with higher mean and std) to all other labels. A signal score at each label position of a simulated molecule is randomly chosen from a frequency table of signal scores.
 
 5) Coverage depth. To obtain enough data, we iterate step I to step IV until reach the predefined coverage depth.
 
 The output of BMSIM is a BNX format (developed by BioNano) text file which contains molecule map length, label positions, and label signal score, ect. (see example BNX file) .

 
 DEPENDENCIES


    Perl module Math::Random. This can be installed using CPAN ;
    Perl module Math::Random::MT. This can be installed using CPAN http://search.cpan.org/~fangly/Math-Random-MT-1.17/MT.pm;

    
USAGE
    
perl BMSIM.pl [options]

 Documentation options:
 	-help	brief help message
	-man	full documentation
 
Required options:
	-cov	Coverage
	-ca	reference genome fasta
	-bnx	bnx output
	-fragile	fragile site file
 
Optional options:
	-lm	length_mean,#possion
	
	-FNp	FN_probality,Binomil~
	-FPi	FP_intensity,possion~lamda
	
	-str	stretch_normal,stretch~N(ave,std),comma-separated lists of values
	
	-Rt	Resolution type
	-RtI	Res_typeI, under1000,under1500
	-RtII	Res_typeII, resolution~N(ave,std)
	
	-snr	SNR, SNR~N(ave,std)
	-Ints	Intensity, Intensity~N(ave,std)
	
	-Chi_pert	Chi_percent
	-Chi_Bi	Chi_Bi
	-Chi_Tr	Chi_Tri
	-Chi_Qua	Chi_Qua
	
	-e	enzyme nicking pattern
	-np1	nick_position1, nicking position of the enzyme pattern in 5' to 3',
             e.g. NtBspQI: 5' ...GCTCTTCN^... 3' the np1 is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the np1 is -1.
                           3' ...CAGAGNN... 5'                              3' ...CGTTAC^NN... 5'
			 
	-np2	nick_position2
	-f	FragileTypeI
	-p	project name for all assemblies
	

    #example：                                
    perl BMSIM.pl -cov 10 -p PLsimCov350 -ca MG1655.fa -bnx PLsimCov350.bnx -fragile MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 >log.txt


