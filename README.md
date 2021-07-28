# BMSIM: BioNano Molecule SIMulator
================

BMSIM is a simulation tool for BioNano molecule data of the BioNano optical mapping platform

---------------


SUMMARY
---------------
BioNano Molecule SIMulator (BMSIM) explicitly incorporated BioNano data models (BioNano molecule length distribution, FN and FP signals, DNA molecules stretching variations, variation in optical resolution, and fragile sites) and the methods to generate chimeric molecules and assign SNR scores for simulated BioNano molecules. We simulated noisy maps from ‘perturbed’ versions of the reference map. Using genomic sequences (.fasta file) as input, BMSIM simulated noisy maps with five main steps:
	I) generate BioNano molecules with random fragmentation and fragile site bias model; 
	II) abel nicking sites for BioNano molecules by in silico restriction digestion.Our program supported all available nicking enzymes currently used in BioNano systerm 	(i.e., Nt.BspQI, Nb.BbvCI, Nb.Bsml and Nb.BsrDI), as well as any artificial nicking sequences that users chose to define; 
	III) incorporate data models for FN sites, FP sites, stretching variations, optical resolution, and chimerism for BioNano molecules; 
	IV) assign SNR and intensity scores for labelling sites; 
	V) iterate for targeted coverage depth. 
	The output of BMSIM is a BNX format text file (.BNX, see example EcoliSim.bnx file) which contains molecule map length, label positions, and label signal score, ect.

 
DEPENDENCIES
---------------
BMSIM requires Perl >= 5.10.0;

Perl module Math::Random. 
This can be installed using CPAN http://search.cpan.org/~grommel/Math-Random-0.70/Random.pm;

Perl module Math::Random::MT. 
This can be installed using CPAN http://search.cpan.org/~fangly/Math-Random-MT-1.17/MT.pm;

    
USAGE
---------------    
perl BMSIM.pl [options]

Example：                                
 
 perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -fragile ./example/MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 -p EcoliSim -bnx out.bnx >log.txt

Documentation options:
 	
	-help	brief help message

Required options:
	
	-cov	Coverage, e.g. 10, represent simulating 10X coverage depth for the simulated genome
	-ca	reference genome fasta, e.g. example/MG1655.fa
	-fragile	fragile site file, e.g. example/MG1655.txt
	-bnx	bnx output, e.g. example/EcoliSim.bnx
	
Optional options:

	-lm	length_mean, #possion
	
	-FNp	FN_probality, Binomil
	-FPi	FP_intensity, possion~lamda
	
	-str	stretch_normal, stretch~N(ave,std),comma-separated lists of values
	
	-Rt	Resolution type
	-RtI	Res_typeI, under1000,under1500
	-RtII	Res_typeII, resolution~N(ave,std)
	
	-snr	SNR, SNR~N(ave,std)
	-Ints	Intensity, Intensity~N(ave,std)
	
	-Chi_pert	Chi_percent
	-Chi_Bi		Chi_Bi
	-Chi_Tr		Chi_Tri
	-Chi_Qua	Chi_Qua
	
	-e	enzyme nicking pattern, e.g. GCTCTTC,GAAGAGC for NtBspQI
	-np1	nick_position1, nicking position of the enzyme pattern in 5' to 3',
             e.g. NtBspQI: 5' ...GCTCTTCN^... 3' the np1 is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the np1 is -1.
                           3' ...CAGAGNN... 5'                              3' ...CGTTAC^NN... 5'
			 
	-np2	nick_position2
	-f	FragileTypeI
	-p	project name for all assemblies
	




