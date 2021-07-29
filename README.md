# BMSIM: BioNano Molecule SIMulator
================

BMSIM is a simulation tool for BioNano molecule data of the BioNano optical mapping platform.

---------------


SUMMARY
---------------
	BioNano Molecule SIMulator (BMSIM) explicitly incorporated BioNano data models (BioNano molecule length distribution, FN and FP signals, DNA molecules stretching variations, variation in optical resolution, and fragile sites) and the methods to generate chimeric molecules and assign SNR scores for simulated BioNano molecules. We simulated noisy maps from ‘perturbed’ versions of the reference map. Using genomic sequences (.fasta file) as input, BMSIM simulated noisy maps with five main steps:
	I) generate BioNano molecules with random fragmentation and fragile site bias model; 
	II) add nicking sites for BioNano molecules by in silico restriction digestion. Our program supported all available nicking enzymes currently used in BioNano systerm 	(i.e., Nt.BspQI, Nb.BbvCI, Nb.Bsml and Nb.BsrDI), as well as any artificial nicking sequences that users chose to define; 
	III) incorporate data models for FN sites, FP sites, stretching variations, optical resolution, and chimerism for BioNano molecules; 
	IV) assign SNR and intensity scores for labelling sites; 
	V) iterate for targeted coverage depth. 
	The output of BMSIM is a BNX format text file (.bnx, see example file "EcoliSim.bnx" ) which contains molecule map length, label positions, and label signal score, ect.

 
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
 
 	Example 1: Use single nicking enzyme NtBspQI, NtBspQI: 5' ...GCTCTTCN^... 3', the nick position is 7
	                                                       3' ...CGAGAAGN... 5'                            
 	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -Fi ./example/MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 -p EcoliSim1 -bnx out1.bnx
 
 
	Example 2: Use dual nicking enzymes: NtBspQI and Nb.BsrDI
            NtBspQI: 5' ...GCTCTTCN^... 3' the nick position is 7;    Nb.BsrDI: 5' ...GCAATGNN... 3'  the nick position is -1
                     3' ...CGAGAAGN... 5'                                       3' ...CGTTAC^NN... 5'
	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -Fi ./example/MG1655.txt -e GCTCTTC,GAAGAGC,GCAATG,CATTGC -np1 7 -np2 -1 -p EcoliSim2 -bnx out2.bnx
 
 
 	Example 3: Change parameters (mean and standard deviation) for the stretch model:
	perl BMSIM.pl -cov 10 -ca ./example/MG1655.fa -Fi ./example/MG1655.txt -e GCTCTTC,GAAGAGC -np1 7 -str 0.95,0.05 -p EcoliSim3 -bnx out3.bnx


Documentation options:

	-h,--help      help message

Required options:

	-cov,--Coverage int          Coverage depth for the simulated genome, e.g. 10, represent simulating 10X coverage depth for the simulated genome
	-ca,--chr_fa string          file with genome sequence(.fa), e.g. example/MG1655.fa
	-Fi,--Fragile_file string    file with fragile site, e.g. example/MG1655.txt
	-bnx,--bnx_file string       output file(.bnx), e.g. example/EcoliSim.bnx

Optional options:
  
	-lm,--length_mean int                  mean length (kb) of BioNano molecule, produce BioNano molecular length according to the exponential distribution (Default:70)
  
	-FNp,--FN_probability float            parameter for the false negative site model, i.e. Binomil distibution, FNp is the probability of getting a sucess in one trial (Default:0.85)
  
	-FPi,--FP_intensity float              the initial interval of false positive site (kb) (Default:100)
  
	-str,--stretch_normal string           parameters (ave,std) for the stretch model,i.e. Normal distribution, Stretch~N(ave,std), ave is mean, std is standard deviation, comma-separated lists of values (Default:0.98,0.08)
  
	-Rt,--Res_type int                     resolution type, "1" for "Res_typeI", "2" for "Res_typeII" (Default:2)
	-RtI,--Res_typeI string                resolution model type I:p1,p2, p1 is the proportion of sites within 1kb distance that are retained, p2 is the proportion of sites within 1.5kb distance that are retained, comma-separated lists of the values p1 and p2 (Default:0,0.9)
	-RtII,--Res_typeII string              parameters (ave,std) for the resolution model type II,i.e. Normal distribution, Res_typeII~N(ave,std), ave is mean, std is standard deviation, comma-separated lists of values (Default:1.2,0.9)
  
	-snr,--SNR string                      parameters (ave,std) for signal noise ratio (SNR) model of BioNano molecule,i.e. Normal distribution, SNR~N(ave,std), ave is mean, std is standard deviation, comma-separated lists of values (Default:3,0.66)
  
	-Ints,--Intensity string               parameters (ave,std) for Intensity model of BioNano molecule, i.e. Normal distribution, Intensity~N(ave,std), ave is mean, std is standard deviation, comma-separated lists of values (Default:1,0.2)
  
	-ChimP,--Chimeric_proportion float     proportion of the total chimeric molecule in raw BioNano molecule (Default:0.1)
	-Bi,--Bimera_proportion float          proportion of Bimera (Default:0.85)
	-Tri,--Trimera_proportion float        proportion of Trimera (Default:0.13)
 	-Quad,--Quadramera_proportion float    proportion of Quadramera (Default:0.02)
  
 	-e,--enzymes_pattern string            nicking pattern of the enzyme, e.g. GCTCTTC,GAAGAGC for NtBspQI; GCAATG,CATTGC for Nb.BsrDI; GCTCTTC,GAAGAGC,GCAATG,CATTGC for dual nicking with NtBspQI and Nb.BsrDI (Default:GCTCTTC,GAAGAGC)
	-np1,--nick_position1 int              nicking position of the enzyme pattern from 5' to 3' for enzyme 1,
                                               e.g. NtBspQI:5' ...GCTCTTCN^... 3' the np1 is 7; Nb.BsrDI:5' ...GCAATGNN... 3'  the np1 is -1.
                                                            3' ...CGAGAAGN... 5'                         3' ...CGTTAC^NN... 5'
	-np2,--nick_position2 int              nicking position of the enzyme pattern from 5' to 3' for enzyme 2,
                                               e.g. NtBspQI:5' ...GCTCTTCN^... 3' the np2 is 7; Nb.BsrDI: 5' ...GCAATGNN... 3'  the np2 is -1.
                                                            3' ...CGAGAAGN... 5'                          3' ...CGTTAC^NN... 5'
  
	-p,--proj string                       project name for the simulation, e.g. EcoliSim
	




