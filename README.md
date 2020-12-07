# VFDetector

An automatic pipeline for the identification, quantification and visualization of virulence factors in metagenomic sequencing data

Dependencies:
Python >=3.7
BOWTIE2
DIAMOND == v0.8.32.94

Install:
Install required python libraries: pip install pandas numpy  matplotlib 

git clone https://github.com/CicyLab/VFDetector.git

Usage:
python VFdetector.py PREDICT --type prot -i simple.fq.gz -o prefix -d DATABASE/ --threads 4

Output:
(1). prefix.VFG.XLS
GID	GeneLength	UniqCount	MultiCount	Abundance
VFG000017	104.000000	57.000000	43.500000	551.252797
VFG000016	102.000000	23.000000	65.000000	492.153507
VFG000014	133.000000	105.000000	0.000000	450.355781
VFG000013	152.000000	105.000000	0.000000	394.061308
VFG000004	204.000000	2.000000	108.333333	308.528051

(2). prefix.VF.XLS
Virulence Factors	Value	VF Name	Subtype	Associated Bacteria	Mechanism	Function
VF0105  2.950607        Lpf     Adherence       Salmonella enterica (serovar typhimurium)       NA      Mediate attachment to the Peyer's patches
VF0102  0.409807        Type 1 fimbriae Adherence       Salmonella enterica (serovar typhimurium)       NA      The adhesin FimH mediates T3SS1-independent uptake in murineDCs.
AI058   0.551186        NA      Adherence|Invasion      NA      NA      NA

(3). alignment files(.sam)
@HD	VN:1.5	SO:query
@PG	PN:DIAMOND
@mm	BlastX
@CO	BlastX-like alignments
@CO	Reporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate
VFG000001|VF0025|Adherence;Secretionsystem_0_0/1	0	VFG000001	2686	255	27M	*	0	0	AGTLDGKMQNLEI*GGSVDAAHTDLSV	*	AS:i:52	NM:i:2	ZL:i:3590	ZR:i:125	ZE:f:7.2e-09	ZI:i:92	ZF:i:3	ZS:i:3	MD:Z:4E8E13
VFG000001|VF0025|Adherence;Secretionsystem_1_0/1	0	VFG000001	1914	255	41M	*	0	0	EYFKTPLPVSLTALDNRAGLSPATWNFQSTYELLDYLLDQN AS:i:89	NM:i:0	ZL:i:3590	ZR:i:220	ZE:f:7.0e-20	ZI:i:100	ZF:i:3	ZS:i:3	MD:Z:41
