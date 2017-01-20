# AmpSeq-SSR
An accurate and efficient method for large-scale SSR genotyping and applications
Copyright Holder:Lun Li, Zhiwei Fang and Weixiong Zhang  
email: weixiong.zhang@wustl.edu.

* Program Requirement
	Running AmpSeq-SSR requires a GNU-like environment. It should be possible to
        run AmpSeq-SSR on a Linux or Mac OS.

	AmpSeq-SSR requires following external tools:
	1. bowtie2
	2. perl
	3. blast tool suite	
* Running 
	cd ./data
	sh AmpSeqSSR.sh sample1.fa ampliconRef.fa primer.fa ../script sample1

* AmpSeqSSR.sh takes five parmeters:
	1. sample1.fa: sequencing reads of amplicons for sample1 in fasta format:              
	2. ampliconRef.fa: corresponding reference sequences of amplicons:              
	3. primer.fa:  primers of amplicons                             
	4. ../script:   the path of script directory
	5. sample1:    sample name
				
* Format of sequencing reads fasta file:
	> &gt;seq343_x100	GCTGCCATCCCTGTACATACACAGACACTATACTTGTTTTACAGTTCTTACACTTGGGCAGCATCATGTCAATGGGTGTTTGGCCAAGCCCAGGACAGCCCCAAAATTTTTCAGCAGAATTCTGGGACTAGTTGTCATTTTTTCTTTTTTCTGTGGTAAATTTAACTTAAGCTTGTTTAGTTATAGTTATGCAGCATCC
	* each item stands for a group of reads with exactly the same sequence, and the group is named like "seqN_xM",where N indicates the unique ID and M represents the number of reads in this group.

* Format of primer fasta file:
	> &gt;AMPL1_F  
	GGGTTCTCGGTGAAGATGGC  
	> &gt;AMPL1_R  
	AGCGCCTCGATGAGCTTG
	* Primer file contains primer sequences of both strand
	
* Output files

	1. sample1.allele_profile: contains allele informations of all simple SSRs within amplicons. And the allele profile file is in following format:  
	ampliconID	ssr_position	motif	stutter_ratio	major_allele:read_count	second_allele:read_count	third_allele:read_count ...  
	e.g.  
	>AMPL1     1_1     cgc     0.043859649122807       CGC8:114        CGC7:5  CGC6:5  CGC4:4  CGC5:3  CGC4_CGC4:2     CGC3:1
	
	2. sample1.stat.addLen.xls: amplicon informations  
	ampliconID      total_reads     average_length_of_major_SSR number_of_invalid_reads        major_allele        major_allele_read_count        second_allele      seconde_allele_read_count       stutter read_count_of_other_alleles  
	e.g.  
	>AMPL1563820     267     199     0       CTAT5,  243     CTAT4,  12      0.0493827160493827  12
