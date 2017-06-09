##PULSAR 1.0 (Phasing Utilizing Lineage Specific Alleles / Rare variants)
This program performs phasing of Whole Genome Sequencing genotypes in pedigrees utilizing lineage specific alleles to establish the IBD pattern within the pedigrees.

###PULSAR_1.0.pl
* This program performs phasing of Whole Genome Sequencing genotypes in pedigrees utilizing lineage specific alleles to establish the IBD pattern within the pedigrees. Input genotypes are in VCF format. Additional files are required containing the pedigree and allele frequencies. See example folder.
* Usage is as follows: >perl PULSAR_1.0.pl $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3] $ARGV[4]
* where $ARGV[0] is the pedigree file in PULSAR format (see example), $ARGV[1] is the genotype file in VCF format, $ARGV[2] is an integer above 0 that is used as the number of adjacent observations that are required to identify a recombination event (we suggest 5), $ARGV[3] is a file containing the MAF file in PULSAR format (this should be sorted and include the same variants as the vcf), and $ARGV[4] is MAF threshold to use as a filter (for example, 0.05 means MAF < 5%) 

###IBD_sharing_estimator.pl
* This program estimates the percentage of genome in which each sequenced individual will share at least one of the two homologous chromosomes they carry IBD with another sequenced individual in the pedigree.
*Usage is as follows: >perl IBD_sharing_estimator.pl $ARGV[0] $ARGV[1]
*where $ARGV[0] is the pedigree file in PULSAR format and $ARGV[1] is the number of MC iterations to perform

###Example/
* This folder provides an example for a small pedigree. Either copy the programs to this folder or call them from the directory the programs are in.
* To run the example: >perl PULSAR_1.0.pl Seven_sib_pedigree.csv Genotypes.vcf 5 Allele_frequency.csv 0.05
* To run the IBD sharing program: >perl IBD_sharing_estimator.pl Seven_sib_pedigree.csv
