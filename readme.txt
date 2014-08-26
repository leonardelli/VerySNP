VerySNP is currently usable for Unix and MAC users only.
VerySNP is based on SVM classifier and discriminates between candidate variants using the informations contained in the VCF files obtained from GATK and SAMtools variant calling.
The package comprises two Python scripts:	- VerySNP_training for automatic training 
						- VerySNP_test for automatic model selection and prediction

VerySNP_training performs a parameter-search during a 10-fold cross validation on the provided training data. Afterwards VerySNP_test selects the best model (in terms of MCC - Matthews Correlation Coefficient), which is then used to predict unknown variants in test data. As output VerySNP_test assigns a class to each variant of the input file adding +1/-1 at the beginning of each variant line in the VCF format: if the variant belongs to the +1 class, it means the variant is likely true, while if it belongs to the -1 class, means that VerySNP predicted that variant as false.

VerySNP is available at https://github.com/leonardell/VerySNP

VerySNP is currently submitted for publication to Bioinformatics:
Lorena Leonardelli, Carmen Maria Livi, Patrice This, Charles Romieu, Claudio Moser and Alessandro Cestaro; "VerySNP: a SVM based tool to get accurate variant calling".


Table of Contents
=================

- Installation and Data Format
- Usage
- Examples
- Additional Information


Installation and Data Format
============================

To run VerySNP it is necessary to download the Libsvm package available at http://www.csie.ntu.edu.tw/~cjlin/libsvm [Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011]. It is enough to copy the binaries of 'svm-predict', 'svm-train' and 'svm-scale' into the same directory of VerySNP.

VerySNP takes as input VCF format files, which are currently limited to GATK (UnifiedGenotyper and VQSR) and SAMtools (mpileup) output files.
We provide another script (select_incommon_vcf_csv.py) to easily select the true and false sets from the variant calling output knowing their relative positions on the chromosomes/contigs.

Usage:
======

-Select the positive set:
python select_incommon_vcf_csv.py filename.vcf filename.csv > output.vcf
	filename.vcf : variant calling output
	filename.csv : your list of true variant positions organized on three columns (0: varaint_id, 1: contig_name, 2: variant_position) in csv format
	output.vcf : vcf file containing the true variants to train VerySNP 

Example:
python select_incommon_vcf_csv.py SAMtoolsVariantCalling.vcf my_TrueVariantList.csv > trueVariants_SAMtools.vcf

-Select the negative set:
python select_incommon_vcf_csv.py filename.vcf filename.csv > output.vcf
	filename.vcf : variant calling output
	filename.csv : your list of false variant positions organized on three columns (1: contig, 2: position) in csv format
	output.vcf : vcf file containing the false variants to train VerySNP

Example:
python select_incommon_vcf_csv.py SAMtoolsVariantCalling.vcf my_FalseVariantList.csv > falseVariant_SAMtools.vcf

-Training:
python VerySNP_training.py -p filename -n filename -k kernel -t type -tn training_name
	-p : vcf file containing the true variants dataset
	-n : vcf file containing the false variants dataset
	-k : kernel type (0: linear, 1: RBF)
	-t : file type (SAM, GATK)
	-tn: name of the folder holding the training output

Example:
python VerySNP_training.py -p SAMpos.vcf -n SAMpneg.vcf -k 0 -t SAM -tn training_sample1

-Prediction:
python VerySNP_test.py -f filename -k kernel -t type -tn training_name
	-f : vcf file for classification
	-k : kernel type (0: linear, 1: RBF) used for the training
	-t : file type (SAM, GATK)
	-tn: name of the folder holding the training output
    
Example:
python VerySNP_test.py -f testfile.vcf -k 1 -t GATK -tn training_sample1

Into the same folder holding the training output files will be created a folder named "output" holding the binary classified data (filename.vcf.results), which is the last VerySNP output.
Attention: To use the classification of VerySNP_test, a VerySNP_training has to be performed before.


Additional Information
======================
VerySNP is able to classify bi-allelic variants only; INDELs and tri-allelic polymorphisms are not taken into consideration neither for training nor for the classification. 

The training data should consist in experimentally validated variants (in case of the true dataset) and in potential variant positions resulted as non-variants by the biological validation (in case of the false dataset). The training data should contain a similar amount of true and false variants (i.e. 400 true variants and 400 false variants) to avoid an unbalanced training of the classification model. A further extension of VerySNP will include an automatic re-balancing if the training datasets should differ in their amount. Larger the training datasets, better VerySNP learns how to recognize true and false variants.

When the training step (VerySNP_training) is completed, 11 files (taking around 300K of memory) are created into the <training_name> folder and are required to accomplish the second step (VerySNP_test), where VerySNP applies the model to the variant calling output. You can delete the 11 files when the test has been successfully accomplished.

Into VerySNP package have been added some examples of training sets named as OrganismName_Tool_True/FalseVariants.vcf, like yeast_gatk_true.vcf and yeast_gatk_false.vcf include a subset of the true and false variants we used to train VerySNP with in order to classify the GATK variant calling on the yeast sample. 

Example: 
python VerySNP_training.py -p yeast_gatk_true.vcf -n yeast_gatk_false.vcf -k 0 -t GATK -tn yeastGATK

Other files included into VeryNSP package are GATK and SAMtools whole variant calling of  the three samples described in the paper (like yeast_gatk_variant_call.vcf and yeast_sam_variant_call.vcf). The whole variant calling VCF files include the unknown variants we then classified with VerySNP_test.py.

Example:
python VerySNP_test.py -f yeast_gatk_variant_call.vcf -k 0 -t GATK -tn yeastGATK
