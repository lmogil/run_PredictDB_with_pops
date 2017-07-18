# run_PredictDB_with_pops
####Lauren S. Mogil#####

###run PredictDB pipeline documentation###

###trial with example set####

###first run with alpha 0.5 then any other alphas needed

###download scrips edited for pop data adapted from https://github.com/hakyimlab/PredictDBPipeline/scripts

#format data according to guidlines in https://github.com/hakyimlab/PredictDBPipeline/wiki/Detailed_Description
#you will need 4 input files that are tab delimited and have headers 
#gene annotation file = gtf format ##this version uses gencode.v18 
#SNP annotation file with 7 columns 
#header = Chr	pos	var_id	ref_a1	alt_a2	snp_id	RSID_dbSNP137	num_alt
#var_id is formatted as chr_pos_refAllele_altAllele_build
#snp_id is formatted as 'snp'_chr_pos
#num_alt = 1
#example of what the file should look like 
	Chr     pos     var_id  ref_a1  alt_a2  snp_id  RSID_dbSNP137   num_alt
	1       719914  1_719914_C_G_b37        C       G       snp_1_719914    rs187772768     1

#genotype file with first column is the var_id and subsequent columns are samples/individuals 
#can use dosage data for genotypes 
#file example:
	
	rsid    26189   25542   24894   25825   25774   25101   26408   25779   26421   24811   26417   24968   24906   25213   25949   26185   25520   26339   26340   
	1_719914_C_G_b37        0.747    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    

#gene expression file with the first column as the ensemble_ID and the rest are individual
#example of expression file
	
	PROBE_ID        26189   25542   24894   25825   25774   25101   26408   25779   26421   24811   26417   24968   24906   25213   25949   26185   25520   26339   
	ENSG00000117620.8       0.045781135559082       0.0341434478759766      -0.226819038391113      -0.14402961730957       0.172945022583008       0.20235538482666

#covariate file (optional) used a PCs file to use to adjust gene expression

#1) run 05_make_files_predictDB.py to format input files correctly

#2) download predictDB from git clone https://github.com/lmogil/run_PredictDB_with_pops 
#run chmod +x joblogs/example/* script/*

#3) run make_dir_tree.sh this will make all the needed directories for the rest of these processes
#place input files in appropriate directories 

#4) make a directory for your population and copy scripts from example file

#5) modify model_parameters.py according to your study parameters and files created

#MESA AFA example the only thing that changes in model parameters is the alpha value

#Change these variables when adapting for different analyses. 

	List of identifiers for each database you'll make:
	STUDY_NAMES = ['AFA']
	File names for gene and snp annotation:
	GENE_ANNOTATION_FN = 'gencode.v18.annotation.gtf'
	SNP_ANNOTATION_FN = 'AFA_annot_chr1-22.txt'
	List of genotype/expression file names:
	GENOTYPE_FNS = ['AFA_snp_chr1-22.txt']
	EXPRESSION_FNS = ['AFA_MESA_Epi_GEX_data_sidno_Nk-20.txt']

	Model metadata/parameters. Keep all as strings:
	SNPSET = '1KG' #or HAPMAP
	ALPHA = '0.5' #change to 1 or 0.05 when changing this parameter
	N_K_FOLDS = '10'
	RSID_LABEL = 'RSID_dbSNP137'
	WINDOW = '1e6'


#6)run python preprocess_example.py
	#this runs a few scripts from the scripts file that separates data and places scripts in correct directories
	#parse_gtf.py 
	#geno_annot_to_RDS.R
	#split_snp_annot_by_chr.py
	#snp_annot_to_RDS.R
	#split_genotype_by_chr.py
	#Rscript expr_to_transposed_RDS.R 
	
##edited  expr_to_transposed_RDS.R to properly process and include covariates
	argv <- commandArgs(trailingOnly = TRUE)
	expressionfile <- argv[1]
	RDSout <- argv[2]
	Presence of covariate file suggests to correct for PEER factors, etc.
	covariatefile <- ifelse(length(argv) == 3, argv[3], NA)

	expression <- read.table(expressionfile, stringsAsFactors = FALSE,
    	header = TRUE, row.names = 1)
	Transpose expression.
	expression <- t(expression)

	if (!is.na(covariatefile)) {
  	Correct expression data for covariates.
  	covariate <- read.table(covariatefile, stringsAsFactors = FALSE,
    	header = TRUE, row.names = 1)
  	pcmat <- data.matrix(covariate[,-1:-3])
    	rownames(pcmat) <- covariate$FID

  	for (i in 1:length(colnames(expression))) {
    	fit <- lm(expression[,i] ~ pcmat)
    	expression[,i] <- fit$residuals
  	}
	}

	saveRDS(expression, RDSout)


##re-ran expr_to_transposed_RDS.R after preprocess_example.py with the following command
##make sure PCs are sorted according to expression file samples

	Rscript expr_to_transposed_RDS.R ../data/input/expression_phenotypes/AFA_MESA_Epi_GEX_data_sidno_Nk-20.txt \
    	../data/intermediate/expression_phenotypes/{expression_file}AFA_MESA_Epi_GEX_data_sidno_Nk-20.RDS \
    	../data/input/expression_phenotypes/afa_PCs_sorted.txt
##this will run a linear regression and account for PC variants 
    
    
#7) run python train_models.py
#see elastic net file for minor changes to the original
#changed output file paths and file name adjustment to the study names etc. 
#outputs all elastic net generated files into /data/intermediate/model_by_chr

	
#8)run python post_process.py
#this will generate concatenated beta, covariate, log, and results files of the population by chromosome
#this will also generate 2 db files one with all of the data and one filtered by FDR < 0.05


#9)transfer all results to separate directory so there's no file override 



