#!/bin/bash 

# This script is the final backend for the fedplink client that runs online with file transfers. It utilizes the options from FILE_IO_UTILS.sh to perform file transfers/waiting etc.
# Any interaction with the shared folder must be done via FILE_IO_UTILS.sh.

# Make sure the number of scored variant blocks do not surpass the number of variants that exist in the data.

if [[ $# -lt 2 ]]
then
	echo "USAGE: $0 [options] [arguments]
	Genotype importing/Variant filtering Steps:
		-separate_VCF_per_chrom [Data config file] [VCF file path] [chromosome id's list path] [Output directory]
		-import_VCF [Data config file] [VCF file path] [chromosome identifier to use] [Imported genotypes directory] : -convert_VCF_2_matbed and save in the imported directory.
		-extract_subset_variants [Data config file] [Imported genotypes directory] [Variants BED file] [Output directory]
		-extract_subset_subjects [Data config file] [Imported genotypes directory] [New subject IDs file] [Output directory]
		-save_GMMAT_genotypes [Data config file] [Imported/filtered genotype directory] [Output GMMAT file]
		-impute_missing_GMMAT_formatted_geno [Data config file] [GMMAT file] [Imputed GMMAT file]
	Data Simulation/Conversion Steps:
		-split_GMMAT_dir [Data config file] [Pooled data GMMAT directory] [Per site output base directory]
		-convert_GMMAT_dir [Data config file] [Pooled data GMMAT directory] [Output directory]
	KeyMaker Steps:
		-generate_keys [Data config file] [Output directory]
	Encryption/Decryption (openssl):
		-symmetric_encrypt_partdec_data [Data config file] [partdec file] [openssl hash file] [Output file]
		-symmetric_decrypt_partdec_data [Data config file] [encrypted partdec file] [openssl hash file] [Output file]
	Null-Model Fitting Steps:
		-validate_input_data [Data config file] [Local data directory]
		-client_calculate_save_XtWX_XtWz [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-add_noise_to_XtWX [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-pool_site_specific_all_site_noise_XtWX [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-collaborative_decrypt_pooled_noisy_XtWX [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-pool_partially_decrypted_pooled_noisy_XtWx_remove_noise [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-pool_XtWX_XtWz_update_beta [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-client_collaborative_decrypt_beta [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-client_pool_partially_decrypted_beta [Data config file] [Epoch index] [# Epochs] [Client index (starts at 0)] [Local data directory]
		-check_convergence_per_updated_beta [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
		-get_last_win_start_i [Data config file]
	P-value Statistic Calculation:
		-client_calculate_save_pvalue_stats [Data config file] [Client index (starts at 0)] [Local data directory]
		-client_pool_pvalue_stats [Data config file] [Client index (starts at 0)] [Local data directory]
	P-value Assignment:
		-pool_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Local data directory]
		-client_collaborative_decrypt_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Local data directory]
		-client_pool_partially_decrypted_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Local data directory]
		-pool_final_plaintext_p_values [Data config file] [Client index (starts at 0)] [Local data directory]
	MetaGMMAT Steps:
		-validate_meta_analysis_input_data [Data config file] [Local data directory]
		-meta_encrypt_save_ST_stats [Data config file] [Client index (starts at 0)] [Output directory]: Encrypt and upload S and T statistics.
		-meta_pool_per_site_ST_stats_add_noise [Data config file] [Client index (starts at 0)] [Output directory] : Pool Gt(Y-Y0) and T from all sites, calculate S for each variant; add single noise to each entry in S and T.
		-meta_pool_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Output directory] : Pool the noise added S and T stats from all sites.
		-meta_collaborative_decrypt_pooled_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Output directory]: Decrypt meta ST stats: Each site adds one set of multiplicative noise.
		-meta_pool_partially_decrypted_pooled_noisy_ST_stats [Data config file] [Client index (starts at 0)] [Local data directory] : Pool results, calculative S/T and p-value for each variant."

	exit
fi

##########################################################################################################################################################################
# DO NOT WRITE ANYTHING ON THE COMMAND LINE OUTSIDE ANY CMD OPTION, SOME OF THE OPTIONS WRITE VALUES THAT ARE READ BY OTHER OPTIONS.
# DO NOT WRITE ANYTHING ON THE COMMAND LINE OUTSIDE ANY CMD OPTION, SOME OF THE OPTIONS WRITE VALUES THAT ARE READ BY OTHER OPTIONS.
# DO NOT WRITE ANYTHING ON THE COMMAND LINE OUTSIDE ANY CMD OPTION, SOME OF THE OPTIONS WRITE VALUES THAT ARE READ BY OTHER OPTIONS.
##########################################################################################################################################################################

cmd_option=$1
data_config_file=$2

# We need the data configuration file in any run of this script.
if [[ ! -f ${data_config_file} ]]
then
	echo "Could not find data config file @ \"${data_config_file}\""
	exit
fi

source ${data_config_file}

if [[ ${cmd_option} == "-separate_VCF_per_chrom" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [VCF file path] [chromosome id's list path] [Output directory]"
		exit
	fi

	vcf_file=$3
	chrom_ids_file=$4
	op_dir=$5

	if [[ ! -f ${vcf_file} ]]
	then
		echo "Could not find VCF file @ ${vcf_file}"
		exit 1
	fi

	if [[ ! -f ${chrom_ids_file} ]]
	then
		echo "Could not find chromosome identifiers file @ ${chrom_ids_file}"
		exit 1
	fi

	if [[ ! -d ${op_dir} ]]
	then
		echo "Could not find output directory @ ${op_dir}"
		exit 1
	fi

	${COLLAGENE_SECURE_EXEC} -separate_VCF_2_chroms ${vcf_file} ${chrom_ids_file} ${op_dir}

	exit
fi

if [[ ${cmd_option} == "-import_VCF" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [VCF file path] [chromosome id for VCF file] [Output directory]"
		exit
	fi

	vcf_file=$3
	chrom_id_2_use=$4
	op_dir=$5

	if [[ ! -f ${vcf_file} ]]
	then
		echo "Could not find ${vcf_file}";
		exit 1
	fi

	if [[ ! -d ${op_dir} ]]
	then
		echo "Could not find output directory @ ${op_dir}"
		exit 1
	fi

	bin_seq_dir="none"
	match_ref_allele="0"
	match_reg_names="1"
	haplo_spec_enc="0"
	op_file=${op_dir}/${chrom_id_2_use}.matbed.gz

	# Extract the variants on this vcf. 
	vcf_filename=$(basename -- "$vcf_file")
	extension="${vcf_filename##*.}"

	cat_vcf_cmd="cat ${vcf_file}"	
	if [[ $extension == "gz" ]]
	then
		echo "VCF: gzipped.."
		cat_vcf_cmd="gzip -cd ${vcf_file}"
	elif [[ $extension == "gzip" ]]
	then
		echo "VCF: gzipped.."
		cat_vcf_cmd="gzip -cd ${vcf_file}"
	elif [[ $extension == "bzip2" ]]
	then
		echo "VCF: bzipped.."
		cat_vcf_cmd="bzip2 -cd ${vcf_file}"
	elif [[ $extension == "bz2" ]]
	then
		echo "VCF: bzipped.."
		cat_vcf_cmd="bzip2 -cd ${vcf_file}"
	else
		echo "VCF: plain.."
		cat_vcf_cmd="cat ${vcf_file}"		
	fi

	# Get the sampel identifiers.
	${cat_vcf_cmd} | head -n 1000 | awk {'if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i}}'} > ${op_dir}/GENO_SAMPLE_IDS.list

	sample_size=`wc -l ${op_dir}/GENO_SAMPLE_IDS.list | awk {'print $1'}`
	echo "Found ${sample_size} subjects in ${vcf_file}."

	# Get the variants.
	${cat_vcf_cmd} | awk 'BEGIN{FS="\t"}{if(substr($1,1,1)!="#"){l_var=length($4);start_pos=$2;ref_all=$4;alt_all=$5;var_id=$3"_"ref_all"_"alt_all;print $1"\t"start_pos-1"\t"start_pos+l_var-1"\t"var_id"\t.\t+"}}' | gzip > ${op_dir}/${chrom_id_2_use}_vars.bed.gz

	n_vars=`gzip -cd ${op_dir}/${chrom_id_2_use}_vars.bed.gz | wc -l | awk {'print $1'}`
	echo "Found ${n_vars} variants."

	# Import the variants.
	VCF_vars_bed_file=${op_dir}/${chrom_id_2_use}_vars.bed.gz
	${COLLAGENE_SECURE_EXEC} -convert_VCF_2_matbed ${vcf_file} ${op_dir}/GENO_SAMPLE_IDS.list ${VCF_vars_bed_file} ${chrom_id_2_use} ${bin_seq_dir} ${match_ref_allele} ${match_reg_names} ${haplo_spec_enc} ${op_file}

	exit
fi

if [[ ${cmd_option} == "-extract_subset_subjects" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [matbed genotypes dir.] [New subject IDs file] [Output directory]"
		exit
	fi

	imported_geno_dir=$3
	sub_subject_ids_file=$4
	op_geno_dir=$5

	if [[ ! -d ${imported_geno_dir} ]]
	then
		echo "Could not find imported genotype directory @ ${imported_geno_dir}"
		exit 1
	fi

	if [[ ! -d ${op_geno_dir} ]]
	then
		echo "Could not find output genotype directory @ %{op_geno_dir}"
		exit 1
	fi

	if [[ ! -f ${sub_subject_ids_file} ]]
	then
		echo "Could not find the subject IDs list file @ ${sub_subject_ids_file}"
		exit 1
	fi

	sample_ids_list_file=${imported_geno_dir}/GENO_SAMPLE_IDS.list
	if [[ ! -f ${sample_ids_list_file} ]]
	then
		echo "Could not find the sample id's list file @ ${sample_ids_list_file}"
		exit 1
	fi

	# Copy the sample identifiers first.
	cp ${sub_subject_ids_file} ${op_geno_dir}/GENO_SAMPLE_IDS.list

	# Process each chromosome.
	for cur_file in ${imported_geno_dir}/*.matbed.gz
	do
		matbed_file=${cur_file}
		file_name=`echo ${matbed_file} | xargs -Ifile basename file`
		op_file=${op_geno_dir}/${file_name}

		echo "Extracting subjects in ${sub_subject_ids_file} from ${matbed_file} and saving to ${op_file}"
		${COLLAGENE_SECURE_EXEC} -extract_subset_subjects ${matbed_file} ${sample_ids_list_file} ${sub_subject_ids_file} ${op_file}
	done

	exit 0
fi

if [[ ${cmd_option} == "-extract_subset_variants" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [matbed genotypes dir.] [variant regions BED path] [Output directory]"
		exit
	fi

	imported_geno_dir=$3
	variants_bed_file=$4
	op_geno_dir=$5

	if [[ ! -d ${imported_geno_dir} ]]
	then
		echo "Could not find imported genotype directory @ ${imported_geno_dir}"
		exit 1
	fi

	if [[ ! -d ${op_geno_dir} ]]
	then
		echo "Could not find output genotype directory @ %{op_geno_dir}"
		exit 1
	fi

	if [[ ! -f ${variants_bed_file} ]]
	then
		echo "Could not find the variants BED file @ ${variants_bed_file}"
		exit 1
	fi

	sample_ids_list_file=${imported_geno_dir}/GENO_SAMPLE_IDS.list
	if [[ ! -f ${sample_ids_list_file} ]]
	then
		echo "Could not find the sample id's list file @ ${sample_ids_list_file}"
		exit 1
	fi

	# Copy the sample identifiers first.
	cp ${sample_ids_list_file} ${op_geno_dir}

	# Process each chromosome.
	for cur_file in ${imported_geno_dir}/*.matbed.gz
	do
		matbed_file=${cur_file}
		file_name=`echo ${matbed_file} | xargs -Ifile basename file`
		op_file=${op_geno_dir}/${file_name}

		echo "Extracting variants in ${variants_bed_file} from ${matbed_file} and saving to ${op_file}"
		${COLLAGENE_SECURE_EXEC} -extract_variant_regions ${matbed_file} ${sample_ids_list_file} ${variants_bed_file} ${op_file}
	done

	exit 0
fi

if [[ ${cmd_option} == "-impute_missing_GMMAT_formatted_geno" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE $0 $1 [Data config file] [GMMAT file] [Imputed GMMAT file]"
		exit
	fi

	gmmat_geno_file=$3
	imp_gmmat_geno_file=$4

	${COLLAGENE_SECURE_EXEC} -impute_missing_GMMAT_formatted_geno ${gmmat_geno_file} ${imp_gmmat_geno_file}
fi
 

if [[ ${cmd_option} == "-save_GMMAT_genotypes" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE $0 $1 [Data config file] [Imported genotypes directory] [Output GMMAT file]"
		exit
	fi

	imported_geno_dir=$3
	op_gmmat_geno_file=$4

	if [[ ! -d ${imported_geno_dir} ]]
	then
		echo "Could not find imported genotypes directory @ ${imported_geno_dir}"
		exit 1
	fi

	geno_sample_ids_list_file=${imported_geno_dir}/GENO_SAMPLE_IDS.list
	if [[ ! -f ${geno_sample_ids_list_file} ]]
	then
		echo "Could not find genotype sample identifiers file @ ${geno_sample_ids_list_file}"
		exit 1
	fi

	matbed_list_file=temp_matbeds.list
	ls ${imported_geno_dir}/*.matbed.gz > ${matbed_list_file}

	${COLLAGENE_SECURE_EXEC} -write_GMMAT_formatted_geno ${matbed_list_file} ${geno_sample_ids_list_file} ${op_gmmat_geno_file}

	exit
fi

if [[ ${cmd_option} == "-get_last_win_start_i" ]]
then
	if [[ $# != 2 ]]
	then
		echo "USAGE: $0 $1 [Data config file]"
		exit
	fi

	source ${data_config_file}
	requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'req_win_end_i=$2 - ($2 % $1);win_start_i=(req_win_end_i-$1);if(req_win_end_i<$1){win_start_i=0};print win_start_i'}`

	echo ${requested_last_block_start_i}
	exit
fi

if [[ ${cmd_option} == "-convert_GMMAT_dir" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Pooled data GMMAT directory] [Output directory]"
		exit
	fi

	source ${data_config_file}

	GMMAT_dir=$3
	INPUT_DATA_DIR=$4

	rm -f -r ${INPUT_DATA_DIR}
	mkdir ${INPUT_DATA_DIR}

	if [[ ! -d ${INPUT_DATA_DIR} ]]
	then
		echo "Could not create data directoy \"${INPUT_DATA_DIR}\""
		exit
	fi

	SIMULATED_GWAS_DATA_GMMAT_DIR=${GMMAT_dir}
	covar_pheno_file=${SIMULATED_GWAS_DATA_GMMAT_DIR}/covar_pheno.txt
	geno_file=${SIMULATED_GWAS_DATA_GMMAT_DIR}/genotypes.txt

	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	echo "WARNING:"
	echo "ASSUMING THAT THE COVAR/PHENO (${covar_pheno_file}) AND GENOTYPE MATRIX (${geno_file}) HAVE MATCHING SAMPLE ORDERING."
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

	if [[ ! -f ${covar_pheno_file} ]]
	then
			echo "Could not find ${covar_pheno_file}"
			exit
	fi

	if [[ ! -f ${geno_file} ]]
	then
			echo "Could not find ${geno_file}"
			exit
	fi

	# Input data and model matrices.
	FEATS_MATRIX_FP=${INPUT_DATA_DIR}/feat_matrix.txt
	PHENO_FP=${INPUT_DATA_DIR}/phenotypes.txt
	MODEL_FP=${INPUT_DATA_DIR}/model.txt
	GENO_FP=${INPUT_DATA_DIR}/genotypes.txt

	echo "Creating GWAS data using GMMAT_dir=${GMMAT_dir}"
	#rm -f -r SIMULATED_GWAS_DATA
	#mkdir SIMULATED_GWAS_DATA
	#INPUT_DATA_DIR=SIMULATED_GWAS_DATA

	head -n ${N_VARS_2_USE} ${geno_file} > ${INPUT_DATA_DIR}/genotypes.txt

	pheno_col_id=`head -n 1 ${covar_pheno_file} | awk {'print $NF'}`
	if [[ "$pheno_col_id" != "PHENOTYPE" ]]
	then
		echo "Could not match the last column id in ${covar_pheno_file} to \"PHENOTYPE\", last column is expected to be the phenotype."
		exit
	fi

	echo "Writing phenotypes and covariates."
	awk {'print $1"\t"$NF'} ${covar_pheno_file} > ${INPUT_DATA_DIR}/phenotypes.txt
	awk {'printf("%s", $1);if(NR==1){printf("\tINTERCEPT");}else{printf("\t1");};for(i=2;i<NF;i++){printf("\t%s", $i);};print "";'} ${covar_pheno_file} > ${INPUT_DATA_DIR}/feat_matrix.txt
	head -n 1 ${INPUT_DATA_DIR}/feat_matrix.txt | awk {'print "N_COVARS\t"NF-2;for(i=2;i<=NF;i++){print $i"\t0.00"}'} > ${INPUT_DATA_DIR}/model.txt

	echo "Copying genotypes."
	cp ${geno_file} ${INPUT_DATA_DIR}/genotypes.txt

	exit
fi

if [[ ${cmd_option} == "-split_GMMAT_dir" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Pooled data GMMAT directory] [Per site output base directory]"
		exit
	fi

	source ${data_config_file}

	GMMAT_dir=$3
	INPUT_DATA_DIR=$4

	echo "Simulating data for ${N_SITES} using ${N_VARS_2_USE} variants."

	rm -f -r ${INPUT_DATA_DIR}
	mkdir ${INPUT_DATA_DIR}

	SIMULATED_GWAS_DATA_GMMAT_DIR=${GMMAT_dir}
	covar_pheno_file=${SIMULATED_GWAS_DATA_GMMAT_DIR}/covar_pheno.txt
	geno_file=${SIMULATED_GWAS_DATA_GMMAT_DIR}/genotypes.txt

	if [[ ! -f ${covar_pheno_file} ]]
	then
		echo "Could not find ${covar_pheno_file}"
		exit
	fi

	if [[ ! -f ${geno_file} ]]
	then
		echo "Could not find ${geno_file}"
		exit
	fi

    # Input data and model matrices.
    FEATS_MATRIX_FP=${INPUT_DATA_DIR}/feat_matrix.txt
    PHENO_FP=${INPUT_DATA_DIR}/phenotypes.txt
    MODEL_FP=${INPUT_DATA_DIR}/model.txt
    GENO_FP=${INPUT_DATA_DIR}/genotypes.txt
	
	echo "Creating GWAS data using GMMAT_dir=${GMMAT_dir}"
	#rm -f -r SIMULATED_GWAS_DATA
	#mkdir SIMULATED_GWAS_DATA
	#INPUT_DATA_DIR=SIMULATED_GWAS_DATA

	head -n ${N_VARS_2_USE} ${geno_file} > ${INPUT_DATA_DIR}/genotypes.txt

	echo "Writing phenotypes and covariates."
	awk {'print $1"\t"$NF'} ${covar_pheno_file} > ${INPUT_DATA_DIR}/phenotypes.txt
	awk {'printf("%s", $1);if(NR==1){printf("\tINTERCEPT");}else{printf("\t1");};for(i=2;i<NF;i++){printf("\t%s", $i);};print "";'} ${covar_pheno_file} > ${INPUT_DATA_DIR}/feat_matrix.txt
    head -n 1 ${INPUT_DATA_DIR}/feat_matrix.txt | awk {'print "N_COVARS\t"NF-2;for(i=2;i<=NF;i++){print $i"\t0.00"}'} > ${INPUT_DATA_DIR}/model.txt

	# Separate the data into N_SITES sites.
	awk {'if(NR>1){print $0}'} ${PHENO_FP} | cut -f1,1 | shuf > shuffled_subject_ids.txt
	awk {'if(NR>1){print $0}'} ${PHENO_FP} | cut -f1,1 > GMMAT_sample_ids.list

	site_i_list=`seq 1 $N_SITES`
	for site_i in ${site_i_list[@]}
	do
		mkdir ${INPUT_DATA_DIR}/SITE_${site_i}
	        echo "Separating data for site ${site_i}"
		split -n l/${site_i}/${N_SITES} shuffled_subject_ids.txt > ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list

	        # Separate phenotypes
		${COLLAGENE_SECURE_EXEC} -extract_rows_per_query_column_preserve_query_order ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/phenotypes.txt 0 1 temp_pheno.txt
		echo -e "SAMPLE\tPHENOTYPE" > pheno_header.txt
		cat pheno_header.txt temp_pheno.txt > ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt

	        # Separate covariates.
		${COLLAGENE_SECURE_EXEC} -extract_rows_per_query_column_preserve_query_order ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/feat_matrix.txt 0 1 temp_feats.txt
		head -n 1 ${INPUT_DATA_DIR}/feat_matrix.txt > feat_header.txt
		cat feat_header.txt temp_feats.txt > ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt

		# Separate genotypes.
		echo "Extracting genotypes for site ${site_i}"
		${COLLAGENE_SECURE_EXEC} -extract_subsample_GMMAT_genotypes_from_GMMAT_genotype_matrix  ${INPUT_DATA_DIR}/genotypes.txt GMMAT_sample_ids.list ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt
	done

	exit 0
fi

if [[ "${cmd_option}" == "-generate_keys" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [data config file] [Output directory]"
		exit 1
	fi

	data_config_file=$2
	KEYS_OP_DIR=$3

	source ${data_config_file}

	echo "Generating ${N_SITES} keys using text params in ${TEXT_PARAMS_PATH}"

	rm -f -r ${KEYS_OP_DIR}
	mkdir ${KEYS_OP_DIR}
	${COLLAGENE_SECURE_EXEC} -generate_per_site_noisy_keys ${N_SITES} ${TEXT_PARAMS_PATH} ${KEYS_OP_DIR}

	exit 0
fi

SIGMOID_TYPE=NATIVE

if [[ "${cmd_option}" == "-validate_input_data" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Local data directory]"
		exit 1
	fi

	data_config_file=$2
	LOCAL_DATA_DIR=$3
	
	rm -f INPUT.ERROR

	collagene_check=`type -P ${COLLAGENE_SECURE_EXEC}`
	if [[ "${collagene_check}" == "" ]]
	then
		echo "Could not find COLLAGENE executable."
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi
	
	if [[ ! -d ${LOCAL_DATA_DIR} ]]
	then
		echo "Could not find local data directory @ ${LOCAL_DATA_DIR}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
			mkdir ${INTERM_DATA_DIR}

			if [[ ! -d ${INTERM_DATA_DIR} ]]
			then
				echo "Could not create the intermediate directory @ \"${INTERM_DATA_DIR}\""
				exit 1
			fi
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt
	LOCAL_GENO_FILE=${LOCAL_DATA_DIR}/genotypes.txt

	if [[ ! -f ${LOCAL_FEAT_FILE} ]]
	then
		echo "Could not find features file @ \"${LOCAL_FEAT_FILE}\""
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${LOCAL_GENO_FILE} ]]
	then
		echo "Could not find genotype file @ \"${LOCAL_GENO_FILE}\""
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${LOCAL_PHENO_FILE} ]]
	then
		echo "Could not find phenotypes file @ \"${LOCAL_PHENO_FILE}\""
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PUBLIC_KEY_FILE} ]]
	then
	echo "Could not find public keys file @ ${PUBLIC_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${RELIN_KEY_FILE} ]]
	then
		echo "Could not find the relin. keys file @ ${RELIN_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${GALOIS_KEY_FILE} ]]
	then
		echo "Could not find the Galois keys file @ ${GALOIS_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PRIVATE_KEY_FILE} ]]
	then
		echo "Could not find the private key file @ ${PRIVATE_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PARTDEC_SYMKEY_FILE} ]]
	then
		echo "Could not find the symmetric encryption key file @ ${PARTDEC_SYMKEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find the security parameters file @ ${TEXT_PARAMS_PATH}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi
	
	# Verify the number of variants to score.
	N_VARS_IN_GENO=`wc -l ${LOCAL_GENO_FILE} | awk {'print $1'}`
	VAR_N_VALID=`echo -e ${N_VARS_IN_GENO}"\t"${N_VARS_2_SCORE}"\t"${VAR_BLOCK_SIZE} | awk {'if($1 < ($2+5)){print "0"}else{print "1"}'}`
	
	if [[ ${VAR_N_VALID} == 0 ]]
	then
		echo "Number of variants is not valid: Total in geno: ${N_VARS_IN_GENO}; Total to score: ${N_VARS_2_SCORE}; Var. block size: ${VAR_BLOCK_SIZE}"
		echo "ERROR" > INPUT.ERROR
		exit 1	
	fi

	# Test a file transfer.
	echo "Testing a file transfer."
	${FILE_IO_UTILS_SCRIPT} -test_file_IO ${data_config_file}
	if [[ $? != 0 ]]
	then
		echo "File transfer test failed."
		echo "ERROR" > INPUT.ERROR
	fi

	# Exit with no errors.
	exit 0
fi # -validate_input_data option.

#-client_calculate_save_XtWX_XtWz [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
if [[ "${cmd_option}" == "-client_calculate_save_XtWX_XtWz" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	if [[ ! -f ${LOCAL_FEAT_FILE} ]]
	then
                echo "Could not find features file @ \"${LOCAL_FEAT_FILE}\""
                exit
        fi

	if [[ ! -f ${LOCAL_PHENO_FILE} ]]
	then
		echo "Could not find phenotypes file @ \"${LOCAL_PHENO_FILE}\""
		exit
	fi

	source ${data_config_file}

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} Calculating XtWX and XtWz at epoch=${cur_epoch} for client ${site_i_per_cmd}"
	${COLLAGENE_SECURE_EXEC} -cryptable_client_calculate_save_XtWX_XtWz ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${LOCAL_FEAT_FILE} ${LOCAL_PHENO_FILE} ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INTERM_DATA_DIR} ${INTERM_DATA_DIR} >& ${INTERM_DATA_DIR}/XtWX_XtWz_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op

	XtWz_FILE=${INTERM_DATA_DIR}/XtWz_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin
	XtWz_PADDED_FILE=${INTERM_DATA_DIR}/XtWz_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_padded.bin
	XtWz_PADDED_ENC_FILE=${INTERM_DATA_DIR}/XtWz_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_padded.bin.enc
	NOISE_MATRIX_FILE=${INTERM_DATA_DIR}/mult_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}.bin
	NOISE_MATRIX_PADDED_FILE=${INTERM_DATA_DIR}/mult_padded_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}.bin

	${COLLAGENE_SECURE_EXEC} -write_plain_matrix_dimensions ${XtWz_FILE} ${INTERM_DATA_DIR}/XtWz_dims.txt
	n_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/XtWz_dims.txt`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} Found ${n_covars} covars in XtWz."

	# Wait for updating the RNG seed.
	sleep 1

	# Replace the noise from COLLAGENE with full noise.
	rm -f ${NOISE_MATRIX_FILE}
	${COLLAGENE_SECURE_EXEC} -generate_mult_full_noise_matrix ${n_covars} ${n_covars} ${NOISE_MATRIX_FILE}
	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${NOISE_MATRIX_FILE} ${NOISE_MATRIX_FILE}.txt

	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${XtWz_FILE} ${XtWz_PADDED_FILE}
	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${NOISE_MATRIX_FILE} ${NOISE_MATRIX_PADDED_FILE}
	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${XtWz_PADDED_FILE} ${INTERM_DATA_DIR}/PADDED_NOISE_DIMS.txt

	# Encrypt padded XtWz matrix.
	${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${XtWz_PADDED_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${XtWz_PADDED_ENC_FILE}

	# Original noise is necessary for inverting.
	#COL_EXP_MULT_NOISE_DIR=${INTERM_DATA_DIR}/col_exp_mult_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}
	#rm -f -r ${COL_EXP_MULT_NOISE_DIR}
	#mkdir ${COL_EXP_MULT_NOISE_DIR}
	ROW_EXP_MULT_NOISE_DIR=${INTERM_DATA_DIR}/row_exp_mult_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}
	rm -f -r ${ROW_EXP_MULT_NOISE_DIR}
	mkdir ${ROW_EXP_MULT_NOISE_DIR}

	#${COLLAGENE_SECURE_EXEC} -col_expand_dense_encrypt_matrix ${NOISE_MATRIX_FILE} ${n_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${COL_EXP_MULT_NOISE_DIR}
	${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${NOISE_MATRIX_FILE} ${n_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${ROW_EXP_MULT_NOISE_DIR}

	# Padded noise is necessary for padding XtWX after inversion.
	COL_EXP_PADDED_MULT_NOISE_DIR=${INTERM_DATA_DIR}/col_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}
	rm -f -r ${COL_EXP_PADDED_MULT_NOISE_DIR}
	mkdir ${COL_EXP_PADDED_MULT_NOISE_DIR}
	#ROW_EXP_PADDED_MULT_NOISE_DIR=${INTERM_DATA_DIR}/row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${site_i_per_cmd}
	#rm -f -r ${ROW_EXP_PADDED_MULT_NOISE_DIR}
	#mkdir ${ROW_EXP_PADDED_MULT_NOISE_DIR}

	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${NOISE_MATRIX_PADDED_FILE} ${INTERM_DATA_DIR}/padded_dims.txt
	n_padded_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/padded_dims.txt`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} Found ${n_padded_covars} padded covars in noise matrix."

	${COLLAGENE_SECURE_EXEC} -col_expand_dense_encrypt_matrix ${NOISE_MATRIX_PADDED_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${COL_EXP_PADDED_MULT_NOISE_DIR}
	#${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${NOISE_MATRIX_PADDED_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${ROW_EXP_PADDED_MULT_NOISE_DIR}

	# SEND ROW/COL EXPANSIONS OF NOISE and enc(XtWz) TO SHARED SPACE
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} Sending row/column expansions of noise and enc(XtWz) to shared space for epoch=${cur_epoch} at client ${site_i_per_cmd}"
	#cp -r ${XtWz_PADDED_ENC_FILE} ${COL_EXP_MULT_NOISE_DIR} ${ROW_EXP_MULT_NOISE_DIR} ${COL_EXP_PADDED_MULT_NOISE_DIR} ${ROW_EXP_PADDED_MULT_NOISE_DIR} ${SHARED_DIR} 
	#echo ${XtWz_PADDED_ENC_FILE} ${COL_EXP_MULT_NOISE_DIR} ${ROW_EXP_MULT_NOISE_DIR} ${COL_EXP_PADDED_MULT_NOISE_DIR} ${ROW_EXP_PADDED_MULT_NOISE_DIR} > ${INTERM_DATA_DIR}/TEMP_FILES.list
	echo ${XtWz_PADDED_ENC_FILE} ${ROW_EXP_MULT_NOISE_DIR} ${COL_EXP_PADDED_MULT_NOISE_DIR} > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list

	exit
fi

#-add_noise_to_XtWX [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]
# Copy the col and row expansions of the encrypted noise levels from all sites, add them up, multiply them with current site's XtWX and send the encrypted result to shared space -- All in ${COLLAGENE_SECURE_EXEC}.
if [[ ${cmd_option} == "-add_noise_to_XtWX" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi

	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	if [[ ! -f ${LOCAL_FEAT_FILE} ]]
	then
		echo "Could not find features file @ \"${LOCAL_FEAT_FILE}\""
		exit
	fi

	if [[ ! -f ${LOCAL_PHENO_FILE} ]]
	then
		echo "Could not find phenotypes file @ \"${LOCAL_PHENO_FILE}\""
		exit
	fi

	source ${data_config_file}

	# Download the ROW/COL Expansions of the noise levels.
	PER_SITE_NOISE_EXP_DIR=${INTERM_DATA_DIR}/PER_SITE_NOISE_EXP_ITER_${cur_epoch}
	PER_SITE_PADDED_NOISE_EXP_DIR=${INTERM_DATA_DIR}/PER_SITE_PADDED_NOISE_EXP_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_NOISE_EXP_DIR}
	rm -f -r ${PER_SITE_PADDED_NOISE_EXP_DIR}
	mkdir ${PER_SITE_NOISE_EXP_DIR}
	mkdir ${PER_SITE_PADDED_NOISE_EXP_DIR}
	
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in local feature matrix on client ${site_i_per_cmd}.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES1.list
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES2.list
	for site_i_f in ${site_iters[@]}
	do
		#echo col_exp_mult_noise_matrix_iter_${cur_epoch}_client_${site_i_f} >> ${INTERM_DATA_DIR}/TEMP_FILES1.list
		echo row_exp_mult_noise_matrix_iter_${cur_epoch}_client_${site_i_f} >> ${INTERM_DATA_DIR}/TEMP_FILES1.list

		echo col_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${site_i_f} >> ${INTERM_DATA_DIR}/TEMP_FILES2.list	
		#echo row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${site_i_f} >> ${INTERM_DATA_DIR}/TEMP_FILES2.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES1.list
	cur_status=$?
	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES2.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES1.list ${PER_SITE_NOISE_EXP_DIR}
	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES2.list ${PER_SITE_PADDED_NOISE_EXP_DIR}
	#cp -r ${SHARED_DIR}/col_exp_mult_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_NOISE_EXP_DIR}
    #cp -r ${SHARED_DIR}/col_exp_mult_padded_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_PADDED_NOISE_EXP_DIR}
	###############################################################################################################################################

	POOLED_ROW_EXP_NOISE_DIR=${INTERM_DATA_DIR}/pooled_row_exp_mult_noise_matrix_iter_${cur_epoch}
	#POOLED_COL_EXP_NOISE_DIR=${INTERM_DATA_DIR}/pooled_col_exp_mult_noise_matrix_iter_${cur_epoch}
    rm -f -r ${POOLED_ROW_EXP_NOISE_DIR}
    mkdir ${POOLED_ROW_EXP_NOISE_DIR}
    #rm -f -r ${POOLED_COL_EXP_NOISE_DIR}
    #mkdir ${POOLED_COL_EXP_NOISE_DIR}

    #POOLED_ROW_EXP_PADDED_NOISE_DIR=${INTERM_DATA_DIR}/pooled_row_exp_mult_padded_noise_matrix_iter_${cur_epoch}
    POOLED_COL_EXP_PADDED_NOISE_DIR=${INTERM_DATA_DIR}/pooled_col_exp_mult_padded_noise_matrix_iter_${cur_epoch}
	#rm -f -r ${POOLED_ROW_EXP_PADDED_NOISE_DIR}
	#mkdir ${POOLED_ROW_EXP_PADDED_NOISE_DIR}
	rm -f -r ${POOLED_COL_EXP_PADDED_NOISE_DIR}
	mkdir ${POOLED_COL_EXP_PADDED_NOISE_DIR}

	for cur_row_i in ${row_col_iters[@]}
	do
		rm -f ${INTERM_DATA_DIR}/temp_row_exp_files.list
		#rm -f ${INTERM_DATA_DIR}/temp_col_exp_files.list
		for cur_site_i in ${site_iters[@]}
		do
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Pooling site ${cur_site_i}'s noise levels for row/col ${cur_row_i}"
			echo "${PER_SITE_NOISE_EXP_DIR}/row_exp_mult_noise_matrix_iter_${cur_epoch}_client_${cur_site_i}/reprow_${cur_row_i}.bin.enc" >> ${INTERM_DATA_DIR}/temp_row_exp_files.list
			#echo "${PER_SITE_NOISE_EXP_DIR}/col_exp_mult_noise_matrix_iter_${cur_epoch}_client_${cur_site_i}/repcol_${cur_row_i}.bin.enc" >> ${INTERM_DATA_DIR}/temp_col_exp_files.list
		done

		# Summate all noise levels from the list of files.
		${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_row_exp_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_ROW_EXP_NOISE_DIR}/reprow_${cur_row_i}.bin.enc >& ${INTERM_DATA_DIR}/noise_row_exp_pooling.txt

		## Summate all noise levels from the list of files.
		#${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_col_exp_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_COL_EXP_NOISE_DIR}/repcol_${cur_row_i}.bin.enc >& ${INTERM_DATA_DIR}/noise_col_exp_pooling.txt
	done

	# Pool the padded noise levels: This is used later after pooling all noisy X'WX from all sites
	for cur_row_i in ${row_col_iters[@]}
	do
		#rm -f ${INTERM_DATA_DIR}/temp_row_exp_files.list
		rm -f ${INTERM_DATA_DIR}/temp_col_exp_files.list
		for cur_site_i in ${site_iters[@]}
		do
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Pooling site ${cur_site_i}'s padded noise levels for row/col ${cur_row_i}"
			#echo "${PER_SITE_PADDED_NOISE_EXP_DIR}/row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${cur_site_i}/reprow_${cur_row_i}.bin.enc" >> ${INTERM_DATA_DIR}/temp_row_exp_files.list
			echo "${PER_SITE_PADDED_NOISE_EXP_DIR}/col_exp_mult_padded_noise_matrix_iter_${cur_epoch}_client_${cur_site_i}/repcol_${cur_row_i}.bin.enc" >> ${INTERM_DATA_DIR}/temp_col_exp_files.list
		done

		## Summate all noise levels from the list of files.
		#${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_row_exp_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_ROW_EXP_PADDED_NOISE_DIR}/reprow_${cur_row_i}.bin.enc >& ${INTERM_DATA_DIR}/padded_noise_row_exp_pooling.txt

		# Summate all noise levels from the list of files.
		${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_col_exp_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_COL_EXP_PADDED_NOISE_DIR}/repcol_${cur_row_i}.bin.enc >& ${INTERM_DATA_DIR}/padded_noise_col_exp_pooling.txt
	done

	# Row/col expand X'WX matrix for this site.
	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${INTERM_DATA_DIR}/XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin ${INTERM_DATA_DIR}/XtWX_dims.txt
	n_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/XtWX_dims.txt`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars} covars in XtWX."

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Col/Row expanding X'WX"

	XtWX_COL_EXP_DIR=${INTERM_DATA_DIR}/col_exp_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${XtWX_COL_EXP_DIR}
	mkdir $XtWX_COL_EXP_DIR
	${COLLAGENE_SECURE_EXEC} -col_expand_dense_encrypt_matrix ${INTERM_DATA_DIR}/XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin ${n_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${XtWX_COL_EXP_DIR} >& ${INTERM_DATA_DIR}/XtWX_col_expansion.op

	XtWX_ROW_EXP_DIR=${INTERM_DATA_DIR}/row_exp_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${XtWX_ROW_EXP_DIR}
	mkdir $XtWX_ROW_EXP_DIR
	${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${INTERM_DATA_DIR}/XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin ${n_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${XtWX_ROW_EXP_DIR} >& ${INTERM_DATA_DIR}/XtWX_row_expansion.op

	# Multiply XtWX with pooled noise using expanded data.
	${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${XtWX_COL_EXP_DIR} ${POOLED_ROW_EXP_NOISE_DIR} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc

	###############################################################################################################################################
	# This upload strips the file names to the basename and uses them to upload, should be fine to upload like this.
	#cp ${INTERM_DATA_DIR}/noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc ${SHARED_DIR}
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading noisy XtWX to shared space."
	echo ${INTERM_DATA_DIR}/noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	###############################################################################################################################################

	exit
fi

#-pool_site_specific_all_site_noise_XtWX
# Pool all noisy and encrypted XtWX from all sites and send it back to shared space.
if [[ ${cmd_option} == "-pool_site_specific_all_site_noise_XtWX" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	n_covars=`awk {'print NF-1'} ${LOCAL_FEAT_FILE} | head -n 1`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in feature matrix.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`


	if [[ ! -f ${LOCAL_FEAT_FILE} ]]
	then
		echo "Could not find features file @ \"${LOCAL_FEAT_FILE}\""
		exit
	fi

	if [[ ! -f ${LOCAL_PHENO_FILE} ]]
	then
		echo "Could not find phenotypes file @ \"${LOCAL_PHENO_FILE}\""
		exit
	fi

	source ${data_config_file}

	rm -f -r ${INTERM_DATA_DIR}/NOISY_PER_SITE_XtWX
	mkdir ${INTERM_DATA_DIR}/NOISY_PER_SITE_XtWX
	#cp -r ${SHARED_DIR}/noisy_XtWX_CLIENT_*_ITER_${cur_epoch}.bin.enc ${INTERM_DATA_DIR}/NOISY_PER_SITE_XtWX
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo noisy_XtWX_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${INTERM_DATA_DIR}/NOISY_PER_SITE_XtWX
	###############################################################################################################################################

	ls ${INTERM_DATA_DIR}/NOISY_PER_SITE_XtWX/* > ${INTERM_DATA_DIR}/temp_noisy_XtWX_files.list

	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_noisy_XtWX_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc >& ${INTERM_DATA_DIR}/noisy_XtWX_pooling.txt

	# Upload to online shared directory.
	#cp ${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc ${SHARED_DIR}
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading pooled onoist XtWX to shared space."
	echo ${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list

	exit 0
fi

#-collaborative_decrypt_pooled_noisy_XtWX
# Collaboratively decrypt all site's sitewise pooled noisy XtWX copies and send them back to shared space.
if [[ ${cmd_option} == "-collaborative_decrypt_pooled_noisy_XtWX" ]]
then
	if [[ $# != 6 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Cur site secret key] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	private_key_share_file=$3
	cur_epoch=$4
	site_i_per_cmd=$5
	LOCAL_DATA_DIR=$6

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi
	
	# Take all other site's XtWX and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in features matrix for client ${site_i_per_cmd} @ epoch: ${cur_epoch}.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	rm -f -r ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX
	mkdir ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX
	#cp -r ${SHARED_DIR}/pooled_noisy_XtWX_CLIENT_*_ITER_${cur_epoch}.bin.enc ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo pooled_noisy_XtWX_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX
	###############################################################################################################################################

	is_site0=0
	if [[ ${site_i_per_cmd} == 0 ]]
	then
		is_site0=1
	fi

	echo "is_site0=${is_site0}"

	# Create directory for pooled noise levels.
	#for cur_XtWX_file in ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX/*
	#do
	for site_i_f in ${site_iters[@]}
	do
		# This file is checked and downloaded above.
		cur_XtWX_file=${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX/pooled_noisy_XtWX_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin.enc
		if [[ ! -f ${cur_XtWX_file} ]]
		then
			echo "Could not find ${cur_XtWX_file} for partial decryption.."
			exit 1
		fi

		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Partially decrypting ${cur_XtWX_file} using Site-${site_i_per_cmd}'s key."

		cur_XtWX_file_partdec_file=${cur_XtWX_file}_partdec_by_${site_i_per_cmd}.partdec
		${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_XtWX_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_XtWX_file_partdec_file}

		# TODO: Use the symmetric key to encrypt the partially decrypted data.
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Symmetric encrypting ${cur_XtWX_file_partdec_file} using shared symetric key."

		cur_XtWX_file_partdec_enc_file=${cur_XtWX_file}_partdec_by_${site_i_per_cmd}.partdec.enc
		$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_XtWX_file_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_XtWX_file_partdec_enc_file}
	done

	###############################################################################################################################################
	# Copy the partially decrypted noisy XtWX.
	#cp ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX/*.partdec ${SHARED_DIR}
	#ls ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX/*.partdec > ${INTERM_DATA_DIR}/TEMP_FILES.list
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading partially decrypted pooled noisy XtWX to shared space from client ${site_i_per_cmd}"
	ls ${INTERM_DATA_DIR}/PER_SITE_POOLED_XtWX/*.partdec.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	###############################################################################################################################################
	exit
fi

#-pool_partially_decrypted_pooled_noisy_XtWx_remove_noise
# Pool the partial decryptions to generate plaintext noisy XtWX, invert it, encrypt it, remove additive noise by multiplying with col/row expansion of the noise levels.
if [[ ${cmd_option} == "-pool_partially_decrypted_pooled_noisy_XtWx_remove_noise" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	# Take all other site's XtWX and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	n_covars=`awk {'print NF-1'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates on client ${site_i_per_cmd} @ epoch=${cur_epoch}.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_ENC_PARTDEC_DIR=${INTERM_DATA_DIR}/PER_SITE_ENC_PARTDEC_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_ENC_PARTDEC_DIR}
    mkdir ${PER_SITE_ENC_PARTDEC_DIR} 
    #cp -r ${SHARED_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_*.partdec ${PER_SITE_PARTDEC_DIR}
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		#echo pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec >> ${INTERM_DATA_DIR}/TEMP_FILES.list	
		echo pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list	
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_ENC_PARTDEC_DIR}
	###############################################################################################################################################

	# Decrypt the partial decrypted files.
	PER_SITE_PARTDEC_DIR=${INTERM_DATA_DIR}/PER_SITE_PARTDEC_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_PARTDEC_DIR}
    mkdir ${PER_SITE_PARTDEC_DIR} 
	for site_i_f in ${site_iters[@]}
	do
		cur_enc_partdec_file=${PER_SITE_ENC_PARTDEC_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec.enc
		cur_partdec_file=${PER_SITE_PARTDEC_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec
		$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}
	done

	ls ${PER_SITE_PARTDEC_DIR}/*.partdec > ${INTERM_DATA_DIR}/temp_partdecs.list

	# TODO:: These need to be renamed to reflect that they are full decrypted noisy XtWX matrices.
	POOLED_NOISY_XtWX_FILE=${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin
	INV_NOISY_XtWX_FILE=${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_inv.bin
	INV_NOISY_XtWX_PADDED_FILE=${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_inv.bin_padded.bin

	${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_partdecs.list ${TEXT_PARAMS_PATH} ${POOLED_NOISY_XtWX_FILE}

	# Invert.
	${COLLAGENE_SECURE_EXEC} -plain_invert_matrix ${POOLED_NOISY_XtWX_FILE} ${INV_NOISY_XtWX_FILE}

	# Scale Z by Z-scaler.
	RESCALED_INV_NOISY_XtWX_FILE=${INV_NOISY_XtWX_FILE}_rescaled.bin
	${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${INV_NOISY_XtWX_FILE} ${Z_SCALER} ${RESCALED_INV_NOISY_XtWX_FILE}

	# Also do 2^n expansion using padding: We need the padding of Z for row2row multiplication to XtWz while updating beta.
	# The signal level is boosted above to ensure beta updates are robust against ct noise in 0 entries of padding.
	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${RESCALED_INV_NOISY_XtWX_FILE} ${INV_NOISY_XtWX_PADDED_FILE}

	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${INV_NOISY_XtWX_PADDED_FILE} ${INTERM_DATA_DIR}/padded_dims.txt
	n_padded_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/padded_dims.txt`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_padded_covars} padded covars in inverted noisy XtWX matrix."

	# Row expand the inverted noisy XtWX matrix and the padded version.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Row expanding noisy inv(XtWX)"
	POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_DIR=${INTERM_DATA_DIR}/POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_DIR=${INTERM_DATA_DIR}/POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}

	rm -f -r ${POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_DIR} ${POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_DIR}
	mkdir ${POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_DIR} ${POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_DIR}
	#${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${INV_NOISY_XtWX_FILE} ${n_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_DIR}
	${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${INV_NOISY_XtWX_PADDED_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_DIR}

	## Multiply with noise: Note that this removes the noise securely from the unpadded matrix, which we do not use at all.
	#${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${INTERM_DATA_DIR}/pooled_col_exp_mult_noise_matrix_iter_${cur_epoch} ${POOLED_INV_NOISY_XtWX_ROW_EXPANSION_CLIENT_DIR} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc

	# Expand padded inv(noisy_XtWX) and multiply with padded noise levels; these are computed in the previous step.
	POOLED_ROW_EXP_PADDED_NOISE_DIR=${INTERM_DATA_DIR}/pooled_row_exp_mult_padded_noise_matrix_iter_${cur_epoch}
	POOLED_COL_EXP_PADDED_NOISE_DIR=${INTERM_DATA_DIR}/pooled_col_exp_mult_padded_noise_matrix_iter_${cur_epoch}

	# Multiple the padded XtWX on the left with padded noise.
	${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${POOLED_COL_EXP_PADDED_NOISE_DIR} ${POOLED_INV_NOISY_XtWX_PADDED_ROW_EXPANSION_CLIENT_DIR} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/padded_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc

	# Row expand the padded inverse XtWX matrix, this will be used in GWX matrix multiplication after pooling them from all sites.
	# Note that this multiplication is done on padded XtWX bc it will then be used to perform row2row multiplication with padded GWX again.
	# We therefore do not need unpadded Z any more.

	exit
fi

# At this stage, each site should have the encrypted inv(XtWX) matrix: Do a check on this matrix to confirm the encrypted operations.

# -pool_XtWX_XtWz_update_beta : Pool encrypted XtWz and update beta.
if [[ ${cmd_option} == "-pool_XtWX_XtWz_update_beta" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi

	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	# Take all other site's XtWX and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	n_covars=`awk {'print NF-1'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates at client ${site_i_per_cmd} @ epoch=${cur_epoch}.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_XtWz_DIR=${INTERM_DATA_DIR}/PER_SITE_XtWz_DIR_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	POOLED_PADDED_XtWz_FILE=${INTERM_DATA_DIR}/pooled_XtWz_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.enc
	rm -f -r ${PER_SITE_XtWz_DIR}
	mkdir ${PER_SITE_XtWz_DIR}
	#cp -r ${SHARED_DIR}/XtWz_CLIENT_*_ITER_${cur_epoch}.bin_padded.bin.enc ${PER_SITE_XtWz_DIR}
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo XtWz_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin_padded.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list	
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_XtWz_DIR}
	###############################################################################################################################################

	ls ${PER_SITE_XtWz_DIR}/* > ${INTERM_DATA_DIR}/temp_XtWz_files.list

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling encrypted padded XtWz from all sites"
	n_files=`wc -l ${INTERM_DATA_DIR}/temp_XtWz_files.list | awk {'print $1'}`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_files} XtWz matrices from all sites"
	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_XtWz_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_PADDED_XtWz_FILE}

	# XtWz is a px1 column vector. We can transpose it by just switching the sizes in the matrix file.
	# inv(XtWX) is a pxp square matrix and it is symmetric.
	inv_XtWX_file=${INTERM_DATA_DIR}/inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc
	inv_padded_XtWX_file=${INTERM_DATA_DIR}/padded_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc
	
	# First, encrypt the local XtWz file.
	enc_trans_PADDED_XtWz_FILE=${INTERM_DATA_DIR}/trans_XtWz_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc
	enc_BETA_FILE=${INTERM_DATA_DIR}/beta_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc

	${COLLAGENE_SECURE_EXEC} -transpose_continuous_encrypted_vector ${POOLED_PADDED_XtWz_FILE} ckks.params ${enc_trans_PADDED_XtWz_FILE}

    ${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${POOLED_PADDED_XtWz_FILE} ${INTERM_DATA_DIR}/padded_dims.txt
    n_padded_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/padded_dims.txt`
    date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_padded_covars} padded covars in XtWz matrix."

	ROW_EXP_TRANS_PADDED_XtWz=${INTERM_DATA_DIR}/ROW_EXP_TRANS_PADDED_XtWZ_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${ROW_EXP_TRANS_PADDED_XtWz}
	mkdir ${ROW_EXP_TRANS_PADDED_XtWz}
	${COLLAGENE_SECURE_EXEC} -row_expand_continuous_encrypted_matrix ${enc_trans_PADDED_XtWz_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${ROW_EXP_TRANS_PADDED_XtWz}

	# Do a row-by-row multiplication of the matrices to get beta.
	${COLLAGENE_SECURE_EXEC} -secure_row2row_inner_prod_continuous_encrypted_matrices ${ROW_EXP_TRANS_PADDED_XtWz}/reprow_0.bin.enc ${inv_padded_XtWX_file} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${enc_BETA_FILE}

	# Copy the beta file to shared directory.
	###############################################################################################################################################
	## This upload strips the file names to the basename and uses them to upload, should be fine to upload like this.
	#cp ${enc_BETA_FILE} ${SHARED_DIR}
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading beta estimate to shared space on client ${site_i_per_cmd} @ epoch=${cur_epoch}"
	echo ${enc_BETA_FILE} > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	###############################################################################################################################################	

	exit
fi

# At this point, we have the beta estimate that is encrypted, we need to decrypt it.
# -client_collaborative_decrypt_beta: Collaboratively decrypt beta and submit partial decryptions to shared space.
# -pool_XtWX_XtWz_update_beta : Pool encrypted XtWz and update beta.
if [[ ${cmd_option} == "-client_collaborative_decrypt_beta" ]]
then
	if [[ $# != 6 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Cur site secret key] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	private_key_share_file=$3
	cur_epoch=$4
	site_i_per_cmd=$5
	LOCAL_DATA_DIR=$6

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	# Take all other site's BETA's and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in local features matrix.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_BETA_DIR=${INTERM_DATA_DIR}/PER_SITE_BETA
	rm -f -r ${PER_SITE_BETA_DIR}
	mkdir ${PER_SITE_BETA_DIR}
	
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	#cp -r ${SHARED_DIR}/beta_CLIENT_*_ITER_${cur_epoch}.bin.enc ${PER_SITE_BETA_DIR}
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo beta_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list	
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_BETA_DIR}
	###############################################################################################################################################

	is_site0=0
	if [[ ${site_i_per_cmd} == 0 ]]
	then
		is_site0=1
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: is_site0=${is_site0}"

	# Create directory for pooled noise levels.
	#for cur_beta_file in ${PER_SITE_BETA_DIR}/*
	#do
	for site_i_f in ${site_iters[@]}
	do
		cur_beta_file=${PER_SITE_BETA_DIR}/beta_CLIENT_${site_i_f}_ITER_${cur_epoch}.bin.enc
		if [[ ! -f ${cur_beta_file} ]]
		then
			echo "Could not find ${cur_beta_file} for partial decryption."
			exit 1
		fi

		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Partially decrypting ${cur_beta_file} using Site-${site_i_per_cmd}'s key."

		cur_beta_file_partdec_file=${cur_beta_file}_partdec_by_${site_i_per_cmd}.partdec
		${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_beta_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_beta_file_partdec_file}

		# Use the symmetric key to encrypt the partially decrypted data.
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Symmetric encrypting ${cur_beta_file_partdec_file} using shared symetric key."

		cur_beta_file_partdec_enc_file=${cur_beta_file}_partdec_by_${site_i_per_cmd}.partdec.enc
		$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_beta_file_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_beta_file_partdec_enc_file}
	done

	# Copy the partially decrypted noisy beta.	
	###############################################################################################################################################
	## This upload strips the file names to the basename and uses them to upload, should be fine to upload like this.
	#cp ${PER_SITE_BETA_DIR}/*.partdec ${SHARED_DIR}
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading partially decrypted and symmetric encrypted beta on site ${site_i_per_cmd} @ epoch=${cur_epoch}"
	ls ${PER_SITE_BETA_DIR}/*.partdec.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	###############################################################################################################################################	

	exit
fi

# -client_pool_partially_decrypted_beta : Download partial decryptions for this site from shared space and pool them to generate plaintext beta.
# At this stage, we should have an updated beta parameter; check to make sure it is reasonable.
if [[ ${cmd_option} == "-client_pool_partially_decrypted_beta" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	# Take all other site's XtWX and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	n_covars=`awk {'print NF-1'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in local features file.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_ENC_PARTDEC_DIR=${INTERM_DATA_DIR}/PER_SITE_ENC_PARTDEC_BETA_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_ENC_PARTDEC_DIR}
	mkdir ${PER_SITE_ENC_PARTDEC_DIR}
	#cp -r ${SHARED_DIR}/beta_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_*.partdec ${PER_SITE_PARTDEC_DIR}
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	#cp -r ${SHARED_DIR}/beta_CLIENT_*_ITER_${cur_epoch}.bin.enc ${PER_SITE_BETA_DIR}
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo beta_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_ENC_PARTDEC_DIR}
	###############################################################################################################################################

	# Decrypt each partdec.
	PER_SITE_PARTDEC_DIR=${INTERM_DATA_DIR}/PER_SITE_PARTDEC_BETA_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_PARTDEC_DIR}
	mkdir ${PER_SITE_PARTDEC_DIR}

	for site_i_f in ${site_iters[@]}
	do
		cur_enc_partdec_file=${PER_SITE_ENC_PARTDEC_DIR}/beta_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec.enc
		cur_partdec_file=${PER_SITE_PARTDEC_DIR}/beta_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc_partdec_by_${site_i_f}.partdec
		$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}
	done

	# Get the list of decrypted partdecs and pool them.
	ls ${PER_SITE_PARTDEC_DIR}/* > ${INTERM_DATA_DIR}/temp_partdecs.list

	POOLED_PADDED_BETA_FILE=${INTERM_DATA_DIR}/pooled_padded_BETA_${site_i_per_cmd}_ITER_${cur_epoch}.bin

	${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_partdecs.list ${TEXT_PARAMS_PATH} ${POOLED_PADDED_BETA_FILE}
	${COLLAGENE_SECURE_EXEC} -plain_unpad_matrix_to_size ${POOLED_PADDED_BETA_FILE} ${n_covars} 1 ${INTERM_DATA_DIR}/full_dec_beta_${cur_epoch}_client_${site_i_per_cmd}_SCALED.bin

	Z_UNSCALER=`echo ${Z_SCALER} | awk '{print 1.0/$1}'`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Z_UNSCALER: ${Z_UNSCALER} on site ${site_i_per_cmd} @ epoch=${cur_epoch}"
	${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${INTERM_DATA_DIR}/full_dec_beta_${cur_epoch}_client_${site_i_per_cmd}_SCALED.bin ${Z_UNSCALER} ${INTERM_DATA_DIR}/full_dec_beta_${cur_epoch}_client_${site_i_per_cmd}.bin

	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${INTERM_DATA_DIR}/full_dec_beta_${cur_epoch}_client_${site_i_per_cmd}.bin ${INTERM_DATA_DIR}/full_dec_beta_${cur_epoch}_client_${site_i_per_cmd}.bin.txt

	exit
fi

# -check_convergence_per_updated_beta : This is done either by beta or by LL estimations.

# Move to next epoch...

################################################################################################################################################
# AFTER THIS, WE WILL DO P-VALUE ASSIGNMENTS. THESE USE PADDED inv(X'WX)=Z and G'WG, G'WX Z X'WG multiplications that require some operations in HE.
################################################################################################################################################

#-client_calculate_save_pvalue_stats [Data config file] [Client index (starts at 0)] [Local data directory]
if [[ ${cmd_option} == "-client_calculate_save_pvalue_stats" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	cur_epoch=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	if [[ "${Z_SCALER}" == "" ]]
	then
		echo "Must set ${Z_SCALER}."
		exit
	fi

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt
	LOCAL_GENO_FILE=${LOCAL_DATA_DIR}/genotypes.txt

	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${INTERM_DATA_DIR}/padded_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc ${INTERM_DATA_DIR}/padded_dims.txt
	n_padded_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/padded_dims.txt`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_padded_covars} padded covars in inv XtWX matrix."

	# Following code rescales noisy Z=inv(XtWX) before padding. This is necessary to decrease the noise levels in the padded row/col of Z.
	INV_NOISY_XtWX_FILE=${INTERM_DATA_DIR}/pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_inv.bin

	# Scale this.
	RESCALED_INV_NOISY_XtWX_FILE=${INTERM_DATA_DIR}/rescaled_pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_inv.bin
	${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${INV_NOISY_XtWX_FILE} ${Z_SCALER} ${RESCALED_INV_NOISY_XtWX_FILE}

	# Pad it.
	PADDED_RESCALED_INV_NOISY_XtWX_FILE=${INTERM_DATA_DIR}/padded_rescaled_pooled_noisy_XtWX_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}.bin_inv.bin
	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${RESCALED_INV_NOISY_XtWX_FILE} ${PADDED_RESCALED_INV_NOISY_XtWX_FILE}

	# Row expand rescaled-padded Z.
	PADDED_RESCALED_INV_XtWX_ROW_EXP_DIR=${INTERM_DATA_DIR}/Row_Exp_Padded_Rescaled_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PADDED_RESCALED_INV_XtWX_ROW_EXP_DIR}
	mkdir ${PADDED_RESCALED_INV_XtWX_ROW_EXP_DIR}
	${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${PADDED_RESCALED_INV_NOISY_XtWX_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${PADDED_RESCALED_INV_XtWX_ROW_EXP_DIR}

	# Remove mult noise.
	# Expand padded inv(noisy_XtWX) and multiply with padded noise levels; these are computed in the previous step.
	POOLED_COL_EXP_PADDED_NOISE_DIR=${INTERM_DATA_DIR}/pooled_col_exp_mult_padded_noise_matrix_iter_${cur_epoch}
	if [[ ! -d ${POOLED_COL_EXP_PADDED_NOISE_DIR} ]]
	then
		echo "Could not find directory \"${POOLED_COL_EXP_PADDED_NOISE_DIR}\""
		exit
	fi

	# Multiple the padded XtWX on the left with padded noise.
	${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${POOLED_COL_EXP_PADDED_NOISE_DIR} ${PADDED_RESCALED_INV_XtWX_ROW_EXP_DIR} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/padded_rescaled_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc

	CUR_invXtWX_FILE=${INTERM_DATA_DIR}/padded_rescaled_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc
	if [[ ! -f ${CUR_invXtWX_FILE} ]]
	then
		echo "Could not find rescaled inverse Z @ \"${CUR_invXtWX_FILE}\""
		exit
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Row expanding pad(inv(XtWX))"
	PADDED_INV_XtWX_ROW_EXP_DIR=${INTERM_DATA_DIR}/Row_Exp_padded_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PADDED_INV_XtWX_ROW_EXP_DIR}
	mkdir ${PADDED_INV_XtWX_ROW_EXP_DIR}
	${COLLAGENE_SECURE_EXEC} -row_expand_continuous_encrypted_matrix ${CUR_invXtWX_FILE} ${VAR_BLOCK_SIZE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${PADDED_INV_XtWX_ROW_EXP_DIR} >& ${INTERM_DATA_DIR}/padded_invXtWX_row_exp.op

	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${INTERM_DATA_DIR}/padded_inv_XtWX_${site_i_per_cmd}_ITER_${cur_epoch}.bin.enc ${INTERM_DATA_DIR}/padded_dims.txt
	n_padded_covars=`awk {'print $1'} ${INTERM_DATA_DIR}/padded_dims.txt`
	echo "Found ${n_padded_covars} padded covars."

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Calculating p-value stats.."

	# Run COLLAGENE to get p-value statistics.
	${COLLAGENE_SECURE_EXEC} -cryptable_client_calculate_save_pvalue_stats ${site_i_per_cmd} ${N_SITES} ${N_EPOCHS} ${VAR_BLOCK_SIZE} ${LOCAL_GENO_FILE} ${LOCAL_FEAT_FILE} ${LOCAL_PHENO_FILE} ${INTERM_DATA_DIR} ${INTERM_DATA_DIR} >& ${INTERM_DATA_DIR}/P_VAL_STAT_CLIENT_${site_i_per_cmd}.op

	# Encrypt and send the p-value statistics.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Encrypting G(Y-Y0), GtWX, GtWG, calculating GtWXZ"

        # Processing 10,000,000 variants max.
		#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
		requested_last_block_start_i=`$0 -get_last_win_start_i $2`
        block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
        for cur_block_start_i in ${block_start_i_list[@]}
        do
            # Set the end of the current variant block.
            cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`
            date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Encrypting/Processing var block [${cur_block_start_i}-${cur_block_end_i}]"
	
			GtWX_FILE=${INTERM_DATA_DIR}/GtWX_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin
			GtWX_PADDED_FILE=${INTERM_DATA_DIR}/GtWX_PADDED_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin
			GtWG_FILE=${INTERM_DATA_DIR}/GtWG_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin
			Gt_YminY0_FILE=${INTERM_DATA_DIR}/P_CHISQ_STAT_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin

			if [[ ! -f ${GtWX_FILE} || ! -f ${GtWG_FILE} || ! -f ${Gt_YminY0_FILE} ]]
			then
				date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
				echo "${date_time_str} @ $1: Could not find one of the files:"
				ls -sort ${GtWX_FILE}
				ls -sort ${GtWG_FILE}
				ls -sort ${Gt_YminY0_FILE}
				exit
			fi

		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Padding GtWX"
		${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${GtWX_FILE} ${GtWX_PADDED_FILE} >& ${INTERM_DATA_DIR}/GtWX_${cur_block_start_i}_${cur_block_end_i}_padding.op

		if [[ ! -f ${GtWX_PADDED_FILE} ]]
		then
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Could not find:"
			ls -sort ${GtWX_PADDED_FILE}
			exit
		fi

		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Encrypting GtWX-padded, GtWG, Gt(Y-Y0)"
		${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${GtWX_PADDED_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${GtWX_PADDED_FILE}.enc >& ${INTERM_DATA_DIR}/GTWX_pad_${cur_block_start_i}_${cur_block_end_i}_enc.op

		# GtWG must also be rescaled to the same scale as Z since it is used in computing T.
		RESCALED_GtWG_FILE=${INTERM_DATA_DIR}/rescaled_GtWG_${cur_block_start_i}_${cur_block_end_i}.bin
		${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${GtWG_FILE} ${Z_SCALER} ${RESCALED_GtWG_FILE}
        ${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${RESCALED_GtWG_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${GtWG_FILE}.enc >& ${INTERM_DATA_DIR}/GtWG_${cur_block_start_i}_${cur_block_end_i}_enc.op

        ${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${Gt_YminY0_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${Gt_YminY0_FILE}.enc >& ${INTERM_DATA_DIR}/Gt_YminY0_${cur_block_start_i}_${cur_block_end_i}_enc.op
		wait

		if [[ ! -f ${GtWX_PADDED_FILE} || ! -f ${GtWG_FILE}.enc || ! -f ${Gt_YminY0_FILE}.enc ]]
		then
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Could not find one of the files:"
			ls -sort ${GtWX_PADDED_FILE} ${GtWG_FILE} ${Gt_YminY0_FILE}
			exit
		fi

		# Multiply G'WX with Z=inv(X'WX) using the row expansion of Z that is previously calculated.
		# We first need the col expansion of G'WX.
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Column expanding padded GtWX"
		GtWX_VAR_BLOCK_COL_EXP_DIR=${INTERM_DATA_DIR}/GtWX_PADDED_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_COL_EXP_CLIENT_${site_i_per_cmd}
		rm -f -r ${GtWX_VAR_BLOCK_COL_EXP_DIR}
		mkdir ${GtWX_VAR_BLOCK_COL_EXP_DIR}
		${COLLAGENE_SECURE_EXEC} -col_expand_dense_encrypt_matrix ${GtWX_PADDED_FILE} ${n_padded_covars} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${GtWX_VAR_BLOCK_COL_EXP_DIR} >& ${INTERM_DATA_DIR}/GtWX_${cur_block_start_i}_${cur_block_end_i}_colexp.op

		n_files=`ls ${GtWX_VAR_BLOCK_COL_EXP_DIR} | wc -l | awk {'print $1'}`
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Col expanded GtWX directory contains ${n_files} files"
		n_files=`ls ${PADDED_INV_XtWX_ROW_EXP_DIR} | wc -l | awk {'print $1'}`
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Row expanded padded inverse XtWX directory contains ${n_files} files"

		# Finally, multiply with Z.
		GtWXZ_FILE=${INTERM_DATA_DIR}/GtWXZ_PADDED_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.enc
		${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${GtWX_VAR_BLOCK_COL_EXP_DIR} ${PADDED_INV_XtWX_ROW_EXP_DIR} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${GtWXZ_FILE} >& ${INTERM_DATA_DIR}/GtWXZ_${cur_block_start_i}_${cur_block_end_i}_multip.op
	done

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading encrypted p-value statistic files to shared space.."
	#cp ${INTERM_DATA_DIR}/GtWG_VAR_BLOCK_*.enc ${SHARED_DIR}
	#cp ${INTERM_DATA_DIR}/P_CHISQ_STAT_VAR_BLOCK_*.enc ${SHARED_DIR}
	#cp ${INTERM_DATA_DIR}/GtWX_PADDED_VAR_BLOCK_*.enc ${SHARED_DIR}
	#cp ${INTERM_DATA_DIR}/GtWXZ_PADDED_VAR_BLOCK_*.enc ${SHARED_DIR}
	###############################################################################################################################################
	## This upload strips the file names to the basename and uses them to upload, should be fine to upload like this.
	ls ${INTERM_DATA_DIR}/GtWG_VAR_BLOCK_*.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	ls ${INTERM_DATA_DIR}/P_CHISQ_STAT_VAR_BLOCK_*.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	ls ${INTERM_DATA_DIR}/GtWX_PADDED_VAR_BLOCK_*.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	ls ${INTERM_DATA_DIR}/GtWXZ_PADDED_VAR_BLOCK_*.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	###############################################################################################################################################	

	exit
fi 

#-client_pool_pvalue_stats [Data config file] [Client index (starts at 0)] [Local data directory]
if [[ ${cmd_option} == "-client_pool_pvalue_stats" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

    data_config_file=$2
    cur_epoch=$3
    site_i_per_cmd=$4
    LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt
	
	# Take all other site's BETA's and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates in local features matrix file.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`

	row_col_iters=`seq 0 ${n_covars_min_one}`
	site_iters=`seq 0 ${n_site_min_one}`

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	# Download and pool G'(y-y0)
	PER_BLOCK_POOLED_Gt_YminY0_DIR=${INTERM_DATA_DIR}/per_block_pooled_Gt_YminY0_${site_i_per_cmd}
	rm -f -r ${PER_BLOCK_POOLED_Gt_YminY0_DIR}
	mkdir ${PER_BLOCK_POOLED_Gt_YminY0_DIR}
	PER_SITE_Gt_YminY0_dir=${INTERM_DATA_DIR}/PER_SITE_Gt_YminY0
	rm -f -r ${PER_SITE_Gt_YminY0_dir}
	mkdir ${PER_SITE_Gt_YminY0_dir}
	#cp ${SHARED_DIR}/P_CHISQ_STAT_VAR_BLOCK_*.enc ${PER_SITE_Gt_YminY0_dir}	

	PER_BLOCK_POOLED_GtWX_DIR=${INTERM_DATA_DIR}/per_block_pooled_GtWX_${site_i_per_cmd}
	rm -f -r ${PER_BLOCK_POOLED_GtWX_DIR}
	mkdir ${PER_BLOCK_POOLED_GtWX_DIR}
	PER_SITE_GtWX_dir=${INTERM_DATA_DIR}/PER_SITE_GtWX
	rm -f -r ${PER_SITE_GtWX_dir}
	mkdir ${PER_SITE_GtWX_dir}
	#cp ${SHARED_DIR}/GtWX_PADDED_VAR_BLOCK_*.enc ${PER_SITE_GtWX_dir}

	PER_BLOCK_POOLED_GtWXZ_DIR=${INTERM_DATA_DIR}/per_block_pooled_GtWXZ_${site_i_per_cmd}
    rm -f -r ${PER_BLOCK_POOLED_GtWXZ_DIR}
    mkdir ${PER_BLOCK_POOLED_GtWXZ_DIR}
    PER_SITE_GtWXZ_dir=${INTERM_DATA_DIR}/PER_SITE_GtWXZ
    rm -f -r ${PER_SITE_GtWXZ_dir}
    mkdir ${PER_SITE_GtWXZ_dir}
    #cp ${SHARED_DIR}/GtWXZ_PADDED_VAR_BLOCK_*.enc ${PER_SITE_GtWXZ_dir}

    PER_BLOCK_POOLED_GtWG_DIR=${INTERM_DATA_DIR}/per_block_pooled_GtWG_${site_i_per_cmd}
    rm -f -r ${PER_BLOCK_POOLED_GtWG_DIR}
    mkdir ${PER_BLOCK_POOLED_GtWG_DIR}
    PER_SITE_GtWG_dir=${INTERM_DATA_DIR}/PER_SITE_GtWG
    rm -f -r ${PER_SITE_GtWG_dir}
    mkdir ${PER_SITE_GtWG_dir}
    #cp ${SHARED_DIR}/GtWG_VAR_BLOCK_*.enc ${PER_SITE_GtWG_dir}

	# Download all files in all blocks.
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
	for cur_block_start_i in ${block_start_i_list[@]}
	do
		# Set the end of the current variant block.
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		###############################################################################################################################################
		# File downloading.
		# Wait for files to be ready: Get the file list first.
		rm -f ${INTERM_DATA_DIR}/Gt_YminY0_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/GtWX_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/GtWXZ_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/GtWG_TEMP_FILES.list	
		for site_i_f in ${site_iters[@]}
		do
			echo P_CHISQ_STAT_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/Gt_YminY0_TEMP_FILES.list	
			echo GtWX_PADDED_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/GtWX_TEMP_FILES.list	
			echo GtWXZ_PADDED_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.enc >> ${INTERM_DATA_DIR}/GtWXZ_TEMP_FILES.list	
			echo GtWG_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/GtWG_TEMP_FILES.list	
		done

		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/Gt_YminY0_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWX_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWXZ_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWG_TEMP_FILES.list
		cur_status=$?

		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/Gt_YminY0_TEMP_FILES.list ${PER_SITE_Gt_YminY0_dir}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWX_TEMP_FILES.list ${PER_SITE_GtWX_dir}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWXZ_TEMP_FILES.list ${PER_SITE_GtWXZ_dir}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/GtWG_TEMP_FILES.list ${PER_SITE_GtWG_dir}
		wait
		###############################################################################################################################################
	done

	# Processing 10,000,000 variants max.
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
	for cur_block_start_i in ${block_start_i_list[@]}
	do
		# Set the end of the current variant block.
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		###############################################
		ls ${PER_SITE_Gt_YminY0_dir}/*_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_*.enc > ${INTERM_DATA_DIR}/temp_stat_files_Gt_YminY0.list
		n_files=`wc -l ${INTERM_DATA_DIR}/temp_stat_files_Gt_YminY0.list | awk {'print $1'}`
		if [[ ${n_files} == 0 ]]
		then
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: No more files."
			break
		fi

		echo "Pooling ${n_files} encrypted Gt(Y-Y0) stats for block [${cur_block_start_i}-${cur_block_end_i}].."
		${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_stat_files_Gt_YminY0.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${PER_BLOCK_POOLED_Gt_YminY0_DIR}/P_CHISQ_STAT_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/Gt_YminY0_pooling.op

		###############################################
                ls ${PER_SITE_GtWX_dir}/*_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_*.enc > ${INTERM_DATA_DIR}/temp_stat_files_GtWX.list
                n_files=`wc -l ${INTERM_DATA_DIR}/temp_stat_files_GtWX.list | awk {'print $1'}`
                if [[ ${n_files} == 0 ]]
                then
                        date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
						echo "${date_time_str} @ $1: No more files."
                        break
                fi

                echo "Pooling ${n_files} encrypted GtWX stats for block [${cur_block_start_i}-${cur_block_end_i}].."
                ${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_stat_files_GtWX.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${PER_BLOCK_POOLED_GtWX_DIR}/GtWX_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/GtWX_pooling.op

		###############################################
		ls ${PER_SITE_GtWXZ_dir}/*_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_*.enc > ${INTERM_DATA_DIR}/temp_stat_files_GtWXZ.list
                n_files=`wc -l ${INTERM_DATA_DIR}/temp_stat_files_GtWXZ.list | awk {'print $1'}`
                if [[ ${n_files} == 0 ]]
                then
                        date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
						echo "${date_time_str} @ $1: No more files."
                        break
                fi

                echo "Pooling ${n_files} encrypted GtWXZ stats for block [${cur_block_start_i}-${cur_block_end_i}].."
                ${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_stat_files_GtWXZ.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${PER_BLOCK_POOLED_GtWXZ_DIR}/GtWXZ_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/GtWX_pooling.op

		###############################################
		ls ${PER_SITE_GtWG_dir}/*_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_*.enc > ${INTERM_DATA_DIR}/temp_stat_files_GtWG.list
                n_files=`wc -l ${INTERM_DATA_DIR}/temp_stat_files_GtWG.list | awk {'print $1'}`
                if [[ ${n_files} == 0 ]]
                then
                        date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
						echo "${date_time_str} @ $1: No more files."
                        break
                fi

                echo "Pooling ${n_files} encrypted GtWG stats for block [${cur_block_start_i}-${cur_block_end_i}].."
                ${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_stat_files_GtWG.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${PER_BLOCK_POOLED_GtWG_DIR}/GtWG_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/GtXG_pooling.op

		# Wait for the jobs to finish.
		wait

		###############################################
		# For the current block, calculate (Gt(Y-Y0)^2, GtWX Z, and P=GtWX Z XtWG (row-row multiplication)
		${COLLAGENE_SECURE_EXEC} -secure_row2row_inner_prod_continuous_encrypted_matrices ${PER_BLOCK_POOLED_Gt_YminY0_DIR}/P_CHISQ_STAT_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${PER_BLOCK_POOLED_Gt_YminY0_DIR}/P_CHISQ_STAT_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/S_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc

		${COLLAGENE_SECURE_EXEC} -secure_row2row_inner_prod_continuous_encrypted_matrices ${PER_BLOCK_POOLED_GtWX_DIR}/GtWX_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${PER_BLOCK_POOLED_GtWXZ_DIR}/GtWXZ_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/T2_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc

		${COLLAGENE_SECURE_EXEC} -secure_sub_cont_ct_matrices ${PER_BLOCK_POOLED_GtWG_DIR}/GtWG_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${INTERM_DATA_DIR}/T2_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${INTERM_DATA_DIR}/T_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc
		
		# Build the final chi-square statistics for this block: Generate random multiplier for S and T vectors.
		ST_NOISE_MATRIX_FILE=${INTERM_DATA_DIR}/ST_noise_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin
		rm -f ${ST_NOISE_MATRIX_FILE}
	    ${COLLAGENE_SECURE_EXEC} -generate_mult_full_noise_matrix ${VAR_BLOCK_SIZE} 1 ${ST_NOISE_MATRIX_FILE}

		NO_ST_NOISE=0
		if [[ ${NO_ST_NOISE} == 1  ]]
		then
			one_over_n_sites=`echo ${N_SITES} | awk {'print 1.0'}`
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Setting 1/n_sites=${one_over_n_sites}"
			COLLAGENE -generate_constant_matrix ${VAR_BLOCK_SIZE} 1 ${one_over_n_sites} ${ST_NOISE_MATRIX_FILE}
		fi

		${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${ST_NOISE_MATRIX_FILE} ${ST_NOISE_MATRIX_FILE}.txt

		${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${ST_NOISE_MATRIX_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${ST_NOISE_MATRIX_FILE}.enc >& ${INTERM_DATA_DIR}/ST_noise_encryption_${cur_block_start_i}_${cur_block_end_i}_enc.op
		
		# Add site-specific noise to the statistics upload them to shared space.
		# Multiply S & T with the ST noise vector.
		NOISY_S_STATS_FILE=${INTERM_DATA_DIR}/noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc
		${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${INTERM_DATA_DIR}/S_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${ST_NOISE_MATRIX_FILE}.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${NOISY_S_STATS_FILE}

		NOISY_T_STATS_FILE=${INTERM_DATA_DIR}/noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc
		${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${INTERM_DATA_DIR}/T_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc ${ST_NOISE_MATRIX_FILE}.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${NOISY_T_STATS_FILE}

		# Submit the noisy S and T stats to shared space.
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Uploading noisy S and T statistics for shared space.."
		#cp ${NOISY_S_STATS_FILE} ${NOISY_T_STATS_FILE} ${SHARED_DIR}
		echo ${NOISY_S_STATS_FILE} ${NOISY_T_STATS_FILE} > ${INTERM_DATA_DIR}/TEMP_FILES.list
		${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	done
fi

if [[ ${cmd_option} == "-pool_noisy_ST_stats" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

        data_config_file=$2
        cur_epoch=$3
        site_i_per_cmd=$4
        LOCAL_DATA_DIR=$5

        INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

        if [[ ! -d ${INTERM_DATA_DIR} ]]
        then
                mkdir ${INTERM_DATA_DIR}
        fi


	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling noisy ST statistics from all sites."

	# Pool the noisy statistics from all sites and save them locally.
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_NOISY_ST_DIR=${INTERM_DATA_DIR}/NOISY_ST_STATS
	rm -f -r ${PER_SITE_NOISY_ST_DIR}
	mkdir ${PER_SITE_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/noisy_S_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/noisy_T_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}

	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
	for cur_block_start_i in ${block_start_i_list[@]}
	do
		# Set the end of the current variant block.
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		###############################################################################################################################################
		# File downloading.
		# Wait for files to be ready: Get the file list first.
		rm -f ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
		for site_i_f in ${site_iters[@]}
		do
			echo noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
			echo noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
		done
		###############################################################################################################################################

		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list
		cur_status=$?

		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list ${PER_SITE_NOISY_ST_DIR}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list ${PER_SITE_NOISY_ST_DIR}
	done
	
	POOLED_NOISY_S_STAT_DIR=${INTERM_DATA_DIR}/POOLED_NOISY_S_STATS
	POOLED_NOISY_T_STAT_DIR=${INTERM_DATA_DIR}/POOLED_NOISY_T_STATS

	rm -f -r ${POOLED_NOISY_S_STAT_DIR}
	rm -f -r ${POOLED_NOISY_T_STAT_DIR}

	mkdir ${POOLED_NOISY_S_STAT_DIR}
	mkdir ${POOLED_NOISY_T_STAT_DIR}

    # Loop over all blocks.
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
    for cur_block_start_i in ${block_start_i_list[@]}
    do
        cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		rm -f ${INTERM_DATA_DIR}/cur_block_S_files.txt
		rm -f ${INTERM_DATA_DIR}/cur_block_T_files.txt
		for cur_site_i in ${site_iters[@]}
		do
                        cur_enc_S_file=${PER_SITE_NOISY_ST_DIR}/noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${cur_site_i}.bin.enc
                        cur_enc_T_file=${PER_SITE_NOISY_ST_DIR}/noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${cur_site_i}.bin.enc

                        if [[ ! -f ${cur_enc_S_file} ]]
                        then
                                echo "Could not find ${cur_enc_S_file}"
                                exit
                        fi

                        if [[ ! -f ${cur_enc_T_file} ]]
                        then
                                echo "Could not find ${cur_enc_T_file}"
                                exit
                        fi

			echo ${cur_enc_S_file} >> ${INTERM_DATA_DIR}/cur_block_S_files.txt
			echo ${cur_enc_T_file} >> ${INTERM_DATA_DIR}/cur_block_T_files.txt
		done

		n_S_files=`wc -l ${INTERM_DATA_DIR}/cur_block_S_files.txt | awk {'print $1'}`
		n_T_files=`wc -l ${INTERM_DATA_DIR}/cur_block_T_files.txt | awk {'print $1'}`

		echo "Block ${cur_block_start_i}-${cur_block_end_i}::Found ${n_S_files} S and ${n_T_files} T files, pooling.."

		# Pool.
		${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/cur_block_S_files.txt ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_NOISY_S_STAT_DIR}/pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/noisy_S_pooling.txt &

		${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/cur_block_T_files.txt ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${POOLED_NOISY_T_STAT_DIR}/pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc >& ${INTERM_DATA_DIR}/noisy_T_pooling.txt &

		# Wait for the current poolings to end.
		wait
	done
	
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading pooled noisy S/T statististics to shared space."
	#cp ${POOLED_NOISY_S_STAT_DIR}/pooled_noisy_S_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc ${SHARED_DIR}
    #cp ${POOLED_NOISY_T_STAT_DIR}/pooled_noisy_T_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc ${SHARED_DIR}
	ls ${POOLED_NOISY_S_STAT_DIR}/pooled_noisy_S_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	ls ${POOLED_NOISY_T_STAT_DIR}/pooled_noisy_T_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list

	exit
fi

# -client_collaborative_decrypt_noisy_ST_stats: Collaboratively decrypt noisy ST statistics.
if [[ ${cmd_option} == "-client_collaborative_decrypt_noisy_ST_stats" ]]
then
	if [[ $# != 6 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Cur site secret key] [Epoch index] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	private_key_share_file=$3
	cur_epoch=$4
	site_i_per_cmd=$5
	LOCAL_DATA_DIR=$6

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi

	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Partially decrypting noisy S/T statistics on site ${site_i_per_cmd}"

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`
		
	PER_SITE_NOISY_ST_DIR=${INTERM_DATA_DIR}/NOISY_ST_STATS
	rm -f -r ${PER_SITE_NOISY_ST_DIR}
	mkdir ${PER_SITE_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/pooled_noisy_S_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/pooled_noisy_T_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
	for cur_block_start_i in ${block_start_i_list[@]}
	do
		# Set the end of the current variant block.
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		###############################################################################################################################################
		# File downloading.
		# Wait for files to be ready: Get the file list first.
		rm -f ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
		for site_i_f in ${site_iters[@]}
		do
			echo pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
			echo pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
		done
		###############################################################################################################################################

		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list
		cur_status=$?

		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list ${PER_SITE_NOISY_ST_DIR}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list ${PER_SITE_NOISY_ST_DIR}
	done
 
	is_site0=0
	if [[ ${site_i_per_cmd} == 0 ]]
	then
		is_site0=1
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: is_site0=${is_site0}"

	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	# Loop over all blocks.
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
    for cur_block_start_i in ${block_start_i_list[@]}
    do
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`
		
		for cur_dec_site_i in ${site_iters[@]}
		do
			cur_enc_S_file=${PER_SITE_NOISY_ST_DIR}/pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${cur_dec_site_i}.bin.enc
			cur_enc_T_file=${PER_SITE_NOISY_ST_DIR}/pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${cur_dec_site_i}.bin.enc

			if [[ ! -f ${cur_enc_S_file} ]]
			then
				echo "Could not find ${cur_enc_S_file}"
				exit
			fi

			if [[ ! -f ${cur_enc_T_file} ]]
			then
				echo "Could not find ${cur_enc_T_file}"
				exit
			fi
				
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Partially decrypting ${cur_S_file} using Site-${site_i_per_cmd}'s key."

			cur_enc_S_file_partdec_file=${cur_enc_S_file}_partdec_by_${site_i_per_cmd}.partdec
			${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_enc_S_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_enc_S_file_partdec_file}

			cur_enc_S_file_partdec_enc_file=${cur_enc_S_file}_partdec_by_${site_i_per_cmd}.partdec.enc
			$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_enc_S_file_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_enc_S_file_partdec_enc_file}
				
			date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
			echo "${date_time_str} @ $1: Partially decrypting ${cur_T_file} using Site-${site_i_per_cmd}'s key."

			cur_enc_T_file_partdec_file=${cur_enc_T_file}_partdec_by_${site_i_per_cmd}.partdec
			${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_enc_T_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_enc_T_file_partdec_file}

			cur_enc_T_file_partdec_enc_file=${cur_enc_T_file}_partdec_by_${site_i_per_cmd}.partdec.enc
			$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_enc_T_file_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_enc_T_file_partdec_enc_file}
		done
	done 

	# Copy the partially decrypted noisy beta.
	#cp ${PER_SITE_NOISY_ST_DIR}/*.partdec ${SHARED_DIR}
	#ls ${PER_SITE_NOISY_ST_DIR}/*.partdec > ${INTERM_DATA_DIR}/TEMP_FILES.list
	ls ${PER_SITE_NOISY_ST_DIR}/*.partdec.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list

	exit 0
fi


if [[ ${cmd_option} == "-client_pool_partially_decrypted_noisy_ST_stats" ]]
then
    if [[ $# != 5 ]]
    then
            echo "USAGE: $0 $1 [Data config file] [Epoch index] [Client index (starts at 0)] [Local data directory]"
            exit
    fi

    data_config_file=$2
    cur_epoch=$3
    site_i_per_cmd=$4
    LOCAL_DATA_DIR=$5

    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

    if [[ ! -d ${INTERM_DATA_DIR} ]]
    then
            mkdir ${INTERM_DATA_DIR}
    fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	if [[ "${Z_SCALER}" == "" ]]
	then
		echo "Need to set Z_SCALER"
		exit
	fi

	LOCAL_FEAT_FILE=${LOCAL_DATA_DIR}/feat_matrix.txt
	LOCAL_PHENO_FILE=${LOCAL_DATA_DIR}/phenotypes.txt

	# Take all other site's XtWX and decrypt them.
	n_covars_min_one=`awk {'print NF-2'} ${LOCAL_FEAT_FILE} | head -n 1`
	n_covars=`awk {'print NF-1'} ${LOCAL_FEAT_FILE} | head -n 1`
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_covars_min_one} covariates on local features matrix file.."
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR=${INTERM_DATA_DIR}/PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
	mkdir ${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
	PER_SITE_PARTDEC_NOISY_ST_STATS_DIR=${INTERM_DATA_DIR}/PER_SITE_PARTDEC_NOISY_ST_STATS_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}
	mkdir ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}

	#cp -r ${SHARED_DIR}/pooled_noisy_S_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}
	#cp -r ${SHARED_DIR}/pooled_noisy_T_VAR_BLOCK_*_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}
	
	###############################################################################################################################################
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
	for cur_block_start_i in ${block_start_i_list[@]}
	do
		# Set the end of the current variant block.
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		###############################################################################################################################################
		# File downloading.
		# Wait for files to be ready: Get the file list first.
		rm -f ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
		rm -f ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
		for site_i_f in ${site_iters[@]}
		do
			echo pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/S_TEMP_FILES.list
			echo pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/T_TEMP_FILES.list
		done
		###############################################################################################################################################

		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list
		cur_status=$?
		${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list
		cur_status=$?

		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list ${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
		${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list ${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}

		for site_i_f in ${site_iters[@]}
		do
			cur_enc_partdec_file=${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc
			cur_partdec_file=${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec
			$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}

			cur_enc_partdec_file=${PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc
			cur_partdec_file=${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec
			$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}
		done
	done
	###############################################################################################################################################

	POOLED_DEC_NOISY_S_MATRIX_DIR=${INTERM_DATA_DIR}/FULL_DEC_POOLED_NOISY_S_MATRICES
	POOLED_DEC_NOISY_T_MATRIX_DIR=${INTERM_DATA_DIR}/FULL_DEC_POOLED_NOISY_T_MATRICES
	rm -f -r ${POOLED_DEC_NOISY_S_MATRIX_DIR}
	rm -f -r ${POOLED_DEC_NOISY_T_MATRIX_DIR}
	mkdir ${POOLED_DEC_NOISY_S_MATRIX_DIR} ${POOLED_DEC_NOISY_T_MATRIX_DIR}

	# Loop over all blocks.
	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`
    for cur_block_start_i in ${block_start_i_list[@]}
    do
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		ls ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_S_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec > ${INTERM_DATA_DIR}/temp_block_S_partdecs.list
		ls ${PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/pooled_noisy_T_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec > ${INTERM_DATA_DIR}/temp_block_T_partdecs.list

		n_S_partdecs=`wc -l ${INTERM_DATA_DIR}/temp_block_S_partdecs.list | awk '{print $1}'`
		n_T_partdecs=`wc -l ${INTERM_DATA_DIR}/temp_block_T_partdecs.list | awk '{print $1}'`

		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Pooling ${n_S_partdecs} partdecs for S and ${n_T_partdecs} partdecs for T"

		CUR_BLOCK_FULLDEC_POOLED_NOISY_S=${POOLED_DEC_NOISY_S_MATRIX_DIR}/full_dec_S_BLOCK_${cur_block_start_i}_${cur_block_end_i}.bin
		${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_block_S_partdecs.list ${TEXT_PARAMS_PATH} ${CUR_BLOCK_FULLDEC_POOLED_NOISY_S}
		${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${CUR_BLOCK_FULLDEC_POOLED_NOISY_S} ${CUR_BLOCK_FULLDEC_POOLED_NOISY_S}.txt

		CUR_BLOCK_FULLDEC_POOLED_NOISY_T=${POOLED_DEC_NOISY_T_MATRIX_DIR}/full_dec_T_BLOCK_${cur_block_start_i}_${cur_block_end_i}.bin
		${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_block_T_partdecs.list ${TEXT_PARAMS_PATH} ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}

		Z_UNSCALER=`echo ${Z_SCALER} | awk '{print 1.0/$1}'`
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Z_UNSCALER: ${Z_UNSCALER}"
		${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T} ${Z_UNSCALER} ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}_unscaled.bin

		${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}_unscaled.bin ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}.txt

		# Finally, calculate and save the p-values.
		cur_block_var_ids=${INTERM_DATA_DIR}/VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_IDS_CLIENT_${site_i_per_cmd}.list
		paste ${CUR_BLOCK_FULLDEC_POOLED_NOISY_S}.txt ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}.txt | awk 'BEGIN{FS="\t"}{print $2"\t"$4}' | ${COLLAGENE_SECURE_EXEC} -assign_chisqr_pvals_per_ST_stats stdin ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.txt
		paste ${cur_block_var_ids} ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.txt > ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.txt_w_var_ids.txt

		# Use R's pchisq function.
		paste ${CUR_BLOCK_FULLDEC_POOLED_NOISY_S}.txt ${CUR_BLOCK_FULLDEC_POOLED_NOISY_T}.txt | awk 'BEGIN{FS="\t"}{print $2"\t"$4}' > ${INTERM_DATA_DIR}/final_ST_stats_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.txt
		Rscript GET_CHISQ_PVAL.R ${INTERM_DATA_DIR}/final_ST_stats_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}.txt ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}_PER_R.txt
		paste ${cur_block_var_ids} ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}_PER_R.txt > ${INTERM_DATA_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}_PER_R_w_var_ids.txt
	done

	# Pool the final p-values to be used as the final list.
	$0 -pool_final_plaintext_p_values data_config.params ${site_i_per_cmd} ${LOCAL_DATA_DIR}

    exit
fi

if [[ ${cmd_option} == "-pool_final_plaintext_p_values" ]] 
then
    if [[ $# != 4 ]]
    then
		echo "USAGE: $0 $1 [Data config file] [Client index (starts at 0)] [Local data directory]"
		exit 1
    fi

	data_config_file=$2
    site_i_per_cmd=$3
    LOCAL_DATA_DIR=$4
	
    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE
	P_VAL_DIR=${INTERM_DATA_DIR}

    if [[ ! -d ${P_VAL_DIR} ]]
    then
		echo "Could not find the final p-values directory ${P_VAL_DIR}"
		exit 1
    fi

	#requested_last_block_start_i=`echo -e ${VAR_BLOCK_SIZE}"\t"${N_VARS_2_SCORE} | awk {'win_start_i=($2-$1);if($2<$1){win_start_i=0};print win_start_i'}`
	requested_last_block_start_i=`$0 -get_last_win_start_i $2`
    block_start_i_list=`seq 0 ${VAR_BLOCK_SIZE} ${requested_last_block_start_i}`

	FINAL_PVALS_FILE=FINAL_P_VALUES_ALL_VAR_BLOCKS_CLIENT_${site_i_per_cmd}_PER_R_w_var_ids.txt
	rm -f ${FINAL_PVALS_FILE}
	for cur_block_start_i in ${block_start_i_list}
	do
		cur_block_end_i=`echo ${cur_block_start_i} | awk -v VAR_BLOCK_SIZE=${VAR_BLOCK_SIZE} {'print $1+VAR_BLOCK_SIZE-1'}`

		CUR_BLOCK_PVAL_FILE=${P_VAL_DIR}/FINAL_P_VALUES_VAR_BLOCK_${cur_block_start_i}_${cur_block_end_i}_CLIENT_${site_i_per_cmd}_PER_R_w_var_ids.txt

		if [[ ! -f ${CUR_BLOCK_PVAL_FILE} ]]
		then
			echo "Could not find ${CUR_BLOCK_PVAL_FILE}"
			exit 1
		fi

		echo "Concatting ${CUR_BLOCK_PVAL_FILE}"
		cat ${CUR_BLOCK_PVAL_FILE} >> ${FINAL_PVALS_FILE}
	done

	exit 0
fi

######################################################################################################################################################
######################################################################################################################################################

if [[ "${cmd_option}" == "-validate_meta_analysis_input_data" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Local data directory]"
		exit 1
	fi

	data_config_file=$2
	LOCAL_DATA_DIR=$3
	
	rm -f INPUT.ERROR
	
	if [[ ! -d ${LOCAL_DATA_DIR} ]]
	then
		echo "Could not find local data directory @ ${LOCAL_DATA_DIR}"
		echo "ERROR" > INPUT.ERROR
		exit
	fi

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
			mkdir ${INTERM_DATA_DIR}

			if [[ ! -d ${INTERM_DATA_DIR} ]]
			then
				echo "Could not create the intermediate directory @ \"${INTERM_DATA_DIR}\""
				exit 1
			fi
	fi

	LOCAL_GMMAT_RESULTS_FILE=${LOCAL_DATA_DIR}/GMMAT_RESULTS.txt
	if [[ ! -f ${LOCAL_GMMAT_RESULTS_FILE} ]]
	then
		echo "Could not find features file @ \"${LOCAL_FEAT_FILE}\""
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PUBLIC_KEY_FILE} ]]
	then
	echo "Could not find public keys file @ ${PUBLIC_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${RELIN_KEY_FILE} ]]
	then
		echo "Could not find the relin. keys file @ ${RELIN_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${GALOIS_KEY_FILE} ]]
	then
		echo "Could not find the Galois keys file @ ${GALOIS_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PRIVATE_KEY_FILE} ]]
	then
		echo "Could not find the private key file @ ${PRIVATE_KEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${PARTDEC_SYMKEY_FILE} ]]
	then
		echo "Could not find the symmetric encryption key file @ ${PARTDEC_SYMKEY_FILE}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find the security parameters file @ ${TEXT_PARAMS_PATH}"
		echo "ERROR" > INPUT.ERROR
		exit 1
	fi

	# Test a file transfer.
	echo "Testing a file transfer."
	${FILE_IO_UTILS_SCRIPT} -test_file_IO ${data_config_file}
	if [[ $? != 0 ]]
	then
		echo "File transfer test failed."
		echo "ERROR" > INPUT.ERROR
	fi

	# Exit with no errors.
	exit 0
fi # -validate_meta_analysis_input_data option.

#-meta_encrypt_save_ST_stats [Data config file] [Output directory]: Encrypt and upload S and T statistics.
if [[ ${cmd_option} == "-meta_encrypt_save_ST_stats" ]]
then
    if [[ $# != 4 ]]
    then
            echo "USAGE: $0 $1 [Data config file] [Client index (starts at 0)] [Local data directory]"
            exit
    fi

    data_config_file=$2
    site_i_per_cmd=$3
    LOCAL_DATA_DIR=$4

    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

    if [[ ! -d ${INTERM_DATA_DIR} ]]
    then
            mkdir ${INTERM_DATA_DIR}
    fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	LOCAL_GMMAT_RESULTS_FP=${LOCAL_DATA_DIR}/GMMAT_RESULTS.txt

	if [[ ! -f ${LOCAL_GMMAT_RESULTS_FP} ]]
	then
		echo "Could not find GMMAT results file @ \"${LOCAL_GMMAT_RESULTS_FP}\""
		exit 1
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Extracting S/T statistics from GMMAT file ${LOCAL_GMMAT_RESULTS_FP}, encrypting them.."

	# Basically, parse the SCORE and VAR columns of the GMMAT results.
	# Separate the S and T statistics from GMMAT results, save them as two vectors.
	META_GMMAT_GtYminY0_FILE=${INTERM_DATA_DIR}/META_GMMAT_Gt_YminY0_CLIENT_${site_i_per_cmd}.bin
	META_GMMAT_T_FILE=${INTERM_DATA_DIR}/META_GMMAT_T_CLIENT_${site_i_per_cmd}.bin
	cat ${LOCAL_GMMAT_RESULTS_FP} | cut -f6,6 | ${COLLAGENE_SECURE_EXEC} -save_matrix_plain_2_bin stdin 0 1 ${META_GMMAT_GtYminY0_FILE}
	cat ${LOCAL_GMMAT_RESULTS_FP} | cut -f7,7 | ${COLLAGENE_SECURE_EXEC} -save_matrix_plain_2_bin stdin 0 1 ${META_GMMAT_T_FILE}

	# Encrypt S and T matrices.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Encrypting S/T statistics."
	${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${META_GMMAT_GtYminY0_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${META_GMMAT_GtYminY0_FILE}.enc
	${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${META_GMMAT_T_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${META_GMMAT_T_FILE}.enc
	
	# Upload the encrypted stats and the noise to the shared space.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading encrypted S/T statistics to shared space."
	echo ${META_GMMAT_GtYminY0_FILE}.enc ${META_GMMAT_T_FILE}.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
fi

# Download all of the noise and ST stats, pool S over all sites, square it, add this site's noise to pooled S and T stats.
if [[ ${cmd_option} == "-meta_pool_per_site_ST_stats_add_noise" ]]
then
    if [[ $# != 4 ]]
    then
            echo "USAGE: $0 $1 [Data config file] [Client index (starts at 0)] [Local data directory]"
            exit
    fi

    data_config_file=$2
    site_i_per_cmd=$3
    LOCAL_DATA_DIR=$4

    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

    if [[ ! -d ${INTERM_DATA_DIR} ]]
    then
            mkdir ${INTERM_DATA_DIR}
    fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi
	
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling encrypted ST statistics and adding noise.."

	# Create the directory.
	PER_SITE_GMMAT_ST_STATS_DIR=${INTERM_DATA_DIR}/PER_SITE_META_ST_STATS
	rm -f -r ${PER_SITE_GMMAT_ST_STATS_DIR}
	mkdir ${PER_SITE_GMMAT_ST_STATS_DIR}

	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo META_GMMAT_Gt_YminY0_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
		echo META_GMMAT_T_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_GMMAT_ST_STATS_DIR}
	#cp -r ${SHARED_DIR}/col_exp_mult_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_NOISE_EXP_DIR}
    #cp -r ${SHARED_DIR}/col_exp_mult_padded_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_PADDED_NOISE_EXP_DIR}
	###############################################################################################################################################

	# Pool the statistics.
	META_POOLED_GtYminY0_STATS_FILE=${INTERM_DATA_DIR}/META_POOLED_GtYminY0_STATS.bin.enc
	ls ${PER_SITE_GMMAT_ST_STATS_DIR}/META_GMMAT_Gt_YminY0_CLIENT_*.bin.enc > ${INTERM_DATA_DIR}/temp_meta_Gt_YminY0_files.list
	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_meta_Gt_YminY0_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_POOLED_GtYminY0_STATS_FILE}

	META_POOLED_S_STATS_FILE=${INTERM_DATA_DIR}/META_POOLED_S_STATS.bin.enc
	${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${META_POOLED_GtYminY0_STATS_FILE} ${META_POOLED_GtYminY0_STATS_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_POOLED_S_STATS_FILE}

	META_POOLED_T_STATS_FILE=${INTERM_DATA_DIR}/META_POOLED_T_STATS.bin.enc
	ls ${PER_SITE_GMMAT_ST_STATS_DIR}/META_GMMAT_T_CLIENT_*.bin.enc > ${INTERM_DATA_DIR}/temp_meta_T_files.list
	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_meta_T_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_POOLED_T_STATS_FILE}

	# Add the noise to S and T.
	# Build the final chi-square statistics for this block: Generate random multiplier for S and T vectors.
	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${META_POOLED_T_STATS_FILE} ${INTERM_DATA_DIR}/meta_T_dims.txt
	n_vars=`cut -f1,1 ${INTERM_DATA_DIR}/meta_T_dims.txt`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Found ${n_vars} variants in GMMAT results file."

	META_ST_NOISE_MATRIX_FILE=${INTERM_DATA_DIR}/META_ST_noise.bin
	rm -f ${META_ST_NOISE_MATRIX_FILE}
	${COLLAGENE_SECURE_EXEC} -generate_mult_full_noise_matrix ${n_vars} 1 ${META_ST_NOISE_MATRIX_FILE}

	NO_ST_NOISE=0
	if [[ ${NO_ST_NOISE} == 1  ]]
	then
		one_over_n_sites=`echo ${N_SITES} | awk {'print 1.0'}`
		echo "Setting 1/n_sites=${one_over_n_sites}"
		COLLAGENE -generate_constant_matrix ${VAR_BLOCK_SIZE} ${n_vars} 1 ${META_ST_NOISE_MATRIX_FILE}
	fi

	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${META_ST_NOISE_MATRIX_FILE} ${META_ST_NOISE_MATRIX_FILE}.txt

	${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${META_ST_NOISE_MATRIX_FILE} ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${META_ST_NOISE_MATRIX_FILE}.enc >& ${INTERM_DATA_DIR}/meta_ST_noise_encryption.op
		
	# Add site-specific noise to the statistics upload them to shared space.
	# Multiply S & T with the ST noise vector.
	META_NOISY_S_STATS_FILE=${INTERM_DATA_DIR}/META_noisy_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc
	${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${META_POOLED_S_STATS_FILE} ${META_ST_NOISE_MATRIX_FILE}.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_NOISY_S_STATS_FILE}

	META_NOISY_T_STATS_FILE=${INTERM_DATA_DIR}/META_noisy_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc
	${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${META_POOLED_T_STATS_FILE} ${META_ST_NOISE_MATRIX_FILE}.enc ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_NOISY_T_STATS_FILE}

	############################################################################################################################################
	# Submit the noisy S and T stats to shared space.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading the per-site-noisy pooled S/T statistics to shared space."
	echo ${META_NOISY_S_STATS_FILE} ${META_NOISY_T_STATS_FILE} > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
fi

if [[ ${cmd_option} == "-meta_pool_noisy_ST_stats" ]]
then
    if [[ $# != 4 ]]
    then
            echo "USAGE: $0 $1 [Data config file] [Client index (starts at 0)] [Local data directory]"
            exit
    fi

    data_config_file=$2
    site_i_per_cmd=$3
    LOCAL_DATA_DIR=$4

    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

    if [[ ! -d ${INTERM_DATA_DIR} ]]
    then
            mkdir ${INTERM_DATA_DIR}
    fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling noisy S/T statistics from all sites."

	# Create the directory.
	PER_SITE_META_NOISY_ST_STATS_DIR=${INTERM_DATA_DIR}/PER_SITE_META_NOISY_ST_STATS
	rm -f -r ${PER_SITE_META_NOISY_ST_STATS_DIR}
	mkdir ${PER_SITE_META_NOISY_ST_STATS_DIR}

	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/TEMP_FILES.list
	for site_i_f in ${site_iters[@]}
	do
		echo META_noisy_S_STATS_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
		echo META_noisy_T_STATS_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/TEMP_FILES.list
	done

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list ${PER_SITE_META_NOISY_ST_STATS_DIR}
	#cp -r ${SHARED_DIR}/col_exp_mult_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_NOISE_EXP_DIR}
    #cp -r ${SHARED_DIR}/col_exp_mult_padded_noise_matrix_iter_* ${SHARED_DIR}/row_exp_mult_padded_noise_matrix_iter_${cur_epoch}_* ${PER_SITE_PADDED_NOISE_EXP_DIR}
	###############################################################################################################################################

	ls ${PER_SITE_META_NOISY_ST_STATS_DIR}/META_noisy_S_STATS_CLIENT_*.bin.enc > ${INTERM_DATA_DIR}/temp_noisy_S_files.list
	META_POOLED_NOISY_S_STATS_FILE=${INTERM_DATA_DIR}/META_POOLED_NOISY_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc
	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_noisy_S_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_POOLED_NOISY_S_STATS_FILE}

	ls ${PER_SITE_META_NOISY_ST_STATS_DIR}/META_noisy_T_STATS_CLIENT_*.bin.enc > ${INTERM_DATA_DIR}/temp_noisy_T_files.list
	META_POOLED_NOISY_T_STATS_FILE=${INTERM_DATA_DIR}/META_POOLED_NOISY_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc
	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${INTERM_DATA_DIR}/temp_noisy_T_files.list ${TEXT_PARAMS_PATH} ${PUBLIC_KEY_FILE} ${RELIN_KEY_FILE} ${GALOIS_KEY_FILE} ${PRIVATE_KEY_FILE} ${META_POOLED_NOISY_T_STATS_FILE}

	############################################################################################################################################
	# Submit the noisy S and T stats to shared space.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading pooled all-site-noisy S/T statistics to shared space."
	echo ${META_POOLED_NOISY_S_STATS_FILE} ${META_POOLED_NOISY_T_STATS_FILE} > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list
fi


if [[ ${cmd_option} == "-meta_collaborative_decrypt_pooled_noisy_ST_stats" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [Cur site secret key] [Client index (starts at 0)] [Local data directory]"
		exit
	fi

	data_config_file=$2
	private_key_share_file=$3
	site_i_per_cmd=$4
	LOCAL_DATA_DIR=$5

	INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

	if [[ ! -d ${INTERM_DATA_DIR} ]]
	then
		mkdir ${INTERM_DATA_DIR}
	fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Partial decrypting and symmetric encrypting all-site noisy S/T statistics."

	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`
		
	META_PER_SITE_POOLED_NOISY_ST_DIR=${INTERM_DATA_DIR}/META_PER_SITE_POOLED_NOISY_ST_STATS
	rm -f -r ${META_PER_SITE_POOLED_NOISY_ST_DIR}
	mkdir ${META_PER_SITE_POOLED_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/pooled_noisy_S_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}
	#cp -r ${SHARED_DIR}/pooled_noisy_T_VAR_BLOCK_*_CLIENT_*.bin.enc ${PER_SITE_NOISY_ST_DIR}
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
	rm -f ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
	for site_i_f in ${site_iters[@]}
	do
		echo META_POOLED_NOISY_S_STATS_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/S_TEMP_FILES.list
		echo META_POOLED_NOISY_T_STATS_CLIENT_${site_i_f}.bin.enc >> ${INTERM_DATA_DIR}/T_TEMP_FILES.list
	done
	###############################################################################################################################################

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list
	cur_status=$?
	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list ${META_PER_SITE_POOLED_NOISY_ST_DIR}
	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list ${META_PER_SITE_POOLED_NOISY_ST_DIR}
 
	is_site0=0
	if [[ ${site_i_per_cmd} == 0 ]]
	then
		is_site0=1
	fi

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: is_site0=${is_site0}"

	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	# Loop over all sites.
	for cur_dec_site_i in ${site_iters[@]}
	do
		cur_enc_S_file=${META_PER_SITE_POOLED_NOISY_ST_DIR}/META_POOLED_NOISY_S_STATS_CLIENT_${cur_dec_site_i}.bin.enc
		cur_enc_T_file=${META_PER_SITE_POOLED_NOISY_ST_DIR}/META_POOLED_NOISY_T_STATS_CLIENT_${cur_dec_site_i}.bin.enc

		if [[ ! -f ${cur_enc_S_file} ]]
		then
			echo "Could not find ${cur_enc_S_file}"
			exit
		fi

		if [[ ! -f ${cur_enc_T_file} ]]
		then
			echo "Could not find ${cur_enc_T_file}"
			exit
		fi
				
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Partially decrypting and symmetric encrypting ${cur_S_file} using Site-${site_i_per_cmd}'s key."
		${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_enc_S_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_enc_S_file}_partdec_by_${site_i_per_cmd}.partdec

		$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_enc_S_file}_partdec_by_${site_i_per_cmd}.partdec ${PARTDEC_SYMKEY_FILE} ${cur_enc_S_file}_partdec_by_${site_i_per_cmd}.partdec.enc
				
		date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
		echo "${date_time_str} @ $1: Partially decrypting and symmetric encrypting ${cur_T_file} using Site-${site_i_per_cmd}'s key."
		${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${cur_enc_T_file} ${is_site0} ${TEXT_PARAMS_PATH} ${private_key_share_file} 0 ${cur_enc_T_file}_partdec_by_${site_i_per_cmd}.partdec

		$0 -symmetric_encrypt_partdec_data ${data_config_file} ${cur_enc_T_file}_partdec_by_${site_i_per_cmd}.partdec ${PARTDEC_SYMKEY_FILE} ${cur_enc_T_file}_partdec_by_${site_i_per_cmd}.partdec.enc
	done

	# Copy the partially decrypted noisy beta.
	#cp ${PER_SITE_NOISY_ST_DIR}/*.partdec ${SHARED_DIR}
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Uploading partially decrypted symmetric encrypted all-site-noisy S/T statistics."
	ls ${META_PER_SITE_POOLED_NOISY_ST_DIR}/*.partdec.enc > ${INTERM_DATA_DIR}/TEMP_FILES.list
	${FILE_IO_UTILS_SCRIPT} -upload_files_to_shared ${data_config_file} ${INTERM_DATA_DIR}/TEMP_FILES.list

	exit 0
fi

if [[ ${cmd_option} == "-meta_pool_partially_decrypted_pooled_noisy_ST_stats" ]]
then
    if [[ $# != 4 ]]
    then
            echo "USAGE: $0 $1 [Data config file] [Client index (starts at 0)] [Local data directory]"
            exit
    fi

    data_config_file=$2
    site_i_per_cmd=$3
    LOCAL_DATA_DIR=$4

    INTERM_DATA_DIR=${LOCAL_DATA_DIR}/INTERMEDIATE

    if [[ ! -d ${INTERM_DATA_DIR} ]]
    then
            mkdir ${INTERM_DATA_DIR}
    fi
	
	if [[ ! -f ${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find \"${TEXT_PARAMS_PATH}\""
		exit
	fi

	# Take all other site's XtWX and decrypt them.
	n_site_min_one=`echo ${N_SITES} | awk '{print $1-1}'`
	site_iters=`seq 0 ${n_site_min_one}`

	META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR=${INTERM_DATA_DIR}/META_PER_SITE_ENC_PARTDEC_POOLED_NOISY_ST_STATS_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
	mkdir ${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling partially decrypted all-site noisy S/T statistics."
	
	###############################################################################################################################################
	# File downloading.
	# Wait for files to be ready: Get the file list first.
	rm -f ${INTERM_DATA_DIR}/T_TEMP_FILES.list	
	rm -f ${INTERM_DATA_DIR}/S_TEMP_FILES.list	
	for site_i_f in ${site_iters[@]}
	do
		echo META_POOLED_NOISY_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/S_TEMP_FILES.list
		echo META_POOLED_NOISY_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc >> ${INTERM_DATA_DIR}/T_TEMP_FILES.list
	done
	###############################################################################################################################################

	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list
	cur_status=$?
	${FILE_IO_UTILS_SCRIPT} -wait_for_files_in_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list
	cur_status=$?

	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/S_TEMP_FILES.list ${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
	${FILE_IO_UTILS_SCRIPT} -download_files_from_shared ${data_config_file} ${INTERM_DATA_DIR}/T_TEMP_FILES.list ${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}
	###############################################################################################################################################

	# Decrypt the files.
	META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR=${INTERM_DATA_DIR}/META_PER_SITE_PARTDEC_POOLED_NOISY_ST_STATS_CLIENT_${site_i_per_cmd}_ITER_${cur_epoch}
	rm -f -r ${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}
	mkdir ${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}

	# Go over each site and decrypt the files.
	for site_i_f in ${site_iters[@]}
	do
		cur_enc_partdec_file=${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc
		cur_partdec_file=${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec
		$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}

		cur_enc_partdec_file=${META_PER_SITE_ENC_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec.enc
		cur_partdec_file=${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_${site_i_f}.partdec
		$0 -symmetric_decrypt_partdec_data ${data_config_file} ${cur_enc_partdec_file} ${PARTDEC_SYMKEY_FILE} ${cur_partdec_file}
	done

	META_POOLED_DEC_NOISY_S_MATRIX_DIR=${INTERM_DATA_DIR}/META_FULL_DEC_POOLED_NOISY_S_MATRICES
	META_POOLED_DEC_NOISY_T_MATRIX_DIR=${INTERM_DATA_DIR}/META_FULL_DEC_POOLED_NOISY_T_MATRICES
	rm -f -r ${META_POOLED_DEC_NOISY_S_MATRIX_DIR}
	rm -f -r ${META_POOLED_DEC_NOISY_T_MATRIX_DIR}
	mkdir ${META_POOLED_DEC_NOISY_S_MATRIX_DIR} ${META_POOLED_DEC_NOISY_T_MATRIX_DIR}

	# Get the file name lists.
	ls ${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_S_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec > ${INTERM_DATA_DIR}/temp_block_S_partdecs.list
	ls ${META_PER_SITE_PARTDEC_NOISY_ST_STATS_DIR}/META_POOLED_NOISY_T_STATS_CLIENT_${site_i_per_cmd}.bin.enc_partdec_by_*.partdec > ${INTERM_DATA_DIR}/temp_block_T_partdecs.list

	n_S_partdecs=`wc -l ${INTERM_DATA_DIR}/temp_block_S_partdecs.list | awk '{print $1}'`
	n_T_partdecs=`wc -l ${INTERM_DATA_DIR}/temp_block_T_partdecs.list | awk '{print $1}'`

	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Pooling ${n_S_partdecs} partdecs for S and ${n_T_partdecs} partdecs for T"

	META_FULLDEC_POOLED_NOISY_S=${META_POOLED_DEC_NOISY_S_MATRIX_DIR}/meta_full_dec_S_BLOCK_CLIENT_${site_i_per_cmd}.bin
	${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_block_S_partdecs.list ${TEXT_PARAMS_PATH} ${META_FULLDEC_POOLED_NOISY_S}
	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${META_FULLDEC_POOLED_NOISY_S} ${META_FULLDEC_POOLED_NOISY_S}.txt

	META_FULLDEC_POOLED_NOISY_T=${META_POOLED_DEC_NOISY_T_MATRIX_DIR}/meta_full_dec_T_BLOCK_CLIENT_${site_i_per_cmd}.bin
	${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${INTERM_DATA_DIR}/temp_block_T_partdecs.list ${TEXT_PARAMS_PATH} ${META_FULLDEC_POOLED_NOISY_T}
	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${META_FULLDEC_POOLED_NOISY_T} ${META_FULLDEC_POOLED_NOISY_T}.txt

	# Finally, calculate and save the p-values.
	date_time_str=`${FILE_IO_UTILS_SCRIPT} -get_date_time_str $2`
	echo "${date_time_str} @ $1: Assigning final p-values."
	paste ${META_FULLDEC_POOLED_NOISY_S}.txt ${META_FULLDEC_POOLED_NOISY_T}.txt | awk 'BEGIN{FS="\t"}{print $2"\t"$4}' | ${COLLAGENE_SECURE_EXEC} -assign_chisqr_pvals_per_ST_stats stdin ${INTERM_DATA_DIR}/FINAL_META_P_VALUES_VAR_BLOCK_CLIENT_${site_i_per_cmd}.txt
	paste ${META_FULLDEC_POOLED_NOISY_S}.txt ${META_FULLDEC_POOLED_NOISY_T}.txt | awk 'BEGIN{FS="\t"}{print $2"\t"$4}' > ${INTERM_DATA_DIR}/META_ALL_SITE_NOISY_ST_STATISTICS.txt
	Rscript GET_CHISQ_PVAL.R ${INTERM_DATA_DIR}/META_ALL_SITE_NOISY_ST_STATISTICS.txt ${INTERM_DATA_DIR}/FINAL_META_P_VALUES_CLIENT_${site_i_per_cmd}_PER_R.txt

	# Paste the variant identifiers.
	cut -f1,1 ${LOCAL_DATA_DIR}/GMMAT_RESULTS.txt | awk {'if(NR>1){print $0}'} > ${INTERM_DATA_DIR}/meta_var_ids.txt
	paste ${INTERM_DATA_DIR}/meta_var_ids.txt ${INTERM_DATA_DIR}/FINAL_META_P_VALUES_CLIENT_${site_i_per_cmd}_PER_R.txt > FINAL_META_P_VALUES_CLIENT_${site_i_per_cmd}_PER_R_w_var_ids.txt

    exit 0
fi

if [[ ${cmd_option} == "-symmetric_encrypt_partdec_data" ]]
then
    if [[ $# != 5 ]]
    then
            echo "USAGE: $0 $1 -symmetric_encrypt_partdec_data [Data config file] [partdec file] [openssl hash file] [Output file]"
            exit
    fi

    data_config_file=$2
    partdec_file=$3
	openssl_hash_file=$4
    output_file=$5

	openssl enc -aes-256-cbc -in ${partdec_file} -out ${output_file} -pass file:${openssl_hash_file}
	exit $?
fi

if [[ ${cmd_option} == "-symmetric_decrypt_partdec_data" ]]
then
    if [[ $# != 5 ]]
    then
            echo "USAGE: $0 $1 -symmetric_decrypt_partdec_data [Data config file] [encrypted partdec file] [openssl hash file] [Output file]"
            exit
    fi

    data_config_file=$2
    enc_partdec_file=$3
	openssl_hash_file=$4
    output_file=$5

	openssl enc -d -aes-256-cbc -in ${enc_partdec_file} -out ${output_file} -pass file:${openssl_hash_file}
	exit $?
fi
