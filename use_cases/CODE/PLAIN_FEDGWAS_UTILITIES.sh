#!/bin/bash

if [[ $# -lt 2 ]]
then
	echo "USAGE: $0 [options] [arguments]
	Options:
		-split_GMMAT_dir [GMMAT directory] [N sites] [N vars 2 use] [Data input directory]
		-run_FedGWAS [Data input directory] [# epochs] [LL epsilon] [plink2 pvals file]"
	exit
fi

cmd_option=$1

COLLAGENE_SECURE_EXEC=COLLAGENE_Release

#SIGMOID_TYPE=KIM_ETAL
SIGMOID_TYPE=NATIVE
#SIGMOID_TYPE=TENSEAL

if [[ ${cmd_option} == "-split_GMMAT_dir" ]]
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 -split_GMMAT_dir [GMMAT directory] [N sites] [N vars 2 use] [Data input directory]"
		exit
	fi

	GMMAT_dir=$2
	N_SITES=$3
	N_VARS_2_USE=$4
	INPUT_DATA_DIR=$5

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
		multicolumn_processing_tools -extract_rows_per_query_column_preserve_query_order ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/phenotypes.txt 0 1 temp_pheno.txt
		echo -e "SAMPLE\tPHENOTYPE" > pheno_header.txt
		cat pheno_header.txt temp_pheno.txt > ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt

	        # Separate covariates.
		multicolumn_processing_tools -extract_rows_per_query_column_preserve_query_order ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/feat_matrix.txt 0 1 temp_feats.txt
		head -n 1 ${INPUT_DATA_DIR}/feat_matrix.txt > feat_header.txt
		cat feat_header.txt temp_feats.txt > ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt

		# Separate genotypes.
		echo "Extracting genotypes for site ${site_i}"
		PhenomeModel -extract_subsample_GMMAT_genotypes_from_GMMAT_genotype_matrix  ${INPUT_DATA_DIR}/genotypes.txt GMMAT_sample_ids.list ${INPUT_DATA_DIR}/SITE_${site_i}/subjects.list ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt

		exit
done
fi

if [[ ${cmd_option} == "-run_FedGWAS" ]]
then
	INPUT_DATA_DIR=$2
	N_EPOCHS=$3
	LL_EPSILON=$4
	plink2_pvals_fp=$5

	if [[ ! -d ${INPUT_DATA_DIR} ]]
	then
		echo "Could not find input data directory @ \"${INPUT_DATA_DIR}\""
		exit
	fi
	
	N_SITES=`ls ${INPUT_DATA_DIR} | grep SITE_ | wc -l | awk {'print $1'}`
	echo "Found ${N_SITES} sites under ${INPUT_DATA_DIR}"

	echo "Running FedGWAS with:
INPUT_DATA_DIR=$INPUT_DATA_DIR
N_SITES=$N_SITES
N_EPOCHS=$N_EPOCHS
LL_EPSILON=$LL_EPSILON
plink2_pvals_fp=$plink2_pvals_fp"

	site_i_list=`seq 1 $N_SITES`

	# Build output directory.
	SHARED_DIR=SHARED_DATA_MODEL_DIR
	rm -f -r ${SHARED_DIR}
	mkdir ${SHARED_DIR}

	########################################################
	# These are the pooled datasets.
	FEATS_MATRIX_FP=${INPUT_DATA_DIR}/feat_matrix.txt
	PHENO_FP=${INPUT_DATA_DIR}/phenotypes.txt
	MODEL_FP=${INPUT_DATA_DIR}/model.txt
	GENO_FP=${INPUT_DATA_DIR}/genotypes.txt
	
	if [[ -f "${FEATS_MATRIX_FP}" ]]
	then
		# Run R-glm 
		Rscript BUILD_LR.R ${INPUT_DATA_DIR} >& CENTRALIZED_R.OP

		# Compare parameters.
		#paste ${OP_DIR}/model.txt trained_model.txt
		awk {'if(NR>1){print $0}'} per_subject_phenotypes.txt > temp_pheno.txt
		paste R_preds.txt temp_pheno.txt | head -n 20
	fi
	########################################################

	epoch_i_list=`seq 0 ${N_EPOCHS}`

	echo "Fitting the null model."
	for cur_epoch in ${epoch_i_list[@]}
	do
		echo "Epoch ${cur_epoch}.."

		# Update the current parameters.
		for site_i in ${site_i_list[@]}
		do
			site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
			${COLLAGENE_SECURE_EXEC} -cryptable_client_calculate_save_XtWX_XtWz ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& XtWX_XtWz_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
		done

		for site_i in ${site_i_list[@]}
	        do
	                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
	                ${COLLAGENE_SECURE_EXEC} -cryptable_client_add_mult_noise_2_XtWX ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& ADD_ALL_SITE_MULT_NOISE_XtWX_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
	        done

	        for site_i in ${site_i_list[@]}
	        do
	                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
	                ${COLLAGENE_SECURE_EXEC} -cryptable_client_pool_site_specific_all_site_noise_XtWX ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& POOL_ALL_SITE_NOISE_XtWX_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
	        done

	        for site_i in ${site_i_list[@]}
	        do
	                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
	                ${COLLAGENE_SECURE_EXEC} -cryptable_collaborative_decrypt_pooled_noisy_XtWX ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& PARTDEC_XtWX_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
        	done

	        for site_i in ${site_i_list[@]}
	        do
	                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
	                ${COLLAGENE_SECURE_EXEC} -cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& POOL_PARTDEC_XtWX_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
	        done

		# Pool and save parameters.
		for site_i in ${site_i_list[@]}
		do
			site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
			${COLLAGENE_SECURE_EXEC} -cryptable_client_pool_XtWX_XtWz_update_beta ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& pool_XtWX_XtWz_CLIENT_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
		done

		for site_i in ${site_i_list[@]}
		do
			site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
			${COLLAGENE_SECURE_EXEC} -cryptable_client_collaborative_decrypt_beta ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& PARTDEC_BETA_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
		done

		for site_i in ${site_i_list[@]}
		do
			site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
			${COLLAGENE_SECURE_EXEC} -cryptable_client_pool_partially_decrypted_beta ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${N_EPOCHS} ${SIGMOID_TYPE} ${LL_EPSILON} ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& POOL_PARTDEC_BETA_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
		done

		for site_i in ${site_i_list[@]}
		do
			site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
			${COLLAGENE_SECURE_EXEC} -client_check_convergence_per_updated_beta ${cur_epoch} ${site_i_per_cmd} ${N_SITES} ${LL_EPSILON} ${SHARED_DIR} >& CHECK_CONVERGENCE_${site_i_per_cmd}_EPOCH_${cur_epoch}.op
		done	
	done # epoch loop.

	echo "CALCULATING P-VALUE STATISTICS FROM EACH SITE"
	VAR_BLOCK_SIZE=1000
	for site_i in ${site_i_list[@]}
	do
		echo "Site ${site_i}"

		site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`	
		${COLLAGENE_SECURE_EXEC} -cryptable_client_calculate_save_pvalue_stats ${site_i_per_cmd} ${N_SITES} ${N_EPOCHS} ${VAR_BLOCK_SIZE} ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& P_VAL_STAT_CLIENT_${site_i_per_cmd}.op
	done

	echo "POOLING P-VALUE STATS FROM ALL SITES"
	for site_i in ${site_i_list[@]}
	do
		echo "Site ${site_i}"

	        site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
	        ${COLLAGENE_SECURE_EXEC} -cryptable_client_pool_pvalue_stats ${site_i_per_cmd} ${N_SITES} ${N_EPOCHS} ${VAR_BLOCK_SIZE} ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& SCALE_UPDATE_CLIENT_${site_i_per_cmd}.op
	done

	echo "PARTIAL DECRYPTING P-VALUE STATS FROM ALL SITES"
        for site_i in ${site_i_list[@]}
        do
                echo "Site ${site_i}"

                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
                ${COLLAGENE_SECURE_EXEC} -cryptable_client_collaborative_decrypt_pval_stats ${site_i_per_cmd} ${N_SITES} ${N_EPOCHS} ${VAR_BLOCK_SIZE} ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& COLLABORATIVE_DECRYPT_P_VAL_STATS_${site_i_per_cmd}.op
        done

        echo "POOLING PARTIALLY DECRYPTED P-VALUE STATS FROM ALL SITES"
        for site_i in ${site_i_list[@]}
        do
                echo "Site ${site_i}"

                site_i_per_cmd=`echo ${site_i} | awk '{print $1-1}'`
                ${COLLAGENE_SECURE_EXEC} -cryptable_client_pool_partially_decrypted_pval_stats ${site_i_per_cmd} ${N_SITES} ${N_EPOCHS} ${VAR_BLOCK_SIZE} ${INPUT_DATA_DIR}/SITE_${site_i}/genotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i}/feat_matrix.txt ${INPUT_DATA_DIR}/SITE_${site_i}/phenotypes.txt ${INPUT_DATA_DIR}/SITE_${site_i} ${SHARED_DIR} >& POOL_PARTDEC_P_VAL_STATS_${site_i_per_cmd}.op

				cut -f4,5 ${SHARED_DIR}/full_dec_pval_matrix_iter_${N_EPOCHS}_client_${site_i_per_cmd}.bin > site_${site_i_per_cmd}_ST_stats.txt
                Rscript GET_CHISQ_PVAL.R site_${site_i_per_cmd}_ST_stats.txt site_${site_i_per_cmd}_p_vals_per_R.txt
                cut -f1,1 ${SHARED_DIR}/full_dec_pval_matrix_iter_${N_EPOCHS}_client_${site_i_per_cmd}.bin > temp_var_ids.txt
                paste temp_var_ids.txt site_${site_i_per_cmd}_p_vals_per_R.txt > site_${site_i_per_cmd}_p_vals_per_R_w_var_ids.txt
        done

		# Calculate p-values.
		chmod 755 GET_CHISQ_PVAL.R
		cut -f4,5 ${SHARED_DIR}/full_dec_pval_matrix_iter_${N_EPOCHS}_client_0.bin > plain_FedLR_ST.txt
		Rscript GET_CHISQ_PVAL.R plain_FedLR_ST.txt plain_FedLR_p_vals.txt
		cut -f1,1 ${SHARED_DIR}/full_dec_pval_matrix_iter_${N_EPOCHS}_client_0.bin > snp_ids.txt
		paste snp_ids.txt plain_FedLR_p_vals.txt > plain_FedLR_p_vals_w_var_ids.txt
fi
