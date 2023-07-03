#!/bin/bash

source data_config.params

if [[ $# -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
	-stop_FedGWAS_processes
	-clean_local_directory
	-validate_run
	-fit_null_model
	-assign_p_values"

	exit
fi

running_site_i=${LOCAL_SITE_I}

cmd_option=$1

TIME_MEM_LOG_FEDGWAS="SECURE_FEDGWAS_TIME_MEM.log"
time_exec=/usr/bin/time

if [[ "${cmd_option}" == "-clean_local_directory" ]]
then
	echo "Cleaning shared data and intermediate files."
	rm -f *.txt
	rm -f *.op
	rm -f *.bin
	rm -f *.enc
	rm -f *.partdec

	# Clean this here only.
	rm -f ${TIME_MEM_LOG_FEDGWAS}

	# Remove the file I/O temp directory, this is necessary to make sure that we can get a fresh copy of the results.
	rm -f -r ${FILE_IO_TEMP_DIR}

	find ${LOCAL_DATA_DIR} -name 'INTERMEDIATE' | xargs -Ifiles rm -f -r files

	exit
fi

epoch_iters=`seq 0 ${N_EPOCHS}`
site_iters=`seq 1 ${N_SITES}`

if [[ ${cmd_option} == "-stop_FedGWAS_processes" ]]
then
	echo "Stopping FedLR processes.."

	n_procs=2
	while [[ $n_procs -gt 1 ]]
	do
		ps aux | grep bash | grep SECURE_FEDGWAS_UTILITIES | awk {'print "kill -9 "$2'} > temp_kill.sh;
		ps aux | grep bash | grep SECURE_FEDGWAS_CLIENT | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep bash | grep RUN_CLIENT_FEDGWAS_COMMANDS | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep bash | grep FILE_IO_UTILITIES | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep scp | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep COLLAGENE | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		n_procs=`wc -l temp_kill.sh | awk {'print $1'}`
		chmod 755 temp_kill.sh;
		./temp_kill.sh
	done

	echo "${n_procs} left.."
fi

if [[ "${cmd_option}" == "-validate_run" ]]
then
	echo "Checking input directories.."

	site_i_per_cmd=${running_site_i}

	echo "Checking ${LOCAL_DATA_DIR} for site-${site_i_per_cmd}"
	${SECURE_FEDGWAS_UTILITIES} -validate_input_data data_config.params ${LOCAL_DATA_DIR}
	
	if [[ -f INPUT.ERROR ]]
	then
		echo "There is a problem with input data."
		exit 1
	fi

	exit
fi

if [[ "${cmd_option}" == "-fit_null_model" ]]
then
	if [[ $# != 1 ]]
	then
		echo "$0 $1"
		exit
	fi

	running_site_i=${LOCAL_SITE_I}

	echo "Running ${N_EPOCHS} iterations for null model fitting at site-${running_site_i}.."

	# Run each epoch.
	for cur_epoch in ${epoch_iters[@]}
	do
		echo "Processing epoch=${cur_epoch}"

		echo "Calculating XtWX and XtWz"
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				# Check the existence of files for this site.

				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_calculate_save_XtWX_XtWz\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -client_calculate_save_XtWX_XtWz data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP1_${cur_site_i}.txt 
			fi
		done
		wait

		echo "Pooling noise to local XtWX matrices."
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "add_noise_to_XtWX\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -add_noise_to_XtWX data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP2_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Pooling XtWX with pooled noise"
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "pool_site_specific_all_site_noise_XtWX\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -pool_site_specific_all_site_noise_XtWX data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP3_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Collaboratively decrypting pooled noisy XtWX"
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "collaborative_decrypt_pooled_noisy_XtWX\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -collaborative_decrypt_pooled_noisy_XtWX data_config.params ${PRIVATE_KEY_FILE} ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP4_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Pooling partial decryptions, inverting, removing noise from XtWX."
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "pool_partially_decrypted_pooled_noisy_XtWx_remove_noise\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -pool_partially_decrypted_pooled_noisy_XtWx_remove_noise data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP5_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Updating beta"
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "pool_XtWX_XtWz_update_beta\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -pool_XtWX_XtWz_update_beta data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP6_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Partially decrypting beta"
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_collaborative_decrypt_beta\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_collaborative_decrypt_beta data_config.params ${PRIVATE_KEY_FILE} ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP7_${cur_site_i}.txt &
			fi
		done
		wait

		echo "Pooling partial decryptions of beta" 
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_pool_partially_decrypted_beta\t${cur_epoch}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_pool_partially_decrypted_beta data_config.params ${cur_epoch} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& EPOCH_${cur_epoch}_STEP8_${cur_site_i}.txt &
			fi
		done
		wait
	done
fi

if [[ "${cmd_option}" == "-assign_p_values" ]]
then
	if [[ $# != 1 ]]
	then
		echo "$0 $1"
		exit
	fi

	running_site_i=${LOCAL_SITE_I}

	echo "Running ${N_EPOCHS} iterations for p-value assignment.."

	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i}"
			${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_calculate_save_pvalue_stats\t${N_EPOCHS}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_calculate_save_pvalue_stats data_config.params ${N_EPOCHS} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& STEP9_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Pooling p-value statistics."
    for cur_site_i in ${site_iters[@]}
    do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i}"
			${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_pool_pvalue_stats\t${N_EPOCHS}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_pool_pvalue_stats data_config.params ${N_EPOCHS} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& STEP10_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Pooling ST Statistics.."
        for cur_site_i in ${site_iters[@]}
        do            
            site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
				echo "Site ${cur_site_i}"
                ${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "pool_noisy_ST_stats\t${N_EPOCHS}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -pool_noisy_ST_stats data_config.params ${N_EPOCHS} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& STEP11_${cur_site_i}.txt &
			fi
        done
        wait

        echo "Partially decrypting noisy ST statistics."
        for cur_site_i in ${site_iters[@]}
        do
            site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
                echo "Site ${cur_site_i} (${site_i_per_cmd})"
                ${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_collaborative_decrypt_noisy_ST_stats\t${N_EPOCHS}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_collaborative_decrypt_noisy_ST_stats data_config.params ${PRIVATE_KEY_FILE} ${N_EPOCHS} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& STEP12_${cur_site_i}.txt &
			fi
        done
        wait

		echo "Pooling partially decrypted noisy ST statistics."
		for cur_site_i in ${site_iters[@]}
		do
			site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

			if [[ ${running_site_i} == ${site_i_per_cmd} ]]
			then
                echo "Site ${cur_site_i} (${site_i_per_cmd})"
				${time_exec} -o ${TIME_MEM_LOG_FEDGWAS} --append -f "client_pool_partially_decrypted_noisy_ST_stats\t${N_EPOCHS}\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES}  -client_pool_partially_decrypted_noisy_ST_stats data_config.params ${N_EPOCHS} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& STEP13_${cur_site_i}.txt &
				${SECURE_FEDGWAS_UTILITIES} -pool_final_plaintext_p_values data_config.params ${site_i_per_cmd} LOCAL_DATA_DIR/

			fi
	done
	wait
fi # -assign_p_values

	
