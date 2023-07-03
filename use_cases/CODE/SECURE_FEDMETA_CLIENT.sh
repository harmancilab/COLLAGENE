#!/bin/bash

source data_config.params

if [[ $# -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
	-stop_FedGWAS_processes
	-clean_local_directory
	-validate_meta_run
	-run_META_analysis"

	exit
fi

running_site_i=${LOCAL_SITE_I}

cmd_option=$1

TIME_MEM_LOG_METAANALYSIS="META_ANALYSIS_TIME_MEM.log"
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
	rm -f ${TIME_MEM_LOG_METAANALYSIS}
	
	# Remove the file I/O temp directory, this is necessary to make sure that we can get a fresh copy of the results.
	rm -f -r ${FILE_IO_TEMP_DIR}

	find ${LOCAL_DATA_DIR} -name 'INTERMEDIATE' | xargs -Ifiles rm -f -r files

	exit
fi

site_iters=`seq 1 ${N_SITES}`

if [[ ${cmd_option} == "-stop_FedGWAS_processes" ]]
then
	echo "Stopping FedLR processes.."

	n_procs=2
	while [[ $n_procs -gt 1 ]]
	do
		ps aux | grep bash | grep SECURE_FEDGWAS_UTILITIES | awk {'print "kill -9 "$2'} > temp_kill.sh;
		ps aux | grep bash | grep SECURE_FEDMETA_CLIENT | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep bash | grep RUN_CLIENT_METAANALYSIS_COMMANDS | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep bash | grep FILE_IO_UTILITIES | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep scp | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		ps aux | grep COLLAGENE | awk {'print "kill -9 "$2'} >> temp_kill.sh;
		n_procs=`wc -l temp_kill.sh | awk {'print $1'}`
		chmod 755 temp_kill.sh;
		./temp_kill.sh
	done

	echo "${n_procs} left.."
fi

if [[ "${cmd_option}" == "-validate_meta_run" ]]
then
	echo "Checking input directories.."

	site_i_per_cmd=${running_site_i}

	echo "Checking ${LOCAL_DATA_DIR} for site-${site_i_per_cmd}"
	${SECURE_FEDGWAS_UTILITIES} -validate_meta_analysis_input_data data_config.params ${LOCAL_DATA_DIR}
	
	if [[ -f INPUT.ERROR ]]
	then
		echo "There is a problem with input data."
		exit 1
	fi

	exit
fi

if [[ "${cmd_option}" == "-run_META_analysis" ]]
then
	if [[ $# != 1 ]]
	then
		echo "$0 $1"
		exit
	fi

	running_site_i=${LOCAL_SITE_I}

	echo "Running meta-analysis of GMMAT results at site-${running_site_i} over ${N_SITES} sites.."

	# Run each epoch.
	echo "Encrypting meta-analysis statistics."
	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			# Check the existence of files for this site.

			echo "Site ${cur_site_i} (${site_i_per_cmd})"
			${time_exec} -o ${TIME_MEM_LOG_METAANALYSIS} --append -f ${cur_iter}"meta_encrypt_save_ST_stats\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -meta_encrypt_save_ST_stats data_config.params ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& META_STEP1_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Pooling S/T statistics and adding site-specific noise."
	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i} (${site_i_per_cmd})"
			${time_exec} -o ${TIME_MEM_LOG_METAANALYSIS} --append -f ${cur_iter}"meta_pool_per_site_ST_stats_add_noise\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -meta_pool_per_site_ST_stats_add_noise data_config.params ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& META_STEP2_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Pooling noisy S/T statistics from all sites."
	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i} (${site_i_per_cmd})"
			${time_exec} -o ${TIME_MEM_LOG_METAANALYSIS} --append -f ${cur_iter}"meta_pool_noisy_ST_stats\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -meta_pool_noisy_ST_stats data_config.params ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& META_STEP3_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Collaboratively decrypting pooled noisy S/T statistics."
	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i} (${site_i_per_cmd})"
			${time_exec} -o ${TIME_MEM_LOG_METAANALYSIS} --append -f ${cur_iter}"meta_collaborative_decrypt_pooled_noisy_ST_stats\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -meta_collaborative_decrypt_pooled_noisy_ST_stats data_config.params ${PRIVATE_KEY_FILE} ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& META_STEP4_${cur_site_i}.txt &
		fi
	done
	wait

	echo "Pooling partial decryptions of noisy S/T statistics."
	for cur_site_i in ${site_iters[@]}
	do
		site_i_per_cmd=`echo ${cur_site_i} | awk '{print $1-1}'`

		if [[ ${running_site_i} == ${site_i_per_cmd} ]]
		then
			echo "Site ${cur_site_i} (${site_i_per_cmd})"
			${time_exec} -o ${TIME_MEM_LOG_METAANALYSIS} --append -f ${cur_iter}"meta_pool_partially_decrypted_pooled_noisy_ST_stats\t"%e"\t"%M ${SECURE_FEDGWAS_UTILITIES} -meta_pool_partially_decrypted_pooled_noisy_ST_stats data_config.params ${site_i_per_cmd} ${LOCAL_DATA_DIR} >& META_STEP5_${cur_site_i}.txt &
		fi
	done
	wait
fi
	
