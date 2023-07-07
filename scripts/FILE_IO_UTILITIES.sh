#!/bin/bash 

# This script implements the secure COLLAGENE pipeline; Runs all sites in sequential order in loops.
# TODO::Use mktemp -d for downloading the metadata::mktemp -d --tmpdir=$PWD
#

COLLAGENE_SECURE_EXEC=COLLAGENE_Release

if [[ $# -lt 2 ]]
then
	echo "USAGE: $0 [options] [arguments]
		-get_resource_size [Data config file] [Resource path]
		-get_date_time_str [Data config file]
		-clean_shared_directory [Data config file]
		-test_file_IO [Data config file] 
		-probe_files_in_shared [Data config file] [File/Directory list]
		-wait_for_files_in_shared [Data config file] [File/Directory list]
		-download_files_from_shared [Data config file] [File/Directory list]
		-upload_files_to_shared [Data config file] [File/Directory list]"

	exit
fi

cmd_option=$1
data_config_file=$2

if [[ ! -f ${data_config_file} ]]
then
	echo "Could not find data config file @ \"${data_config_file}\""
	exit
fi

source ${data_config_file}

# Check to make sure the temp IO directory exists.
if [[ ! -d ${FILE_IO_TEMP_DIR} ]]
then
	echo "Could not find ${FILE_IO_TEMP_DIR}, creating it.."
	mkdir ${FILE_IO_TEMP_DIR}

	if [[ ! -d ${FILE_IO_TEMP_DIR} ]]
	then
		echo "Could not create ${FILE_IO_TEMP_DIR}, exiting.."
		exit
	fi
fi

# Verify COLLAGENE executable.
collagene_check=`type -P ${COLLAGENE_SECURE_EXEC}`
if [[ "${collagene_check}" == "" ]]
then
	echo "Could not find COLLAGENE executable."
	exit 1
fi

if [[ ${cmd_option} == "-get_date_time_str" ]]
then
	date_time_str=`date "+[%m/%d/%y - %H:%M:%S]"`
	echo $date_time_str

	exit 0
fi

if [[ ${cmd_option} == "-get_resource_size" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [data config path] [Resource path (file/directory)]"
		exit 1
	fi

	resource_path=$3

	if [[ -f ${resource_path} ]]
	then
		file_size=`stat --printf=%s ${resource_path}`
		echo ${file_size}
		exit 0
	fi

	if [[ -d ${resource_path} ]]
	then
		file_size=`find ${resource_path} |  xargs -Ifiles stat --printf='%s\n' files | awk 'BEGIN{tot_size=0}{tot_size+=$1}END{print tot_size}'`
		echo ${file_size}
		exit 0
	fi

	# Could not figure out the type of resource.
	exho "-1"
	exit 1
fi

# -test_file_IO with the selected IO option.
if [[ ${cmd_option} == "-test_file_IO" ]]
then
	if [[ $# != 2 ]]
	then
		echo "USAGE: $0 $1 [Data config file]"
		exit 1
	fi

	data_config_file=$2

	# Verify COLLAGENE executable.
	collagene_check=`type -P ${COLLAGENE_SECURE_EXEC}`
	if [[ "${collagene_check}" == "" ]]
	then
		echo "Could not find COLLAGENE executable."
		exit 1
	fi

	# Upload/download a test file.
	test_file=random_matrix_${RANDOM}.bin
	${COLLAGENE_SECURE_EXEC} -generate_mult_full_noise_matrix 10 10 ${test_file}

	# Check key permissions if SCP is used.
	if [[ ${IO_TYPE} == "SCP" ]]
	then
		if [[ ! -f ${SCP_FEDERATION_SECRET_KEY} ]]
		then
			echo "Could not find federation IO key @ ${SCP_FEDERATION_SECRET_KEY}"
			exit 1
		fi

		fed_key_permission=`stat -c "%a" ${SCP_FEDERATION_SECRET_KEY}`
		if [[ ${fed_key_permission} != "600" ]]
		then
			echo "Key permission is not set correctly: ${fed_key_permission}, should be set to 600: \"chmod 600 ${SCP_FEDERATION_SECRET_KEY}\""
			exit 1
		fi

		# Do not use the root directories.
		if [[ "${SCP_REMOTE_SHARED_DIR}" == "/" ]]
		then
			echo "The remote shared directory seems to be the root directory, cannot use the root directory as the shared directory."
			exit 1
		fi

	elif [[ ${IO_TYPE} == "LOCAL" ]]
	then
		if [[ ! -d ${LOCAL_REMOTE_SHARED_DIR} ]]
		then
			echo "Could not find local ${LOCAL_REMOTE_SHARED_DIR}"
			exit 1
		fi

		# Do not use the root directories.
		if [[ "${LOCAL_REMOTE_SHARED_DIR}" == "/" ]]
		then
			echo "The local shared directory seems to be the root directory, cannot use the root directory as the shared directory."
			exit 1
		fi
	elif [[ ${IO_TYPE} == "S3" ]]
	then
		aws_exec_check=`type -P ${AWS_CLI_PATH}`
		if [[ "${aws_exec_check}" == "" ]]
		then
			echo "Could not find AWS executable @ ${AWS_CLI_PATH}"
			exit 1
		fi
	else
		echo "Could not figure out the IO type: ${IO_TYPE}"
		exit 1
	fi


	echo ${test_file} > temp.list

	$0 -upload_files_to_shared ${data_config_file} temp.list
	last_ret=$?
	echo "Result of upload is ${last_ret}"

	if [[ $last_ret != 0 ]]
	then
		exit $last_ret
	fi

	$0 -probe_files_in_shared ${data_config_file} temp.list
	last_ret=$?
	echo "Result of probing is ${last_ret}"

	if [[ $last_ret != 0 ]]
	then
		exit $last_ret
	fi

	mkdir temp_test_dir
	echo ${test_file} 
	$0 -download_files_from_shared ${data_config_file} temp.list temp_test_dir
	last_ret=$?
	echo "Result of download is ${last_ret}"

	if [[ $last_ret != 0 ]]
	then
		exit $last_ret
	fi

	echo "Difference between original and downloaded files: (Should be empty)"
	diff temp_test_dir/${test_file} ${test_file}

	echo "Everything looks good.."

	exit $last_ret
fi

if [[ ${cmd_option} == "-clean_shared_directory" ]]
then
    if [[ $# != 2 ]]
    then
        echo "USAGE: $0 $1 [Data config file]"
        exit
    fi

	source ${data_config_file}

	echo "Cleaning all files on the shared directory -- ${SCP_REMOTE_SHARED_DIR}. Run the following command:"

	if [[ ${IO_TYPE} == "SCP" ]]
	then
		if [[ "${LOCAL_REMOTE_SHARED_DIR}" == "/" ]]
		then
			echo "The local shared directory seems to be the root directory, cannot clean the root directory."
			exit 1
		fi

		echo "BE VERY CAREFUL WHILE RUNNING THE FOLLOWING COMMAND, MAKE SURE IT DOES NOT DELETE ANY IMPORTANT DIRECTORIES:"

		echo ${SSH_PREAMBLE} ${SCP_HOST} "rm -f -r ${SCP_REMOTE_SHARED_DIR}/*"
		last_ret=$?
		exit ${last_ret}
	elif [[ ${IO_TYPE} == "LOCAL" ]]
	then
		if [[ "${LOCAL_REMOTE_SHARED_DIR}" == "/" ]]
		then
			echo "The local shared directory seems to be the root directory, cannot clean the root directory."
			exit 1
		fi

		echo "BE VERY CAREFUL WHILE RUNNING THE FOLLOWING COMMAND, MAKE SURE IT DOES NOT DELETE ANY IMPORTANT DIRECTORIES:"

		echo "rm -f -r ${LOCAL_REMOTE_SHARED_DIR}/*"
		last_ret=$?
		exit ${last_ret}
	elif [[ ${IO_TYPE} == "S3" ]]
	then
		echo "${AWS_CLI_PATH} s3 rm --recursive ${S3_REMOTE_SHARED_DIR_KEY}" 
		last_ret=$?
		exit ${last_ret}
	fi
fi

if [[ ${cmd_option} == "-probe_files_in_shared" ]]
then
    if [[ $# != 3 ]]
    then
        echo "USAGE: $0 $1 [Data config file] [file/directory list]"
        exit 1
    fi

    source ${data_config_file}

    data_config_file=$2
    FILE_DIR_LIST_FILE=$3

	if [[ ${IO_TYPE} == "SCP" ]]
	then
		echo "Using SCP-based file I/O" > /dev/null
	elif [[ ${IO_TYPE} == "LOCAL" ]]
	then
		if [[ ! -d ${LOCAL_REMOTE_SHARED_DIR} ]]
		then
			echo "Could not find local ${LOCAL_REMOTE_SHARED_DIR}"
			exit 1
		fi
	elif [[ ${IO_TYPE} == "S3" ]]
	then
		echo "Using S3 file I/O" >> /dev/null
	else 
		echo "Could not figure out the IO type: ${IO_TYPE}"
		exit
	fi

	# Get a new temp directory for this request.
	CUR_FILE_IO_TEMP_DIR=`mktemp -d --tmpdir=${FILE_IO_TEMP_DIR} --suffix=${cmd_option}`
	if [[ ! -d ${CUR_FILE_IO_TEMP_DIR} ]]
	then
		echo "Could not create the temp directory: ${CUR_FILE_IO_TEMP_DIR}"
		exit 1
	fi

	files=`cat ${FILE_DIR_LIST_FILE}`

	# Check to make sure we have the list file.
	if [[ ! -f ${FILE_DIR_LIST_FILE} ]]
	then
		echo "Could not find the list of files to probe @ \"${FILE_DIR_LIST_FILE}\""
		exit 1
	fi # file list check.


	for cur_file in ${files[@]}
	do
		date_time_str=`$0 -get_date_time_str $2`
		echo "${date_time_str} Probing ${cur_file}"

		base_filename=`echo ${cur_file} | xargs -Ifile basename file`
		if [[ "${base_filename}" != "${cur_file}" ]]
		then
			echo "The file name ${cur_file} are not stripped, they must be stripped."
			exit 1
		fi

		echo "Downloading manifest for ${cur_file}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		manifest_file=${base_filename}.upload_manifest			
		# Get download command.
		download_cmd=""
		if [[ ${IO_TYPE} == "SCP" ]]
		then
			download_cmd="${SCP_PREAMBLE} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "LOCAL" ]]
		then
			download_cmd="cp -r ${LOCAL_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "S3" ]]
		then
			download_cmd="${AWS_CLI_PATH} s3 cp ${S3_REMOTE_SHARED_DIR_KEY}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		fi
	
		echo "Trying: ${download_cmd}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		${download_cmd} >& /dev/null

		echo "Success of manifest check for file ${cur_file}: $?" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op

		# Check the existence of the manifest file.
		if [[ ! -f ${CUR_FILE_IO_TEMP_DIR}/${manifest_file} ]]
		then
			exit 1
		fi
	done # file loop.
     
	# All files are found!
	exit 0
fi

# This option blocks running.
if [[ ${cmd_option} == "-wait_for_files_in_shared" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [File list to check]"
		exit 1
	fi

	data_config_file=$2
	FILE_DIR_LIST_FILE=$3

	# Check to make sure we have the list file.
	if [[ ! -f ${FILE_DIR_LIST_FILE} ]]
	then
		echo "Could not find the list of files to probe @ \"${FILE_DIR_LIST_FILE}\""
		exit
	fi # file list check.

	cur_status=1
	while [[ $cur_status == 1 ]]
	do
		# Update the current status.
		$0 -probe_files_in_shared ${data_config_file} ${FILE_DIR_LIST_FILE}
		cur_status=$?

		# Wait for 2 second.
		sleep 2
	done

	# return 0 error code.
	exit 0
fi

if [[ ${cmd_option} == "-upload_files_to_shared" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [File list to upload]"
		exit 1
	fi

	data_config_file=$2
	file_list_2_upload=$3

	if [[ ${IO_TYPE} == "SCP" ]]
	then
		echo "Using SCP-based file I/O" > /dev/null
	elif [[ ${IO_TYPE} == "LOCAL" ]]
	then
		if [[ ! -d ${LOCAL_REMOTE_SHARED_DIR} ]]
		then
			echo "Could not find local ${LOCAL_REMOTE_SHARED_DIR}"
			exit 1
		fi
	elif [[ ${IO_TYPE} == "S3" ]]
	then
		echo "Using S3 file I/O" >> /dev/null
	else
		echo "Could not figure out the IO type: ${IO_TYPE}"
		exit
	fi

	# Get a new temp directory for this request.
	CUR_FILE_IO_TEMP_DIR=`mktemp -d --tmpdir=${FILE_IO_TEMP_DIR} --suffix=${cmd_option}`
	if [[ ! -d ${CUR_FILE_IO_TEMP_DIR} ]]
	then
		echo "Could not create the temp directory: ${CUR_FILE_IO_TEMP_DIR}"
		exit 1
	fi

	files=`cat ${file_list_2_upload}`

	# The idea is to list all files and upload them one-by-one.
	for cur_file in ${files[@]}
	do
		date_time_str=`$0 -get_date_time_str $2`
		echo "${date_time_str} Uploading ${cur_file} -- ${CUR_FILE_IO_TEMP_DIR}"

		# Get the base name of the file to upload to the base. The local directory structure is not restored in the remote host.
		base_filename=`echo ${cur_file} | xargs -Ifile basename file`

		echo "Downloading manifest for ${cur_file}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		manifest_file=${base_filename}.upload_manifest

		# Get download command.
		download_cmd=""
		if [[ ${IO_TYPE} == "SCP" ]]
		then
			download_cmd="${SCP_PREAMBLE} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "LOCAL" ]]
		then
			download_cmd="cp ${LOCAL_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "S3" ]]
		then
			download_cmd="${AWS_CLI_PATH} s3 cp ${S3_REMOTE_SHARED_DIR_KEY}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		fi
			
		echo "Trying: ${download_cmd}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		${download_cmd} >& /dev/null
		manifest_download_res=$?

		# We don't check the status here since the file may just not be there. The status is checked later on.
		#if [[ $manifest_download_res != 0 ]]
		#then
		#	exit 1
		#fi

		echo "Success of manifest check for uploading file ${cur_file}: $?" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op

		# Check the existence of the downloaded manifest file, which may not be there, this means file can be uploaded.
		if [[ ! -f ${CUR_FILE_IO_TEMP_DIR}/${manifest_file} ]]
		then
			echo "Processing ${cur_file}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
			upload_cmd=""
			if [[ ${IO_TYPE} == "SCP" ]]
			then
				upload_cmd="${SCP_PREAMBLE} -r ${cur_file} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}/${base_filename}"
			elif [[ ${IO_TYPE} == "LOCAL" ]]
			then
				upload_cmd="cp -r ${cur_file} ${LOCAL_REMOTE_SHARED_DIR}/${base_filename}"
			elif [[ ${IO_TYPE} == "S3" ]]
			then
				rec_flag=""
				if [[ -d ${cur_file} ]]
				then
					rec_flag="--recursive"
				fi

				upload_cmd="${AWS_CLI_PATH} s3 cp ${rec_flag} ${cur_file} ${S3_REMOTE_SHARED_DIR_KEY}/${base_filename}"
			fi
	
			echo "Uploading file: " ${upload_cmd} >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
			${upload_cmd} >& /dev/null
			file_upload_res=$?
			echo "cmd result: ${file_upload_res}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op

			# If the file is uploaded, upload the manifest. Note that this depends on the success of the upload.
			if [[ ${file_upload_res} == 0 ]]
			then
				# Upload the manifest.
				manifest_file=${CUR_FILE_IO_TEMP_DIR}/${base_filename}.upload_manifest

				# Write file information to the manifest. This is used to evaluate the amount of network traffice we generate.
				#file_size=`stat --printf="%s" ${cur_file}`
				file_size=`$0 -get_resource_size $2 ${cur_file}`
				echo -e "${base_filename}\t${file_size}" > ${manifest_file}

				manifest_upload_cmd=""
				if [[ ${IO_TYPE} == "SCP" ]]
				then
					manifest_upload_cmd="${SCP_PREAMBLE} -r ${manifest_file} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}"
				elif [[ ${IO_TYPE} == "LOCAL" ]]
				then
					manifest_upload_cmd="cp -r ${manifest_file} ${LOCAL_REMOTE_SHARED_DIR}"
				elif [[ ${IO_TYPE} == "S3" ]]
				then
					manifest_upload_cmd="${AWS_CLI_PATH} s3 cp ${manifest_file} ${S3_REMOTE_SHARED_DIR_KEY}/${base_filename}.upload_manifest"
				fi

				echo "Uploading manifest: " ${upload_cmd} >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
				${manifest_upload_cmd} >& /dev/null
				manifest_upload_res=$?

				# This is the final manifest upload check that we skipped above to make sure manifest is uploaded correctly.
				if [[ ${manifest_upload_res} != 0 ]]
				then
					exit 1
				fi
			else
				exit 1
			fi
		else
			echo "Manifest already exists for ${cur_file}.."
			exit 1
		fi
	done # file list loop.


	# No errors, everything is upploaded.
	exit 0
fi


if [[ ${cmd_option} == "-download_files_from_shared" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Data config file] [File list to download] [Local directory to save the files]"
		exit 1
	fi

	data_config_file=$2
	file_list_2_download=$3
	local_directory_2_download=$4

	if [[ ${IO_TYPE} == "SCP" ]]
	then
		echo "Using SCP-based file I/O" > /dev/null
	elif [[ ${IO_TYPE} == "LOCAL" ]]
	then
		if [[ ! -d ${LOCAL_REMOTE_SHARED_DIR} ]]
		then
			echo "Could not find local ${LOCAL_REMOTE_SHARED_DIR}"
			exit 1
		fi
	elif [[ ${IO_TYPE} == "S3" ]]
	then
		echo "Using S3 file I/O" >> /dev/null
	else
		echo "Could not figure out the IO type: ${IO_TYPE}"
		exit
	fi

	# Get a new temp directory for this request.
	CUR_FILE_IO_TEMP_DIR=`mktemp -d --tmpdir=${FILE_IO_TEMP_DIR} --suffix=${cmd_option}`
	if [[ ! -d ${CUR_FILE_IO_TEMP_DIR} ]]
	then
		echo "Could not create the temp directory: ${CUR_FILE_IO_TEMP_DIR}"
		exit 1
	fi

	files=`cat ${file_list_2_download}`

	# Go over all the files and download them.
	for cur_file in ${files[@]}
	do
		date_time_str=`$0 -get_date_time_str $2`
		echo "${date_time_str} Downloading ${cur_file} -- ${CUR_FILE_IO_TEMP_DIR}"

		base_filename=`echo ${cur_file} | xargs -Ifile basename file`

		echo "Downloading manifest for ${cur_file}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		manifest_file=${base_filename}.upload_manifest

		# Get the download command.
		download_cmd=""
		if [[ ${IO_TYPE} == "SCP" ]]
		then
			download_cmd="${SCP_PREAMBLE} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "LOCAL" ]]
		then
			download_cmd="cp -r ${LOCAL_REMOTE_SHARED_DIR}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		elif [[ ${IO_TYPE} == "S3" ]]
		then
			download_cmd="${AWS_CLI_PATH} s3 cp ${S3_REMOTE_SHARED_DIR_KEY}/${manifest_file} ${CUR_FILE_IO_TEMP_DIR}"
		fi

		echo "Trying manifest download: ${download_cmd}" >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op
		${download_cmd} >& /dev/null
		manifest_download_res=$?

		if [[ ${manifest_download_res} != 0 ]]
		then
			echo "Manifest download failed for ${cur_file}"
			exit 1
		fi

		# If the manifest file is downloaded, download the actual file.
		if [[ -f ${CUR_FILE_IO_TEMP_DIR}/${manifest_file} ]]
		then
			# Get the download command.
			if [[ ${IO_TYPE} == "SCP" ]]
			then
				download_cmd="${SCP_PREAMBLE} ${SCP_HOST}:${SCP_REMOTE_SHARED_DIR}/${base_filename} ${local_directory_2_download}"
			elif [[ ${IO_TYPE} == "LOCAL" ]]
			then
				download_cmd="cp -r ${LOCAL_REMOTE_SHARED_DIR}/${base_filename} ${local_directory_2_download}"
			elif [[ ${IO_TYPE} == "S3" ]]
			then
				rec_flag=""
				#if [[ -d ${cur_file} ]]
				#then
				#	rec_flag="--recursive"
				#fi

				download_cmd="${AWS_CLI_PATH} s3 cp ${rec_flag} ${S3_REMOTE_SHARED_DIR_KEY}/${base_filename} ${local_directory_2_download}"
			fi

			echo "FOUND MANIFEST, DOWNLOADING: ${download_cmd}.." >> ${CUR_FILE_IO_TEMP_DIR}/STATS.op

			${download_cmd} >& /dev/null
			file_download_res=$?

			# This is a special case for S3 directories, it tries to download now as a directory.
			if [[ ${file_download_res} != 0 ]]
			then
				if [[ ${IO_TYPE} == "S3" ]]
				then
					date_time_str=`$0 -get_date_time_str $2`
					echo "${date_time_str} Note: Trying S3 download as a directory."
					rec_flag="--recursive"

					# Create the directory from the file name, we do not need to really create the directory.
					mkdir ${local_directory_2_download}/${base_filename}

					# Download.
					download_cmd="${AWS_CLI_PATH} s3 cp ${rec_flag} ${S3_REMOTE_SHARED_DIR_KEY}/${base_filename} ${local_directory_2_download}/${base_filename}"

					${download_cmd} >& /dev/null
					file_download_res=$?
				fi
			fi

			# Check for errors.
			if [[ ${file_download_res} != 0 ]]
			then
				echo "file download failed for ${cur_file}"
				exit 1
			fi
		else
			echo "Could not find the manifest for ${cur_file}"
			exit 1
		fi
	done # file list loop.

	exit 0
fi

echo "Unknown option "
exit 1