#!/bin/bash

#./SECURE_FEDGWAS_CLIENT.sh -stop_FedGWAS_processes

./SECURE_FEDMETA_CLIENT.sh -clean_local_directory
./SECURE_FEDMETA_CLIENT.sh -validate_meta_run

if [[ $? != 0 ]]
then
	echo "Input validity check failed.."
	exit 1
fi

./SECURE_FEDMETA_CLIENT.sh -run_META_analysis >& ALL_META.OP

