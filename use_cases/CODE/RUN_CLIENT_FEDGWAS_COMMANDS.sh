#!/bin/bash

#./SECURE_FEDGWAS_CLIENT.sh -stop_FedGWAS_processes

./SECURE_FEDGWAS_CLIENT.sh -clean_local_directory
./SECURE_FEDGWAS_CLIENT.sh -validate_run

if [[ $? != 0 ]]
then
	echo "Input validity check failed.."
	exit 1
fi

./SECURE_FEDGWAS_CLIENT.sh -fit_null_model >& ALL.OP
./SECURE_FEDGWAS_CLIENT.sh -assign_p_values >& ALL_p_vals.OP

