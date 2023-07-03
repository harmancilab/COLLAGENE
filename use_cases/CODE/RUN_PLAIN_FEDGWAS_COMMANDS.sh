if [[ ! -d "LOCAL_DATA_DIR" ]]
then
	echo "Could not find the local per-site data."
fi

./PLAIN_FEDGWAS_UTILITIES.sh -run_FedGWAS LOCAL_DATA_DIR 5 0.00001 none



