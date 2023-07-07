#!/bin/bash

if [[ ! -f "COLLAGENE.sh" ]]
then
	echo "Could not fine COLLAGENE.sh script, copy it from under scripts directory."
	exit 1
fi

# You can skip data generation to reuse the same data.
GENERATE_DATA=1
if [[ ${GENERATE_DATA} == 1 ]]
then
	n_K_clust=4
	n_dim=2
	./KMEANS_COMMANDS.sh -generate_data ${n_K_clust} ${n_dim}
fi

./KMEANS_COMMANDS.sh -setup_directories
./KMEANS_COMMANDS.sh -run_protocol
./KMEANS_COMMANDS.sh -plot_results

