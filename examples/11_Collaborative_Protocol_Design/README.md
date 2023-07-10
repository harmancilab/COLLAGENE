## Structuring Scripts for Collaborative Analysis
One of the challanges in implementing collaborative analysis pipelines is synchronizing the sites and implementing the blocking calls. Producer-consumer design pattern can help with this challange.

The code snippet in network I/O can serve as a template for building protocols that use file-based communication. Generally, the protocols require all sites perform same operations on their local data and share the encrypted results with other sites.

Thus, from the perspective of any site, the script can be same without any hardcoded parameters except for the site index (0-based) and the number of sites.

This protocol script can be implemented on one of the sites and be reviewed and approved by other sites as the running protocol to be used by all sites.

Following is a pseudo-bash-script that demonstrates this idea:
```
i_site=[Input site index]
N_sites=[# sites]

# Generate the intermediate matrix using site-specific local data.
GENERATE_INTERMEDIATE_MATRIX.sh ${i_site} matrix_from_${i_site}.enc

# Upload the local data so that it can be processed by other sites.
echo "matrix_from_${i_site}.enc" > files.list
./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params files.list

# All sites should have uploaded the intermediate matrix to be processed by other sites.
# Processe all the intermediate files from other sites.
for i_producer_site in {1..$N_sites}
do
	echo "matrix_from_${i_site}.enc_processed_by_${i_consumer_site}.enc" > files.list
	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list

	./FILE_IO_UTILITIES.sh -download_files_to_shared data_config.params files.list

	PROCESS_INTERMEDIATE_MATRIX.sh matrix_from_${i_producer_site}.enc matrix_from_${i_producer_site}.enc_processed_by_${i_site}.enc

	echo "matrix_from_${i_producer_site}.enc_processed_by_${i_site}.enc" > files.list
	./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params files.list
done

# After files are processed, wait for other sites to process this site's data, download, and use for further processing. 
# Note below how the site indexing is changed with respect to previous loop.
for i_consumer_site in {1..$N_sites}
do
	# Wait for the files to be processed by other sites.
	echo "matrix_from_${i_site}.enc_processed_by_${i_consumer_site}.enc" > files.list
	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list

	# Download the files.
	./FILE_IO_UTILITIES.sh -download_files_to_shared data_config.params files.list
done

# Further processing on the results using the files that are downloaded.
```

Note that it is generally a good idea to add error checks after all calls. All *FILE_IO_UTILS.sh* commands should return 0 on success and 1 on an error. After running any of the options, it is generally a good practice to check for errors. For example:
```
...
	echo "matrix_from_${i_producer_site}.enc_processed_by_${i_site}.enc" > files.list
	./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params files.list
	
	if [[ $? == 1 ]]
	then
		echo "Failed upload of the processed file."

		... Handle error condition; the upload can be retried or an fatal error message can be written ...
	fi

...
```

This is implemented in network I/O example script under previous example folder.

__Redundancy of Network I/O__: Above design idea generates redundant network I/O since each producer site downloads its produced files. These cases can be avoided by adding a conditional to make sure the sites do not download/upload their own data. Note that this does not create security concerns; it is only necessary to decrease the network I/O.

