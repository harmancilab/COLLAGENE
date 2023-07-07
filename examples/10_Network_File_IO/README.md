# Network I/O (network_file_io_example.sh, FILE_IO_UTILITIES.sh)

Network I/O interface aims to simplify the network communication between parties when they are exchanging data sets using file-based communications. 

While these communications can be established by other approaches such as RPC, IPC, sockets etc. over network, these require real-time online presence of all processes and complicates the software design, application development, and software maintanence.

COLLAGENE currently supports sharing of encrypted data files to share intermediate data among collaborating sites.

The encrypted intermediate data is stored in a remote server (AWS, SSH) on a shared directory. COLLAGENE's *FILE_IO_UTILITIES.sh* script implements the options to download/upload/probe files in the remote server.

COLLAGENE supports 3 options for file I/O:
<ol>
<li> Local: Data is stored at a local folder. This option can be used for simulating and testing collaborative protocols by running all sites on the same client (Clients run in their folders and the shared directory is on the same computer.)</li> 
<li> SCP Server: Data is stored at an SFTP server that is accessed via SCP command with accession keys. COLLAGENE requires a key-based authentication to automate file transfers. This requires the sites to have authorized access to a shared ssh server. These keys can be setup using *ssh-keygen* and *ssh-copy-id* commands:</li> 

	# Run ssh-keygen to generate a public/private key pair for authentication, note that the passphrase should be skipped at this step.
	ssh-keygen -f test_key
	
	# Now, setup the key:
	ssh-copy-id -i test_key user_id@127.0.0.1
	
	# Replace "user_id" and "127.0.0.1" with username and server address.

<li> AWS S3 storage: Data is stored at an S3 bucket. It is necessary to setup the S3 bucket to transfer data to the bucket. The sites must have access to the bucket with write/read permissions. AWS CLI must be configured at each site (https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html).</li> 
</ol>

There are examples of using these options in *network_file_io_example.sh* script. This script help setup the basic functions. Also, *API/* directory contains the documentation for the options that *FILE_IO_UTILS.sh* implements.

## Design Pattern for Network I/O in Collaborative Data Analysis Pipelines
Network I/O is established using a producer-consumer type of model where files are produced by different sites and "consumed" by other sites without being deleted.

The interface implements the necessary upload/download, probing and waiting options. Some of these options has to block to make sure the producer generated the necessary file.

## Basic Working Idea of Network I/O
FILE_IO_UTILS.sh implements the basic functionalities of file-based communications. 

Each file is uploaded with a corresponding manifest file, which indicates existence of the file in the shared storage.

The communications must block for certain functions, e.g., when the sites need a file from a site, they wait for a file to become available. 

## Example Protocol ("networked_partial_decryption_example.sh"):
When sites are doing collaborative decryption, sites must exchange partial decryptions using the network I/O functions.

We include below one pseudo-script with explanations to demonstrate how one common script can be used to perform networked partial decryption among 3 sites.
```
# All sites perform the same set of steps in this protocol except for the different site indices below.
i_site=$1
N_sites=3

# The random matrix can represent an intermediate result (e.g., masked feature correlation matrix) that will be used by sites to perform an 
./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 matrix_from_site${i_site}.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} matrix_from_site0.bin matrix_from_site${i_site}.bin.enc

# Site uploads the matrix to be decrypted so other sites can download and partially decrypt it.
echo "matrix_from_site0.bin.enc" > files.list 
./FILE_IO_UTILITIES.sh -upload_files_to_shared files.list 

# Partially decrypt all sites: Loop over 0-{N_sites-1} (This )
for i_source_site in {0..${N_sites_min_one}}
do
	# Site partially decrypts and uploads the data to shared directory.
	echo "matrix_from_site${i_site}.bin.enc" > files.list 
	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list 
	./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params files.list 

	./COLLAGENE.sh -partial_decrypt_matrix matrix_from_site0.enc SITE_${i_site} matrix_from_site${i_source_site}.bin.enc 0 40 matrix_from_site${i_source_site}.bin.enc_partdec_by_${i_site}.partdec
	./COLLAGENE.sh -symmetric_encrypt_partdec_data SITE_${i_site} matrix_from_site${i_source_site}.bin.enc_partdec_by_${i_site}.partdec SITE_${i_site}/partdec_data_enc_hash.symmetric_key matrix_from_site${i_source_site}.bin.enc_partdec_by_${i_site}.partdec.enc
done

# At this point, site download the partial decryptions for its matrix:
for i_decrypting_site in {0..${N_sites_min_one}}
do
	echo "matrix_from_site${i_site}.bin.enc_partdec_by_${i_site}.partdec.enc" > files.list
	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list

	# Download the file.
	./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params files.list

	# Remove symmetric encryption from the sites using partdec encryption key.
	./COLLAGENE.sh -symmetric_decrypt_partdec_data SITE_${i_site} matrix_from_site${i_site}.bin.enc_partdec_by_${i_decrypting_site}.partdec.enc SITE_${i_site}/partdec_data_enc_hash.symmetric_key matrix_from_site${i_site}.bin.enc_partdec_by_${i_decrypting_site}.partdec
done

# Pool partial decryptions:
ls matrix_from_site${i_site}.enc_partdec_by_*.partdec > site_${i_site}.partdecs
./COLLAGENE.sh -pool_partially_decrypted_matrices SITE_${i_site} site_${i_site}.partdecs full_decrypted_matrix_from_site${i_site}.bin

# Write the matrix to text.
./COLLAGENE.sh -save_matrix_text full_decrypted_matrix_from_site${i_site}.bin full_decrypted_matrix_from_site${i_site}.bin.txt
./COLLAGENE.sh -save_matrix_text matrix_from_site${i_site}.bin matrix_from_site${i_site}.bin.txt

echo "Original matrix:"
cat matrix_from_site${i_site}.bin.txt

echo "Decrypted matrix:"
cat full_decrypted_matrix_from_site${i_site}.bin.txt
```

The full implementation of this script can be found locally with name "networked_partial_decryption_example.sh" (including error checks). 

To run this script, it is necessary to put it in different directories (including the key directories) for each site.

```
# Copy the scripts and keys.
rm -f -r SITE_0_DIR SITE_1_DIR SITE_2_DIR
# Setup the working directories for the sites:
mkdir SITE_0_DIR SITE_1_DIR SITE_2_DIR
cp -r networked_partial_decryption_example.sh COLLAGENE.sh FILE_IO_UTILITIES.sh data_config.params SITE_0 SITE_0_DIR
cp -r networked_partial_decryption_example.sh COLLAGENE.sh FILE_IO_UTILITIES.sh data_config.params SITE_1 SITE_1_DIR
cp -r networked_partial_decryption_example.sh COLLAGENE.sh FILE_IO_UTILITIES.sh data_config.params SITE_2 SITE_2_DIR
```

Before starting the scripts, make sure the shared directory is clean:
```
./FILE_IO_UTILITIES.sh -clean_shared_directory data_config.params
```
This option does not immediately delete the files in shared directory. Instead, it write the command to run for cleaning the shared directory. This is done on purpose to make sure users do not accidentally delete important files.

After cleaning the shared directory, we can start all site's scripts:
```
# Start all scripts.
cd SITE_0_DIR
nohup ./networked_partial_decryption_example.sh 0 >& SITE0.op &
cd ../SITE_1_DIR
nohup ./networked_partial_decryption_example.sh 1 >& SITE1.op &
cd ../SITE_2_DIR
nohup ./networked_partial_decryption_example.sh 2 >& SITE2.op &
cd ..
```

You can try commenting one of the site's to see its effect on the protocol; the protocol blocks at the partial decryption stage since the sites will wait for the encrypted matrices from all sites to be available.

This can be simply resolved by starting the site.

---

*network_file_io_example.sh*, *FILE_IO_UTILITIES.sh* contain more extensive examples and documentation of network file transfers that can be used for sharing data among collaborating sites.



