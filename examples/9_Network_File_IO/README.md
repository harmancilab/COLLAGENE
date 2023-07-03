# Network I/O (network_file_io_example.sh, FILE_IO_UTILITIES.sh)

COLLAGENE currently supports sharing of encrypted data files to share intermediate data among collaborating sites.

The encrypted intermediate data is stored in a remote server. COLLAGENE's *FILE_IO_UTILITIES.sh* script implements the options to download/upload/probe files in the remote server.

COLLAGENE supports 3 options for file I/O:
<ol>
<li> Local: Data is stored at a local folder. This option can be used for simulating and testing collaborative protocols by running all sites on the same client (Clients run in their folders and the shared directory is on the same computer.)</li> 
<li> SCP Server: Data is stored at an sFTP server that is access via SCP command with accession keys. COLLAGENE requires a key-based authentication to automate file transfers. This requires the sites to have authorized access to a shared ssh server. These keys can be setup using *ssh-keygen* and *ssh-copy-id* commands:</li> 

	# Run ssh-keygen to generate a public/private key pair for authentication, note that the passphrase should be skipped at this step.
	ssh-keygen -f test_key
	
	# Now, setup the key:
	ssh-copy-id -i test_key user_id@127.0.0.1
	
	# Replace "user_id" and "127.0.0.1" with username and server address.

<li> AWS S3 storage: Data is stored at an S3 bucket. It is necessary to setup the S3 bucket to transfer data to the bucket. The sites must have access to the bucket with write/read permissions. AWS CLI must be configured at each site (https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html).</li> 
</ol>

---

*network_file_io_example.sh*, *FILE_IO_UTILITIES.sh* contain more extensive examples and documentation of network file transfers that can be used for sharing data among collaborating sites.



