#########################################################################
# COLLAGENE utilizes file transfers to pass encrypted data between sites.
#
# This makes it easier to track network transfer without need to setup interprocess communication and synchronization.
# It also makes the federated/collaborative implementations more fault tolerant as protocols can be resumed in case of network failures.
#
# COLLAGENE's file I/O is implemented in FILE_IO_UTILITIES.sh and wraps SCP and AWS CLI system calls.
# This script uses several variables to set the storage option (local/SCP/S3), which are stored in data_config.params file.
#
# COLLAGENE's file I/O can use 3 options:
#########################################################################
#
# 1) Local -- This is for simulation purposes where all sites are running on the same machine and use the local file system to share encrypted files.
# This option is selected by setting: "IO_TYPE=LOCAL" in data_config.prams file
#
# Requires only a folder on the local computer:
# LOCAL_REMOTE_SHARED_DIR=/internal/aharmanci1/dir/SHARED_DATA_MODEL_DIR
#
# This folder should be existing otherwise the script will give error and exit.
#
#########################################################################
#
# 2) SCP -- SCP-based file transfers to a remote file server.
# This option is selected by setting: "IO_TYPE=SCP" in data_config.params file
# 
# Requires setting up the key-based authentication for each site to access the remote file server.
# You can use following commands to generate keys using ssh-keygen and set them up on the remote file server:
# 
# username=aharmanci1
# hostname=129.106.31.208
# mkdir FEDERATION_KEYS
# ssh-keygen -t rsa -f FEDERATION_KEYS/federation -P ""
# ssh-copy-id -i FEDERATION_KEYS/federation ${username}@${hostname}
#
# After these are run successfully, it is also necessary to setup the folder that will be used to store the encrypted datasets at the remote file server.
# For this, one site can login to the server and create a folder at a known path.
# The remote directory path is set by: "SCP_REMOTE_SHARED_DIR=/data/aharmanci1/FEDERATION_DATA" in data_config.params file
#
#########################################################################
#
# 3) AWS/S3 bucket -- AWS/S3 bucket is used to share encrypted data matrices among collaborating sites.
# It is strongly recommended to setup an IAM account with read/write access to the bucket that will store the encrypted data matrices.
# This option is selected by setting: "IO_TYPE=S3" in data_config.params file
# 
# Requires downloading and configuring the AWS CLI: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
# After setting AWS CLI path:
# AWS_CLI_PATH=/home/aharmanci1/binaries/aws
# Shared data bucket:
# SHARED_DATA_BUCKET=s3://secureomics
# The key to the directory that will be used as the prefix:
# S3_REMOTE_SHARED_DIR_KEY=${SHARED_DATA_BUCKET}/SHARED_FEDGWAS_DATA_MODEL_DIR
#########################################################################

if [[ ! -f "./COLLAGENE.sh" ]]
then
	echo "Could not find COLLAGENE script, copy it from source code scripts/ directory."
	exit 1
fi

if [[ ! -f "data_config.params" ]]
then
	echo "Could not find data_config.params script, copy it from source code scripts/ directory."
	exit 1
fi

if [[ ! -f "/FILE_IO_UTILITIES.sh" ]]
then
	echo "Could not find /FILE_IO_UTILITIES.sh script, copy it from source code scripts/ directory."
	exit 1
fi

# Make all scripts executable.
#dos2unix *.sh
chmod 755 *.sh

# After setting up the option to share files among sites, each site should run the test option that uploads a random matrix, probes it, and downloads it:
echo "TESTING FILE TRANSFER AND PROBING"
./FILE_IO_UTILITIES.sh -test_file_IO data_config.params

# Make sure the testing succeeded.
if [[ $? -ne 0 ]]
then
	echo "File I/O failed, check I/O type and options in data_config.params file"
	exit 1
fi

# After this command runs successfully, it writes file I/O status information under a folder named TEMP_FILE_IO/. This folder keeps track of all file I/O operations.

# Following commands tests uploading files and directories.

# Upload the matrices.
mat_i_list=`seq 1 10`
rm -f matrix_files.list
for mat_i in ${mat_i_list[@]}
do
	cur_mat_file=random_matrix_${RANDOM}.bin
	./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 ${cur_mat_file}

	echo ${cur_mat_file} >> matrix_files.list
done

echo "UPLOADING 10 MATRICES..."
./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params matrix_files.list

echo "DOWNLOADING MATRICES..."
rm -f -r downloaded_matrices
mkdir downloaded_matrices
./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params matrix_files.list downloaded_matrices

echo "UPLOADING MATRIX DIRECTORY..."
echo "downloaded_matrices" > matrix_files.list
./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params matrix_files.list

echo "DOWNLOADING MATRIX DIRECTORY..."
rm -f -r downloaded_matrices_per_directory
mkdir downloaded_matrices_per_directory
./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params matrix_files.list downloaded_matrices_per_directory


