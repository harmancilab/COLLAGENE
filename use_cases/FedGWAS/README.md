This directory contains the instructions for running secure FedGWAS client.

In order to run the client, you need to first setup:
1) Download the data from zenodo, decompress, and navigate to the federated GWAS data directory:
```
wget -O data.tar.bz2 -c https://zenodo.org/record/8106630/files/CLLGN_DATA_DIR_07_02_23_20_00_40.tar.bz2?download=1
tar -xvjf data.tar.bz2 
cd CLLGN_DATA_DIR_07_02_23_20_00_40/Mega/PER_CLIENT_LOCAL_DIRS
```

2) Within the data directory, there are FedGWAS (Mega) and FedMeta (Meta) directories. 

In both of these, there are 3 directories corresponding to each site. Under, for example, *SITE_1/*, there are 3 more directories:

<ul>
<li> *LOCAL_DATA_DIR/*: Genotype, phenotype, covariates data.</li> 
<li> *SITE_KEYS/*: This directory holds the keys for this site, including the secret key share, public key etc. There should be preset keys in this folder for each site.</li> 
<li> *FEDERATION_KEYS/*: This directory holds the ssh keys, if the sites share encrypted datasets via a directory on an scp server.</li> 
</ul>

There are also numerous scripts that are used to run the protocols. One of the important files is *data_config.params*, which sets up the parameters. The parameters in this file should match between the sites, except for *LOCAL_SITE_I* parameter, which is the site's 0-based index. It is important to set this correctly so that protocol can keep track of files while aggregating and decrypting.

For federated GWAS, three scripts are used:
<ul>
<li> RUN_CLIENT_FEDGWAS_COMMANDS.sh : This script calls the client on the site, this is used to start FedGWAS client at each site.</li> 
<li> SECURE_FEDGWAS_CLIENT.sh : This is the wrapper for the utilities (13 steps of GWAS) </li> 
<li> SECURE_FEDGWAS_UTILITIES.sh : This is the script that implements each GWAS step. You can explore each file.</li> 
</ul>

Note that the pipeline uses a lot of intermediate files and file naming may get overwhelming since all data processing is done by passing matrices between different options rather than optimizing everything in memory.

3) Setup the keys from *KeyMaker*. Note that there should be preset keys for each site under, for example, *SITE_1/SITE_KEYS/*. These keys can be used for a test run.

If you want to create your own keys, you can refer to the key generation example under *examples/* folder of the source code, which explains step by step how to generate the keys. After you run this example, you just need to copy each set of keys to each site under *SITE_1/SITE_KEYS*, *SITE_2/SITE_KEYS*, *SITE_3/SITE_KEYS*.

4) Review the setting in data_config.params. Make sure to set COLLAGENE executable path in this configuration file unless it is included in *PATH*.

The default options should work without any modifications for a test run.

5) IO_TYPE should be set in *data_config.params* file and tested using FILE_IO_UTILITIES.sh script. You can use *FILE_IO_UTILS.sh* for this:

```
./FILE_IO_UTILITIES.sh -test_file_IO data_config.params
```

Note that by default, I/O type is set to a local directory, i.e., you can see the lines in *data_config.params* file:
```
IO_TYPE=LOCAL
LOCAL_REMOTE_SHARED_DIR=/data/SHARED_DATA_MODEL_DIR
```
Other options are commented out. To use AWS or SCP options, uncomment them after configuring AWS or SSH keys.

6) Validate the input files using; note that the pipeline does this check but one can explicitly validate data and some of the config options.
```
./SECURE_FEDGWAS_CLIENT.sh -validate_run
```

If there are no errors, this script should return without any messages.

7) Clean the shared working directory under the storage server before running any commands.

You can use following option for this:
```
./FILE_IO_UTILITIES.sh -clean_shared_directory data_config.params
```

8) You can now start the computations in the background using:
```
nohup ./RUN_CLIENT_FEDGWAS_COMMANDS.sh &
```

Note that this command should be run in the directory for each SITE.

When you are at *PER_CLIENT_LOCAL_DIRS/*, you can do this as:
```
cd SITE_1
nohup ./RUN_CLIENT_FEDGWAS_COMMANDS.sh &
cd ../SITE_2
nohup ./RUN_CLIENT_FEDGWAS_COMMANDS.sh &
cd ../SITE_3
nohup ./RUN_CLIENT_FEDGWAS_COMMANDS.sh &
```

These should have started 3 processes that run each client's commands. If you do not run or stop one of the sites, other sites will block at the current step of calculations.

9) To check progress, you can look for the time report that is used by each step in the time log file. Also, the output log should be written in each output file. These files are useful for troubleshooting errors.

10) If you need to stop all of the running protocol executables and scripts, you use following:
```
./SECURE_FEDGWAS_CLIENT.sh -stop_FedGWAS_processes
```
This option searches for COLLAGENE processes and tries to stop them. After this command ends, you should ensure all processes ended before starting the protocol again.

11) Final results are written to a tab-delimited file in each site's working directory, e.g., *FINAL_P_VALUES_ALL_VAR_BLOCKS_CLIENT_0_PER_R_w_var_ids.txt*