# Key Generation 
This example shows an example of generating keys using *KeyMaker* server (https://www.secureomics.org/KeyMaker) in a collaborative setting of 3 sites using the *COLLAGENE.sh* script.

Note that sites are numbered in the key-generation setup. The sites need to be assigned unique and consecutive ids (starting from 0) and adhere to the usage of these id's.

There are 4 steps in this example and the script key_generation_example.sh implement these steps in separate options.

## Setting parameters file *ckks.params*
CKKS keys are generated specifically for the set of encryption parameters (polynomial modulus, coefficient modulus, etc) that are described in *ckks.params* file. 

It is therefore necessary to agree on setting the CKKS parameters prior to setting up the keys. This requires the sites to know the depth and ciphertext sizes of algorithms.

Although the default parameters should be sufficient for many analysis, it may still be necessary to tune them depending on the applications. Example 2 (inside folder named *2_parameter_selection/*) includes a more detailed description of how parameters should be selected.

The sites agree on the CKKS parameters and save in the file named *ckks.params*. The default parameters can be used as a starting point.

## DSK Encryption Key Generation
Before key generation, the sites generate the keys *DSK encryption keys*, which are used for encrypting the key shares of each site by the KeyMaker. 

*DSK encryption keys* are simple openssl public/private key pairs that are used to encrypt/decrypt the actual distributed secret keys (DSKs). *DSK encryption keys* are only used for securely retrieving the DSKs and are not used for anything else. Server also uses them to determine the number of sites.

After generating the DSK encryption keys, a pair of keys (public and private) are generated at the site. The public key is sent to KeyMaker and private key is kept confidential at the site.

DSK encryption keys are generated using "-generate_DSK_encryption_key" option of COLLAGENE.sh script. 

For example, site-0 would use following to generate the keys:
```
COLLAGENE.sh -generate_DSK_encryption_key 0
```
This command generates the file site_0.dsk_enc_public_key and site_0.dsk_enc_private_key files. Similarly, sites 1 and 2 would run locally:
```
COLLAGENE.sh -generate_DSK_encryption_key 1
COLLAGENE.sh -generate_DSK_encryption_key 2
```
Private DSK decryption keys are protected with a password that is asked three times at the command line.

__IMPORTANT:__ The names of the key files should not be modified because they are used by *KeyMaker* while encrypting the key shares. Also, *COLLAGENE.sh* script uses each site's index while decrypting the encrypted key shares from *KeyMaker*.

After DSK encryption keys are generated, each site will share their public key with the site that will initiate the key request. The sharing of the public key does not create any risks since it is only used for encryption.

Also, the DSK encryption keys are very small in size and can be shared using the shared folder.

After the public DSK encryption keys are retrieved from all sites, these are copied into a folder and archived in a *.tar* file:
```
DSK_ENCRYPTION_KEYS_DIR=CLIENT_DSK_ENCRYPTION_KEYS
mkdir ${DSK_ENCRYPTION_KEYS_DIR}
cp *.dsk_enc_public_key ${DSK_ENCRYPTION_KEYS_DIR}
cp ckks.params ${DSK_ENCRYPTION_KEYS_DIR}
tar -cvf ${DSK_ENCRYPTION_KEYS_DIR}.tar ${DSK_ENCRYPTION_KEYS_DIR}
```
Note that the ckks.params file must be uploaded with the public DSK encryption keys because the HE keys will be generated for the specific set of parameters in the ckks parameters file.

After the archive file is created, it is uploaded to KeyMaker web site and submitted as a key generation task.

## Using KeyMaker Portal
Any on of the sites open a web browser and navigates to *KeyMaker* website at https://www.secureomics.org/KeyMaker. 

### Step 1: Generate DSK Encryption Keys
Step 1 of the analysis is described above for generating the DSK encryption keys. It includes the steps to generate the keys that are described above. We do not need this step in this tutorial since we already completed this step.

After uploading is complete, a table entry should be visible under uploads table.

### Step 2: Upload DSK Encryption Keys
Scroll down to step 2 and upload the tar file named ${DSK_ENCRYPTION_KEYS_DIR}.tar. This tar archive contains only the public DSK encryption keys (and ckks parameters file), does not leak any information other than the number of sites.

### Step 3: Start DSK Generation Task
After uploading, move onto step 3, where key generation task is submitted. Click on the key generation task submit button. It asks to verify submission and starts the task. 

There may be several updates to task status.

The key generation should finish within 10 minutes. After finished, a link should appear in the tasks table. Clicking on this will save the generated keys.

### Final Step: Download the keys
The result is a new archive file that is located on AWS cloud and the pre-signed URL can be shared with all sites (Time limit for this link is 2 weeks from the time it is finished, after that it becomes invalid.). 

This file can also be downloaded programmatically. For this, right click on the download button and save link to clipboard. Next, use wget to download the DSK tar file similar to following:
```
wget -O DSK_KEYS_FROM_SERVER.tar "http://s3.amazonaws.com/secureomics/KeyMaker_Data/c44d4f5b6fdcdf9772f59c5acd639d05eb7534fe/SITE_KEYS.tar?AWSAccessKeyId=AKIA4RBEEGAGVWXEAA4H&Expires=1669768655&Signature=H0LXeKIu%2BgDuIwjvzzVMgQovUbA%3D"
```
(The long link should be replaced with the copied link.)

Each site can use this to download the DSK tar file.

Although tar file includes the DSK for all sites, each site's DSK is encrypted with the corresponding public key and can only be decrypted using the site's private key. Also, DSK archive contains other keys that are shared among sites (public, Galois, relinearization) that do not leak any information. 

## Setting up HE-Key Shares
After each site downloads the archive file, they decompress and setup their own keys:

Site-0 performs following:
```
tar -xvf DSK_KEYS_FROM_SERVER.tar --strip-components 2 -C RECEIVED_KEYS
COLLAGENE.sh -decrypt_site_DSK 0 RECEIVED_KEYS
```

This command decrypts the key share for site-0. DSK decryption key file, e.g., *site_0.dsk_enc_private_key*, must be located in the same directory. The password for decryption will be asked at the command line and this must match the password provided while DSK encryption key was generated.

The received keys directory contains several other common keys (relinearization, Galois, partdec encryption keys) that will be used for executing protocols. It also contains the encrypted key-shares of other sites but they cannot be decrypted by other 3rd sites (each site can only decrypt their own key share).

Final step is setting up the keys for COLLAGENE script:
```
mkdir SITE_0
COLLAGENE.sh -set_params ckks.params RECEIVED_KEYS/pooled.public_key RECEIVED_KEYS/pooled.relin_keys RECEIVED_KEYS/pooled.galois_keys RECEIVED_KEYS/site_0.secret_key SITE_0
```
After this step, the site is ready to perform secure computations. This step simply copies the specified keys to the key directory.

Note that the partdec keys (See below) should be copied manually if they will be used:
```
cp RECEIVED_KEYS/partdec_data_enc_hash.symmetric_key SITE_${i_site}
```

## Partial decrypted data (partdec) Encryption Keys
KeyMaker generates symmetric keys that the sites can use to encrypt partially decrypted data matrices. These files are named respectively for each site (e.g., site_0_partdec_data_enc_hash.symmetric_key.enc) and they are encrypted by the site's DSK encryption key by KeyMaker.

After decryption of these partdec keys, these keys are named to partdec_data_enc_hash.symmetric_key.enc.

These files are the same for all sites and should only be used to encrypt partially encrypted matrices so that a third party does not have cleartext access to partially decrypted data on the shared space, i.e., SSH server or AWS cloud.

## Testing the Keys
*key_generation_example.sh* contains an example of generating a random matrix, encrypting it, and collectively decrypting it using the key shares that are just set up at each site's directory.

This option is implemented by "-example_test_ckks_keys". It generates random matrices on each site, then encrypts them using the shared public key, then decrypts each matrix.


