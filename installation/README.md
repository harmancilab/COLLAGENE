# Installation

## Dependencies
COLLAGENE runs on Linux systems. It requires working installations of *openssl* and *gcc* (version 9+).

Extensions require *git* and *curl*, which are available in most Linux installations.

## Steps to build COLLAGENE

### Source Code
You can clone COLLAGENE source code via git:
```
git clone https://github.com/harmancilab/COLLAGENE.git
```
After cloning, navigate into COLLAGENE directory:
```
cd COLLAGENE
```

### Install SEAL
It is first necessary to install SEAL that is located under *src/SEAL-main*. Note that COLLAGENE was initially developed with SEAL version 4.0.

Depending on the system, it may be necessary to install gcc-9 or above to ensure SEAL installation works with optimizations. If you install gcc manually, make sure to set the shared library paths correctly (LD_LIBRARY_PATH, and LIBRARY_PATH variables).

SEAL4.0 requires cmake with version greater than 3. If you don't have it installed on your system, you can download it using following:
```
cd src/SEAL-main
wget -c https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh
chmod 755 cmake-3.23.1-linux-x86_64.sh
./cmake-3.23.1-linux-x86_64.sh
```
Last command installs cmake in the local directory after confirming the license.

After cmake is installed, use following commands to setup cmake directories:
```
mkdir build
mkdir install
```

Finally, run cmake (assuming you installed it locally):
```
./cmake-3.23.1-linux-x86_64/bin/cmake -S . -B build -DSEAL_USE_CXX17=ON -DSEAL_BUILD_DEPS=ON -DSEAL_USE_MSGSL=ON -DSEAL_USE_ZSTD=ON -DSEAL_USE_ZLIB=OFF -DCMAKE_INSTALL_PREFIX=install -DSEAL_USE_INTEL_HEXL=ON
```

If gcc is locally installed, you can specify its location using following:
```
./cmake-3.23.1-linux-x86_64/bin/cmake -S . -B build -DSEAL_USE_CXX17=ON -DSEAL_BUILD_DEPS=ON -DSEAL_USE_MSGSL=ON -DSEAL_USE_ZSTD=ON -DSEAL_USE_ZLIB=OFF -DCMAKE_INSTALL_PREFIX=install -DSEAL_USE_INTEL_HEXL=ON -DCMAKE_C_COMPILER=gcc10-gcc -DCMAKE_CXX_COMPILER=gcc10-g++
```
Or you can set environment variables named *CXX* and *CC*:
```
export CXX=gcc10-c++
export CC=gcc10-cc
./cmake-3.23.1-linux-x86_64/bin/cmake -S . -B build -DSEAL_USE_CXX17=ON -DSEAL_BUILD_DEPS=ON -DSEAL_USE_MSGSL=ON -DSEAL_USE_ZSTD=ON -DSEAL_USE_ZLIB=OFF -DCMAKE_INSTALL_PREFIX=install -DSEAL_USE_INTEL_HEXL=ON
```

Finally, you can build and install SEAL:
```
./cmake-3.23.1-linux-x86_64/bin/cmake --build build
./cmake-3.23.1-linux-x86_64/bin/cmake --install build
```

cmake3.13 has a different syntax for installation (Second command)
```
cmake3 --build build --target install
```

These commands build and install SEAL locally under *install/* directory. 

### Install COLLAGENE
After installing SEAL, you can install COLLAGENE. 

At this stage, you need to make sure g++ version you used for building SEAL is the same as for building COLLAGENE. If you have a local version of gcc, open *Makefile* under base directory (above *src/*) and replace the line:
```
CXX=g++
```
with 
```
CXX=[path to local g++ executable] 
```

Finally, COLLAGENE can be built by running *make* at the main directory (above *src/*):
```
make clean
make
```

After running successfully, This should build *COLLAGENE_Release* executable under *bin/* directory.
You need to make sure that *bin/COLLAGENE_Release* is in your *PATH* variable so that it can be called by the COLLAGENE scripts under *scripts/* directory:

For the current bash session, you can set the PATH variable using:
```
export PATH=$PATH:${PWD}/bin
```

## Network Setup 

COLLAGENE's network I/O uses SCP and AWS transfers. More details on setting these up can be found under *examples/10_Network_File_IO/* example folder.

__Note__: Local folders can be used for testing, development, and debugging purposes without any requirement of configuration. This can be done by setting the local folder usage option in the configuration file for *FILE_IO_UTILS.sh* script.

## Docker Usage
In case your environment does not support installation of COLLAGENE, you can use the preset COLLAGENE image, in which you can setup the code.

To download and start docker, run following commands:
```
project_id="collagene_public"
image_tag="v3"

sudo docker pull secureomics/${project_id}:${image_tag}

sudo docker run --name COLLAGENE_client_app -v $PWD:/host -i -t secureomics/${project_id}:${image_tag} /bin/bash
```

This should start the docker container. You can install cmake and gcc-10 using following:
```
yum -y install cmake3
yum -y install gcc10.x86_64
yum -y install gcc10-c++
yum -y install bzip2
```

This version of cmake should install a version greater than 3.13 and can be used to install SEAL-4.0. 

After installing the SEAL library, you can build COLLAGENE executable using above commands under the main directory. Please make sure to set the g++ version in the Makefile that was used to build SEAL library.

Note that container has preinstalled AWS command line interface but it must be configured with user's credentials before it can be used. 

More information about configuring AWS client can be found here: https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html

Container should also have an R-4.2 installation and gcc10 installation (gcc10-c++).

# Running COLLAGENE
After COLLAGENE is installed, the best place to start is to go over the *Examples/* folder that contains tutorials starting from key generation and collaborative analysis examples.

You can use the executable (*bin/COLLAGENE_Release*) or the driver script (*scripts/COLLAGENE.sh*) to use COLLAGENE's functions. Driver script is a wrapper for the executable's functions.

__IMPORTANT__: If *bin/COLLAGENE_Release* executable is not visible in the directories listed in your *PATH* variable, you need to set following entry in *COLLAGENE.sh* script:
```
COLLAGENE_SECURE_EXEC=[Path to COLLAGENE_Release executable under bin/ directory]
```

If *COLLAGENE_Release* is not found by *COLLAGENE.sh*, it issues an error and exit.

 