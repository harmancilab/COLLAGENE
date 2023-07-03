# Docker Usage
This directory contains the script for pulling the COLLAGENE public client image from DockerHub and starting a container.

Image can be pulled from DockerHub using:
```
project_id="collagene_public"
image_tag="v3"

sudo docker pull secureomics/${project_id}:${image_tag}
```

You can start the container using:
```
sudo docker run --name COLLAGENE_client_app -v $PWD:/host -i -t secureomics/${project_id}:${image_tag} /bin/bash
```

The container should have an installation of gcc10, aws (not configured), and R4.2, and openssl which are needed for running COLLAGENE use cases and examples.

You can checkout the code using git and build COLLAGENE:
```
git clone https://github.com/harmancilab/COLLAGENE.git
```