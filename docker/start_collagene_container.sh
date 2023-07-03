project_id="collagene_public"
image_tag="v3"

sudo docker pull secureomics/${project_id}:${image_tag}

sudo docker run --name COLLAGENE_client_app -v $PWD:/host -i -t secureomics/${project_id}:${image_tag} /bin/bash
