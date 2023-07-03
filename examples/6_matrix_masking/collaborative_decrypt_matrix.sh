#!/bin/bash

if [[ $# != 2 ]]
then
	echo "USAGE: $0 [Encrypted matrix file] [# sites]"
	exit 1
fi

enc_matrix_file=$1
n_sites=$2

n_sites_min_one=`echo ${n_sites} | awk {'print $1-1}'`
site_iters=`seq 0 ${n_sites_min_one}`

rm -f *.partdec
for decrypting_i_site in ${site_iters[@]}
do
	echo "Partially decrypting ${i_site}.'s data matrix by site-${decrypting_i_site}.."
	smdg_noise_var_bit_size=0
	./COLLAGENE.sh -partial_decrypt_matrix SITE_${decrypting_i_site} ${enc_matrix_file} ${decrypting_i_site} ${smdg_noise_var_bit_size} ${enc_matrix_file}_partdec_by_${decrypting_i_site}.partdec
done

i_site=0

# Final step: Pool the matrix at all sites.
echo "Pooling partial decryption for site-${i_site}'s partial decryptions.."
ls *.partdec > PARTDECS.list
n_partdecs=`wc -l PARTDECS.list | awk {'print $1'}`
echo "Pooling ${n_partdec} parts.."
./COLLAGENE.sh -pool_partially_decrypted_matrices SITE_${i_site} PARTDECS.list ${enc_matrix_file}_collaborative_dec.bin

echo "Converting decrypting matrix to text.."
./COLLAGENE.sh -save_matrix_text ${enc_matrix_file}_collaborative_dec.bin ${enc_matrix_file}_collaborative_dec.bin.txt

