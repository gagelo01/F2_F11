#!/bin/bash
myArray=(1c_compute_wald.R 1e_vep.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 1c_compute_wald.R' && ./1c_compute_wald.R && echo 'Initializing 1e_vep.R' && ./1e_vep.R && echo 'The master script finished without errors'
