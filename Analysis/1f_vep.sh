var=$(pwd)
cd /mnt/sda/boujer01/Transcriptomique/ensembl-vep
export PERL5LIB=/home/boujer01/perl5/lib/perl5/x86_64-linux-gnu-thread-multi/
output_name="${var}/Data/Modified/VEPoutput.txt"
input_vcf="${var}/Data/Modified/VEPinput.vcf"
./vep -i $input_vcf --o $output_name --assembly GRCh37 --cache --symbol --sift b --tab --fields Uploaded_variation,SYMBOL,Consequence,SIFT --dir_cache ./ --offline --force_overwrite
