read -p "Please enter the name of the taxon you want to generate a pangenome for: " taxon_name

./datasets download genome taxon "$taxon_name" --dehydrated --assembly-source RefSeq --assembly-level complete --include genome,seq-report --filename genomes.zip
unzip -qq genomes.zip
./datasets rehydrate --directory ./ 


echo "Sorting FNA files..."
mkdir -p fna_data
for file in ncbi_dataset/data/*/*.fna
do
    mv "${file}" "fna_data/"
done
#remove temp files
rm README.md md5sum.txt