echo "Running Prokka..."

# Directory containing .fna files
fna_dir="./fna_data"
# Create a new directory for .gff3 files
mkdir -p ./prokka
# Iterate over .fna files, run Prokka
genome_total=$(ls -1 "${fna_dir}"/*.fna | wc -l)
genome_number=1
for file in "${fna_dir}"/*.fna
do
    echo "Annotating genome ${genome_number} of ${genome_total}..."
    filename=$(basename "${file}")
    filename_without_extension="${filename%.*}"
    $HOME/prokka/bin/prokka "${file}" --quiet --cpus 10 --centre X --compliant --outdir "./prokka/${filename_without_extension}"
    echo "Genome ${genome_number} of ${genome_total} annotated!"
    genome_number=$((genome_number+1))
done

#move .gff3 files from "prokka" to new "gff3_data" directory
mkdir -p "./gff3_data"
for file in prokka/*/*.gff
do
    directory_name=$(dirname $file)
    accession=$(basename $directory_name)
    new_file="${directory_name}/${accession}.gff"
    mv "${file}" "${new_file}"
    mv "${new_file}" "gff3_data/"
done
echo "Genome annotation complete!"