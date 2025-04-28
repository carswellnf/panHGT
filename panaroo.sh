# Make Pangenome using 10 CPU threads
echo "Running panaroo..."
panaroo -i gff3_filtered/*.gff -o panaroo_results \
    --threads 10 \
    -a core \
    --clean-mode sensitive \
    --remove-invalid-genes \

echo "Panaroo Complete!"
