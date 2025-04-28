echo "Running BLAST..."

mkdir -p ./HGT_evidence/blast_results

blastn -query ./HGT_evidence/cloud_genome/cloudGeneSeqs.fa \
       -db ./nt_prok_db/nt_prok \
       -outfmt "6 qseqid stitle staxid pident length evalue bitscore sseq" \
       -num_threads 10 \
       -max_target_seqs 10 \
       > ./HGT_evidence/blast_results/blast_results.tsv

HOST_GENUS="Aeromonas"  # Set the host genus to filter out

awk -F'\t' -v host_genus="$HOST_GENUS" '
  {
    split($2, a, " ");          # Split stitle into words
    current_genus = a[1];       # Extract the first word (genus)
    if (tolower(current_genus) != tolower(host_genus)) {  # Case-insensitive comparison
      print $1 "\t" a[1] " " a[2];  # Print gene ID + full species name
    }
  }
' ./HGT_evidence/blast_results/blast_results.tsv > ./HGT_evidence/blast_results/non_host_hits.tsv

awk -F'\t' '
  {
    split($2, a, " ");          # Split stitle into words
      print $1 "\t" a[1] " " a[2];  # Print gene ID + full species name
  }
' ./HGT_evidence/blast_results/blast_results.tsv > ./HGT_evidence/blast_results/blast_results_species.tsv

sort -u ./HGT_evidence/blast_results/non_host_hits.tsv > ./HGT_evidence/blast_results/non_host_hits_dedup.tsv
sort -u ./HGT_evidence/blast_results/blast_results_species.tsv > ./HGT_evidence/blast_results/blast_results_species_dedup.tsv