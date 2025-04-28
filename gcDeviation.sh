echo "Finding GC% outliers in cloud genome sequences..."

#Step 1: Calculate mean and standard deviation of GC% from pan_genome_reference.fa
mkdir -p ./HGT_evidence/gcDeviation
mkdir -p ./HGT_evidence/gcDeviation/gc_per_genome

# Loop through all FNA files and compute GC%
for fna in ./fna_data/*.fna; do
    genome_name=$(basename "$fna" .fna)
    grep -v ">" "$fna" | tr -d '\n' | awk '
        {len = length($0); gc = gsub(/[GC]/, ""); print (gc / len) * 100}
    ' > "./HGT_evidence/gcDeviation/gc_per_genome/${genome_name}_gc.txt"
done

# Combine all GC% values
cat ./HGT_evidence/gcDeviation/gc_per_genome/*_gc.txt > ./HGT_evidence/gcDeviation/combined_gc_values.txt.tmp

# Calculate mean and standard deviation
read mean stdev <<< $(awk '
    {sum += $1; sumsq += $1^2; total++} 
    END {
        mean = sum/total
        stdev = sqrt(sumsq/total - (sum/total)^2)
        print mean, stdev
    }' ./HGT_evidence/gcDeviation/combined_gc_values.txt.tmp)

# Print mean and standard deviation
echo "Mean GC%: $mean"
echo "Standard deviation: $stdev"

echo -e "Mean GC%: $mean\nStandard deviation: $stdev" > ./HGT_evidence/gcDeviation/gc_mean_stdev.txt

# Step 2: Filter outliers
awk -v mean="$mean" -v stdev="$stdev" '
    BEGIN {FS = "\t"; OFS = "\t"; print "Gene_ID", "GC%"}
    /^>/ {gene = substr($0,2); next}
    {gc = 0; len = length($0);
     for (i=1; i<=len; i++) {char = substr($0,i,1); 
     if (char ~ /[GCgc]/) gc++}
     gc_pct = (gc/len)*100;
     if (gc_pct < (mean - 3*stdev) || gc_pct > (mean + 3*stdev)) 
     print gene, gc_pct}
' ./HGT_evidence/cloud_genome/cloudGeneSeqs.fa > ./HGT_evidence/gcDeviation/gcOutliersUnfilterd.tsv.tmp

# Remove duplicate Gene_IDs from gc_outliers.tsv
awk '!seen[$1]++' ./HGT_evidence/gcDeviation/gcOutliersUnfilterd.tsv.tmp > ./HGT_evidence/gcDeviation/gc_outliers.tsv

# Remove temp files
rm ./HGT_evidence/gcDeviation/combined_gc_values.txt.tmp
rm ./HGT_evidence/gcDeviation/gcOutliersUnfilterd.tsv.tmp


echo "Results saved!"