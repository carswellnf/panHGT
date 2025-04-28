mkdir -p ./HGT_evidence
mkdir -p ./HGT_evidence/cloud_genome
# Extract unique genes
awk -F "," '$4 <= 15 {print $0}' ./panaroo_results/gene_presence_absence_roary.csv > ./HGT_evidence/cloud_genome/cloudGenes.list
awk -F ',' 'NR>15 {print $1}' ./HGT_evidence/cloud_genome/cloudGenes.list > ./HGT_evidence/cloud_genome/cloudGeneIds.list
seqkit grep -f ./HGT_evidence/cloud_genome/cloudGeneIds.list ./panaroo_results/pan_genome_reference.fa > ./HGT_evidence/cloud_genome/cloudGeneSeqs.fa