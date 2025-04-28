awk -f calcScoreHGT.awk \
    ./panaroo_results/gene_presence_absence_roary.csv \
    ./HGT_evidence/mgeProximity/panaroo_gene_map.tsv \
    ./HGT_evidence/gcDeviation/gc_outliers.tsv \
    ./HGT_evidence/mgeProximity/mge_overlaps.tsv \
    ./HGT_evidence/blast_results/non_host_hits.tsv
