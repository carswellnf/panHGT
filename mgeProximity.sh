#!/bin/bash
set -euo pipefail

# ========================
# Configuration
# ========================
HGT_EVIDENCE_DIR="./HGT_evidence"
ANNOTATIONS_DIR="./gff3_filtered"
PANAROO_RESULTS="./panaroo_results"
OUTPUT_DIR="${HGT_EVIDENCE_DIR}/mgeProximity"
WINDOW_SIZE=5000
GENE_PRESENCE_FILE="${PANAROO_RESULTS}/gene_presence_absence_roary.csv"
UNIQUE_GENES_FILE="${HGT_EVIDENCE_DIR}/cloud_genome/cloudGeneSeqs.fa"
FIRST_ISOLATE_COL=15  # Column where isolate data starts

# ========================
# Initialization
# ========================
echo "Initializing MGE filtering pipeline..."
mkdir -p "$OUTPUT_DIR"
mkdir -p "${OUTPUT_DIR}/MGE_beds"

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check file existence and content
check_file() {
    if [ ! -f "$1" ]; then
        log "ERROR: File $1 not found!"
        exit 1
    fi
    if [ ! -s "$1" ]; then
        log "ERROR: File $1 is empty!"
        exit 1
    fi
}

# ========================
# Step 1: Extract outlier gene IDs
# ========================
log "Step 1: Extracting gene IDs..."
check_file "${UNIQUE_GENES_FILE}"
awk '/^>/ {print substr($0, 2)}' "${UNIQUE_GENES_FILE}" > "${OUTPUT_DIR}/gene_ids.list.tmp"
check_file "${OUTPUT_DIR}/gene_ids.list.tmp"
log "Extracted $(wc -l < "${OUTPUT_DIR}/gene_ids.list.tmp") genes"

# ========================
# Step 2: Create gene locations file
# ========================
log "Step 2: Creating gene locations file..."
for gff in "${ANNOTATIONS_DIR}"/*.gff; do
    check_file "$gff"
    awk -F'\t' '
        $3 == "gene" {
            gene_id = "";
            split($9, attributes, ";");
            for (i in attributes) {
                if (attributes[i] ~ /^ID=/) {
                    split(attributes[i], id_field, "=");
                    gene_id = id_field[2];
                    sub(/_gene$/, "", gene_id);  # Remove "_gene" suffix
                    break;
                }
            }
            if (gene_id != "") {
                print gene_id "\t" $1 "\t" $4 "\t" $5
            }
        }' "$gff"
done > "${OUTPUT_DIR}/gene_locations.tsv.tmp"
check_file "${OUTPUT_DIR}/gene_locations.tsv.tmp"
log "Gene locations file created with $(wc -l < "${OUTPUT_DIR}/gene_locations.tsv.tmp") entries"

# ========================
# Step 3: Create genome size file
# ========================
log "Step 3: Creating genome size file..."
for gff in "${ANNOTATIONS_DIR}"/*.gff; do
    genome=$(basename "$gff" .gff)
    awk -F'\t' '
        $1 !~ /^#/ && $5 ~ /^[0-9]+$/ {
            if ($5 > max[$1]) max[$1] = $5
        }
        END {
            for (c in max) print c "\t" max[c]
        }' "$gff" > "${OUTPUT_DIR}/MGE_beds/${genome}_sizes.tmp"
    check_file "${OUTPUT_DIR}/MGE_beds/${genome}_sizes.tmp"
done
cat "${OUTPUT_DIR}"/MGE_beds/*_sizes.tmp | sort -u > "${OUTPUT_DIR}/genome_sizes.txt.tmp"
check_file "${OUTPUT_DIR}/genome_sizes.txt.tmp"
rm "${OUTPUT_DIR}"/MGE_beds/*_sizes.tmp
log "Genome sizes file created with $(wc -l < "${OUTPUT_DIR}/genome_sizes.txt.tmp") contigs"

# ========================
# Step 4: Strict 1:1 Gene Mapping (First Isolate Only)
# ========================
log "Step 4: Running strict 1:1 gene mapping..."
ORIGINAL_GENES=$(wc -l < "${OUTPUT_DIR}/gene_ids.list.tmp")

awk -v FIRST_ISOLATE_COL="$FIRST_ISOLATE_COL" -v OUTPUT_DIR="$OUTPUT_DIR" '
    BEGIN {
        FS = ","
        OFS = "\t"
        total_mapped = 0
        
        # Load genes
        while ((getline < (OUTPUT_DIR "/gene_ids.list.tmp")) > 0) {
            outlier_genes[$1] = 1
        }
        close(OUTPUT_DIR "/gene_ids.list.tmp")
    }
    FNR == 1 {
        # Skip header
        next
    }
    $1 in outlier_genes {
        # Process only if exact match with cluster ID
        for (i = FIRST_ISOLATE_COL; i <= NF; i++) {
            if ($i != "" && $i != "NA") {
                # Take the first gene ID if multiple are present in cell
                split($i, isolate_genes, ";")
                gsub(/"/, "", isolate_genes[1])  # Remove quotes from first gene
                if (isolate_genes[1] != "") {
                    print $1, isolate_genes[1]
                    total_mapped++
                    break  # Exit after first found isolate
                }
            }
        }
    }
    END {
        print "Original outlier genes: " '"$ORIGINAL_GENES"' > (OUTPUT_DIR "/mapping_stats.log")
        print "Genes successfully mapped: " total_mapped > (OUTPUT_DIR "/mapping_stats.log")
    }' "$GENE_PRESENCE_FILE" > "${OUTPUT_DIR}/panaroo_gene_map.tsv"

check_file "${OUTPUT_DIR}/panaroo_gene_map.tsv"
MAPPED_GENES=$(cut -f1 "${OUTPUT_DIR}/panaroo_gene_map.tsv" | sort -u | wc -l)
log "Mapped $MAPPED_GENES genes out of $ORIGINAL_GENES"

rm "${OUTPUT_DIR}/gene_ids.list.tmp"

# ========================
# Step 5: Create genes.bed file
# ========================
log "Step 5: Creating genes.bed file..."
# Create a mapping from isolate gene ID to panaroo cluster ID
awk -F'\t' 'NR==FNR {cluster[$2]=$1; next} 
    {if ($1 in cluster) print $2 "\t" $3 "\t" $4 "\t" $1 "\t" cluster[$1]}' \
    "${OUTPUT_DIR}/panaroo_gene_map.tsv" \
    "${OUTPUT_DIR}/gene_locations.tsv.tmp" > "${OUTPUT_DIR}/genes.bed"

check_file "${OUTPUT_DIR}/genes.bed"
log "Created genes.bed with $(wc -l < "${OUTPUT_DIR}/genes.bed") entries"
rm "${OUTPUT_DIR}/gene_locations.tsv.tmp"

# ========================
# Step 6: Expand regions and find MGEs
# ========================
log "Step 6: Expanding regions around genes..."
bedtools slop -i <(cut -f1-4 "${OUTPUT_DIR}/genes.bed") -g "${OUTPUT_DIR}/genome_sizes.txt.tmp" -b $WINDOW_SIZE > "${OUTPUT_DIR}/expandedRegions.bed"
check_file "${OUTPUT_DIR}/expandedRegions.bed"
log "Expanded to $(wc -l < "${OUTPUT_DIR}/expandedRegions.bed") regions"
rm "${OUTPUT_DIR}/genome_sizes.txt.tmp"

log "Identifying MGEs in annotations..."
for gff in "${ANNOTATIONS_DIR}"/*.gff; do
    genome=$(basename "$gff" .gff)
    awk -F'\t' 'BEGIN{OFS="\t"} 
        tolower($0) ~ /mobile|transpos|phage|prophage|integrase|conjug|mobiliz|plasmid|virus|repfib|insertion|segrega|origin/{
            print $1,$4,$5,$3
        }' "$gff" > "${OUTPUT_DIR}/MGE_beds/${genome}_mges.bed"
    check_file "${OUTPUT_DIR}/MGE_beds/${genome}_mges.bed"
done

MGE_TMP="${OUTPUT_DIR}/mge_combined.bed.tmp"
cat "${OUTPUT_DIR}"/MGE_beds/*_mges.bed > "$MGE_TMP"
check_file "$MGE_TMP"
log "Found $(wc -l < "$MGE_TMP") total MGE features"

# ========================
# Step 7: Find overlaps between genes and MGEs
# ========================
log "Step 7: Finding gene-MGE overlaps..."
bedtools intersect -a "${OUTPUT_DIR}/expandedRegions.bed" -b "$MGE_TMP" -wo > "${OUTPUT_DIR}/gene_mge_overlaps.tsv.tmp"
check_file "${OUTPUT_DIR}/gene_mge_overlaps.tsv.tmp"
log "Found $(wc -l < "${OUTPUT_DIR}/gene_mge_overlaps.tsv.tmp") overlaps"
rm "${MGE_TMP}"
rm "${OUTPUT_DIR}/expandedRegions.bed"

# ========================
# Step 8: Create final output with cluster IDs
# ========================
log "Step 8: Creating final output TSV..."
# Join with panaroo_gene_map.tsv to get cluster IDs
awk -F'\t' 'NR==FNR {cluster[$2]=$1; next} 
    {if ($4 in cluster) print cluster[$4] "\t" $0}' \
    "${OUTPUT_DIR}/panaroo_gene_map.tsv" \
    "${OUTPUT_DIR}/gene_mge_overlaps.tsv.tmp" > "${OUTPUT_DIR}/mge_overlaps.tsv"
rm "${OUTPUT_DIR}/gene_mge_overlaps.tsv.tmp"
# Add header
echo -e "PanarooClusterID\tGeneChrom\tGeneStart\tGeneEnd\tGeneID\tMGEChrom\tMGEStart\tMGEEnd\tMGEType\tOverlapLength" | \
    cat - "${OUTPUT_DIR}/mge_overlaps.tsv" > "${OUTPUT_DIR}/mgeOverlapsHeader.tsv.tmp" && \
    mv "${OUTPUT_DIR}/mgeOverlapsHeader.tsv.tmp" "${OUTPUT_DIR}/mge_overlaps.tsv"

check_file "${OUTPUT_DIR}/mge_overlaps.tsv"
log "Final output created with $(($(wc -l < "${OUTPUT_DIR}/mge_overlaps.tsv") - 1)) gene-MGE associations"

# ========================
# Cleanup and completion
# ========================
log  "Results in ${OUTPUT_DIR}/mge_overlaps.tsv"