BEGIN {
    FS = OFS = "\t";
    print "Gene", "Annotation", "GC_outlier", "MGE_proximity", "Non_host_hits", "HGT_score" > "./HGT_evidence/gene_score_table.tsv";
    print "Gene", "Annotation", "Confidence" > "./HGT_evidence/all_hgt_candidates.tsv";
    print "Strong", "Moderate", "Weak", "Total"  > "./HGT_evidence/hgt_gene_counts.tsv"
}

# First load annotations from Panaroo file
FILENAME == ARGV[1] {
    if (FNR == 1) {
        # Skip header
        next;
    }
    
    # Manually parse CSV line to handle quoted fields with commas
    line = $0;
    gene = "";
    annotation = "";
    
    # Extract gene (first field)
    match(line, /^([^,]*),/, arr);
    gene = arr[1];
    line = substr(line, RLENGTH + 1);
    
    # Extract non-unique gene name (second field, discard)
    match(line, /^([^,]*),/, arr);
    line = substr(line, RLENGTH + 1);
    
    # Extract annotation (third field)
    if (substr(line, 1, 1) == "\"") {
        # Handle quoted annotation
        match(line, /^"([^"]*)"/, arr);
        annotation = arr[1];
    } else {
        # Handle unquoted annotation
        match(line, /^([^,]*)/, arr);
        annotation = arr[1];
    }
    
    annotations[gene] = annotation;
    next;
}

# Load all gene IDs from panaroo gene map
FILENAME == ARGV[2] && FNR > 1 {
    all_genes[$1] = 1;
    next;
}

# GC outliers
FILENAME == ARGV[3] {
    gc_outliers[$1] = 1;
    next;
}

# MGE proximity
FILENAME == ARGV[4] {
    mge_genes[$1] = 1;
    next;
}

# Non-host hits
FILENAME == ARGV[5] {
    non_host_counts[$1]++;
    next;
}

END {
    for (gene in all_genes) {
        gc_flag = (gene in gc_outliers) ? 1 : 0;
        mge_flag = (gene in mge_genes) ? 1 : 0;
        non_host_flag = 0;
        non_host_score = 0;
        
        if (gene in non_host_counts) {
            non_host_flag = 1;
            non_host_score = (non_host_counts[gene] >= 5) ? 40 : 30;
        }
        
        score = (gc_flag * 20) + (mge_flag * 40) + non_host_score;
        annotation = (gene in annotations) ? annotations[gene] : "NA";
        
        # Write to main priority table
        print gene, annotation, gc_flag, mge_flag, non_host_flag, score > "./HGT_evidence/gene_score_table.tsv";
        
        if (score == 100) {
            confidence = "strong"
            strong_count ++
        }
        if (score >= 80 && score < 100) {
            confidence = "moderate"
            moderate_count ++
        }
        if (score >= 60 && score < 80) {
            confidence = "weak"
            weak_count ++
        }

        # Write to candidates list if score is >= 70
        if (score >= 60) {
            print gene, annotation, confidence > "./HGT_evidence/all_hgt_candidates.tsv";
        }
    }
    print strong_count, moderate_count, weak_count, strong_count + moderate_count +  weak_count > "./HGT_evidence/hgt_gene_counts.tsv"
}