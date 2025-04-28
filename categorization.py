import pandas as pd
import re
import json

# 1. Load data
df = pd.read_csv("./HGT_evidence/all_hgt_candidates.tsv", sep="\t")

# 2. Load rules from JSON file
with open("./gene_category_dictionary.json", "r") as file:
    rules = json.load(file)

# 3. Categorize function
def categorize(annotation):
    if pd.isna(annotation):
        return "Other"
    
    categories = []
    annotation_lower = str(annotation).lower()
    
    for category, keywords in rules.items():
        for keyword in keywords:
            # Use regex to match keywords as substrings, ignoring case
            if re.search(rf"\b{keyword.lower()}\b", annotation_lower, re.IGNORECASE):
                categories.append(category)
                break
                
    return ", ".join(categories)

# 4. Apply categorization
df["Functional_Category"] = df["Annotation"].apply(categorize)
df.to_csv("./HGT_evidence/categorized_genes.tsv", sep="\t", index=False)