import csv
from collections import Counter
import subprocess
import sys

def most_common_tsv(file_path, column_index, n=1):
    """
    Finds the nth most common string in a specified column of a TSV file.
    
    Args:
        file_path (str): Path to the TSV file
        column_index (int): Zero-based index of the column to analyze
        n (int): Rank of the most common string to retrieve (default is 1)
    """
    try:
        with open(file_path, 'r', newline='', encoding='utf-8') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
    
            # Skip header if your file has one (remove this line if no header)
            next(reader, None)
            
            # Extract the specified column
            column_values = [row[column_index] for row in reader if len(row) > column_index]
            
            if not column_values:
                print("No data found in the specified column.")
                return
            
            # Count occurrences and find most common
            counter = Counter(column_values)
            most_common_list = counter.most_common()
            
            if len(most_common_list) < n:
                print(f"Not enough items to find the {n}\U000000BA most common.")
                return
            
            nth_common = most_common_list[n-1]
            
            print(f"{n}\U000000BA most common string in column {column_index}: '{nth_common[0]}'")
            print(f"Occurrences: {nth_common[1]}")
            
    except Exception as e:
        print(f"An error occurred: {e}")



def extract_first_two_words(input_file, output_file, source_column):
    """
    Creates a new TSV with:
    - Column 1: Original value from column 1 of input TSV
    - Column 2: First two words from specified source column
    
    Args:
        input_file (str): Path to input TSV file
        output_file (str): Path to output TSV file
        source_column (int): Zero-based index of column to extract words from
    """
    try:
        with open(input_file, 'r', newline='', encoding='utf-8') as infile, \
             open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            
            # Write header for output file
            writer.writerow(['Original_ID', 'First_Two_Words'])
            
            # Skip header in input file if it exists
            next(reader, None)
            
            for row in reader:
                if len(row) > max(0, source_column):
                    original_id = row[0]
                    source_text = row[source_column]
                    
                    # Get first two words
                    words = source_text.split()[:2]
                    first_two_words = ' '.join(words)
                    
                    writer.writerow([original_id, first_two_words])
        
    except Exception as e:
        print(f"An error occurred: {e}")



with open("./HGT_evidence/summary.txt", "w", encoding="utf-8") as summary_file:
    original_stdout = sys.stdout
    sys.stdout = summary_file
    try:    
        most_common_tsv("./HGT_evidence/blast_results/non_host_hits_dedup.tsv", 1, 1)
    finally:
        sys.stdout = original_stdout

    summary_file.write("\n")
    sys.stdout = summary_file
    try:
        most_common_tsv("./HGT_evidence/blast_results/blast_results_species_dedup.tsv", 1, 1)
    finally:
        sys.stdout = original_stdout
