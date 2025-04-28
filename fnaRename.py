import json
import os
import glob

jsonl_file = "./ncbi_dataset/data/assembly_data_report.jsonl"

main_directory = "./ncbi_dataset/data"

# Possible genome file extensions
extensions = [".fna", ".fasta", ".fa"]

# Process each line in the JSONL file
with open(jsonl_file, "r") as f:
    for line in f:
        data = json.loads(line.strip())

        # Extract relevant info
        accession = data.get("accession", "unknown_accession")
        organism_name = data.get("organism", {}).get("organismName", "Unknown_Organism").replace(" ", "_")

        # Walk through all subdirectories
        file_found = False
        for subdir, _, _ in os.walk(main_directory):
            # Use glob to search for files that match the accession in each subdirectory
            matching_files = []
            for ext in extensions:
                matching_files.extend(glob.glob(os.path.join(subdir, f"{accession}*{ext}")))  # Search in subfolder

            if matching_files:
                current_filename = matching_files[0]  # Take the first match
                file_extension = os.path.splitext(current_filename)[1]  # Preserve original extension
                new_filename = os.path.join(subdir, f"{organism_name}-{accession}{file_extension}")

                # Rename the file
                os.rename(current_filename, new_filename)
                print(f"Renamed: {current_filename} → {new_filename}")
                file_found = True
                break  # Stop searching once the file is found and renamed

        if not file_found:
            print(f"File not found for: {accession} (Checked extensions: {', '.join(extensions)})")
