import pandas as pd
import csv

# Input and Output Files
pr_sr_file = "HG002_pr_sr_percentage.tsv"
pr_sr_df = pd.read_csv(pr_sr_file, sep="\t", dtype=str)

header_data = {
    'CHROMOSOME': pr_sr_df['Chromosome'].str.strip().tolist(),
    'BED_START': pr_sr_df['BED_Start'].str.strip().tolist(),
    'BED_END': pr_sr_df['BED_End'].str.strip().tolist(),
    'MANTA_START': pr_sr_df['Manta_Start'].str.strip().tolist(),
    'MANTA_END': pr_sr_df['Manta_End'].str.strip().tolist(),
    'PR_PERCENTAGE': pr_sr_df['PR_Percentage'].str.strip().tolist(),
    'SR_PERCENTAGE': pr_sr_df['SR_Percentage'].str.strip().tolist(),
    'MANTA_COVERAGE': pr_sr_df['Manta_Coverage'].str.strip().tolist(),
    'SVTYPE': pr_sr_df['SVTYPE'].str.strip().tolist()  # Include SVTYPE
}

print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(key + ": " + str(values[:5]) + "...")

match_count = 0
updated_rows = []

# Input and Output BED Files
with open("benchmarking_HG002_header_smoove.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    
    header_row = next(reader)  
    header_row.extend(["PR_PERCENTAGE", "SR_PERCENTAGE", "MANTA_COVERAGE"])  # Add new columns to the output header
    
    for row in reader:
        # Fill missing columns with empty strings
        row = row + [''] * (len(header_row) - len(row) - 3)
        
        chrom = row[0].strip() if len(row) > 0 else ""  # Chromosome value
        start = int(row[1].strip()) if len(row) > 1 else 0  # Start position
        end = int(row[2].strip()) if len(row) > 2 else 0  # End position
        #sv_type = row[3].strip() if len(row) > 3 else ""  # Benchmark SVTYPE value
        manta_binary = int(row[6].strip()) if len(row) > 6 else 0  # Manta binary value (usually column 5)
        
        pr_percentage = '.'
        sr_percentage = '.'
        manta_coverage = '.'

        # Only consider rows where MANTA == 1
        if manta_binary == 1:
            # Filter matching rows from the PR/SR data based on the Chromosome
            matching_rows = pr_sr_df[pr_sr_df['Chromosome'].str.strip() == chrom]
            
            # Iterate over matching rows and compare coordinates
            for _, pr_row in matching_rows.iterrows():
                pr_manta_start = int(pr_row['Manta_Start'].strip())
                pr_manta_end = int(pr_row['Manta_End'].strip())
                
                # Fallback: Match based on Manta coordinates
                if (start - 100 <= pr_manta_start <= end) or (start <= pr_manta_end <= end + 100):
                    match_count += 1
                    pr_percentage = float(pr_row['PR_Percentage'].strip())  # PR_Percentage
                    sr_percentage = float(pr_row['SR_Percentage'].strip()) if pd.notna(pr_row['SR_Percentage']) else 0.0
                    manta_coverage = float(pr_row['Manta_Coverage'].strip())  # Manta_Coverage value
                    break  

        row_with_percentages = row + [pr_percentage, sr_percentage, manta_coverage]
        updated_rows.append(row_with_percentages)

# Write Updated BED File
output_file = "benchmarking_HG002_header_smoove_manta.bed"
with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")
    
    writer.writerow(header_row)
    writer.writerows(updated_rows)

print("Total number of matches:", match_count)
print("Updated BED file saved to:", output_file)

# Verify Column Alignment
print("\nVerifying column alignment:")
with open(output_file, "r") as file:
    reader = csv.reader(file, delimiter="\t")
    headers = next(reader)
    print("Headers:", headers)
    for row in reader:
        if len(row) != len(headers):
            print("Column mismatch found in row:", row)
            break
    else:
        print("All rows are correctly aligned with the header.")

