import pandas as pd
import csv

pr_sr_file = "HG002_pr_sr_mantaonly_with_coverage.tsv"
pr_sr_df = pd.read_csv(pr_sr_file, sep="\t", dtype=str)

header_data = {
    'Chromosome': pr_sr_df['Chromosome'].str.strip().tolist(),
    'Manta_Start': pr_sr_df['Manta_Start'].str.strip().tolist(),
    'Manta_End': pr_sr_df['Manta_End'].str.strip().tolist(),
    'PR_Percentage': pr_sr_df['PR_Percentage'].str.strip().tolist(),
    'SR_Percentage': pr_sr_df['SR_Percentage'].str.strip().tolist(),
    'Manta_Coverage': pr_sr_df['Manta_Coverage'].str.strip().tolist()  
}

print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(key + ": " + str(values[:5]) + "...")

match_count = 0
updated_rows = []

#total number of sample calls = 9111; total no. of MANTA appearances = 9111
with open("benchmarking_HG002_header.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    
    header_row = next(reader)  
    
    header_row.extend(["PR_Percentage", "SR_Percentage", "Manta_Coverage"])
    
    for row in reader:
        # Filling missing columns with empty strings:
        row = row + [''] * (len(header_row) - len(row) - 3)  
        
        chrom = row[0].strip() if len(row) > 0 else ""  # Chromosome value
        start = int(row[1].strip()) if len(row) > 1 else 0  # Start position
        end = int(row[2].strip()) if len(row) > 2 else 0  # End position
        manta_binary = int(row[5].strip()) if len(row) > 5 else 0  # Manta binary value (usually column 5)

        pr_percentage = 0.0
        sr_percentage = 0.0
        manta_coverage = 0.0

        # Filter matching rows from the PR/SR data based on the Chromosome and Start position
        matching_rows = pr_sr_df[pr_sr_df['Chromosome'].str.strip() == chrom]
        
        # Iterate over matching rows and compare Manta_Start and Manta_End positions
        for _, pr_row in matching_rows.iterrows():
            pr_manta_start = int(pr_row['Manta_Start'].strip())
            pr_manta_end = int(pr_row['Manta_End'].strip())
            
            # Compare only if manta == 1 (i.e., Manta is active)
            if manta_binary == 1:
                # Check if Manta_Start and Manta_End fall within the Benchmark's start and end range
                if pr_manta_start >= start and pr_manta_end <= end:
                    match_count += 1

                    # Assign values from the PR/SR data
                    pr_percentage = float(pr_row['PR_Percentage'].strip())  # PR_Percentage
                    sr_percentage = float(pr_row['SR_Percentage'].strip()) if pd.notna(pr_row['SR_Percentage']) else 0.0
                    manta_coverage = float(pr_row['Manta_Coverage'].strip())  # Manta_Coverage value
                    break  

        row_with_percentages = row + [pr_percentage, sr_percentage, manta_coverage]
        
        updated_rows.append(row_with_percentages)

output_file = "benchmarking_updated_with_percentages_coverage.bed"
with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")
    
    writer.writerow(header_row)
    
    writer.writerows(updated_rows)

print("Total number of matches:", match_count)
print("Updated BED file saved to:", output_file)

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












