import pandas as pd
import csv

# Load PR/SR data from the TSV file into a DataFrame
pr_sr_file = "HG002_pe_sr_dellyonly.tsv"
pr_sr_df = pd.read_csv(pr_sr_file, sep="\t", dtype=str)  # Read the file as strings to handle any non-numeric values

# Create a dictionary to store stripped and cleaned header data
header_data = {
    'Chromosome': pr_sr_df['Chromosome'].str.strip().tolist(),  # Chromosome column
    'BED_START': pr_sr_df['BED_Start'].str.strip().tolist(),
    'BED_END': pr_sr_df['BED_End'].str.strip().tolist(),
    'DELLY_START': pr_sr_df['Delly_Start'].str.strip().tolist(),  # Delly Start positions
    'DELLY_END': pr_sr_df['Delly_End'].str.strip().tolist(),  # Delly End positions
    'PE_DELLY': pr_sr_df['PE_Delly'].str.strip().tolist(),  # Paired-End support values
    'SR_DELLY': pr_sr_df['SR_Delly'].str.strip().tolist(),  # Split-Read support values
    'DELLY_COVERAGE': pr_sr_df['Delly_Coverage'].str.strip().tolist(),  # Total coverage values
    'SVTYPE': pr_sr_df['SVTYPE'].str.strip().tolist()  # SV type column
}

# Print the first 5 entries of the header data for verification
print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(key + ": " + str(values[:5]) + "...")

match_count = 0  # Counter for the number of matching rows found
updated_rows = []  # List to store rows with added percentages and coverage data

# Open the benchmarking BED file for reading
with open("benchmarking_HG002_header_smoove_manta.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")  # Read the file as a tab-delimited file
    
    header_row = next(reader)  # Read the header row
    
    # Extend the header with new columns for PE, SR, Delly Coverage, and SVTYPE
    header_row.extend(["PE_DELLY", "SR_DELLY", "DELLY_COVERAGE"])
    
    for row in reader:
        # Ensure each row has enough columns; fill missing ones with empty strings
        row = row + [''] * (len(header_row) - len(row) - 3)
        
        # Extract Chromosome, Start, and End positions from the current row
        chrom = row[0].strip() if len(row) > 0 else ""  # Chromosome
        start = int(row[1].strip()) if len(row) > 1 else 0  # Start position
        end = int(row[2].strip()) if len(row) > 2 else 0  # End position
        sv_type = row[3].strip() if len(row) > 3 else ""  # Benchmark SVTYPE value
        delly_binary = int(row[7].strip()) if len(row) > 7 else 0  # Delly binary value (1 indicates active)

        # Initialize default values for PR, SR, coverage, and SVTYPE
        pr_delly = '.'
        sr_delly = '.'
        delly_coverage = '.'
        
        if delly_binary == 1:
            # Filter PR/SR data for rows matching the Chromosome value
            matching_rows = pr_sr_df[pr_sr_df['Chromosome'].str.strip() == chrom]
        
            # Iterate through matching rows and compare start and end positions
            for _, pr_row in matching_rows.iterrows():
                pr_delly_start = int(pr_row['Delly_Start'].strip())
                pr_delly_end = int(pr_row['Delly_End'].strip())
                bed_start = int(pr_row['BED_Start'].strip())
                bed_end = int(pr_row['BED_End'].strip())
                
                if (start - 100 <= pr_delly_start <= end) or (start <= pr_delly_end <= end + 100):
                    match_count += 1
                    pr_delly = pr_row['PE_Delly'].strip()  # PE_Delly
                    sr_delly = pr_row['SR_Delly'].strip() if pd.notna(pr_row['SR_Delly']) else 0.0
                    delly_coverage = pr_row['Delly_Coverage'].strip()  # Delly_Coverage
                    break 
                
                elif sv_type == pr_row['SVTYPE'].strip():
                    if (start == bed_start) and (end == bed_end):
                        pr_delly = pr_row['PE_Delly'].strip()  # PE_Delly
                        sr_delly = pr_row['SR_Delly'].strip() if pd.notna(pr_row['SR_Delly']) else 0.0
                        delly_coverage = pr_row['Delly_Coverage'].strip()  # Delly_Coverage
                        break 

        # Add the PR, SR, Delly Coverage, and SVTYPE values to the current row
        row_with_percentages = row + [pr_delly, sr_delly, delly_coverage]
        updated_rows.append(row_with_percentages)  # Add updated row to the list

# Write the updated data to a new BED file
output_file = "benchmarking_HG002_header_smoove_manta_delly.bed"
with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")  # Create a writer object
    
    writer.writerow(header_row)  # Write the updated header
    writer.writerows(updated_rows)  # Write all updated rows

print("Total number of matches:", match_count)  # Print the total number of matches
print("Updated BED file saved to:", output_file)  # Print the output file name

# Verify that all rows in the output file align correctly with the header
print("\nVerifying column alignment:")
with open(output_file, "r") as file:
    reader = csv.reader(file, delimiter="\t")  # Read the output file
    headers = next(reader)  # Read the header row
    print("Headers:", headers)  # Print the header for verification
    for row in reader:
        if len(row) != len(headers):  # Check if the number of columns in a row matches the header
            print("Column mismatch found in row:", row)  # Print any mismatched rows
            break
    else:
        print("All rows are correctly aligned with the header.")  # Print success message if all rows match

