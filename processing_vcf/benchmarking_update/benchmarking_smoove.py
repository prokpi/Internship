import pandas as pd
import csv

# Input file for Smoove data
pe_sr_file = "HG002_pe_sr_smooveonly.tsv"
pe_sr_df = pd.read_csv(pe_sr_file, sep="\t", dtype=str)  # Read the file as strings to handle any non-numeric values

# Create a dictionary to store stripped and cleaned header data
header_data = {
    'Chromosome': pe_sr_df['Chromosome'].str.strip().tolist(),  # Chromosome column
    'SMOOVE_START': pe_sr_df['Smoove_Start'].str.strip().tolist(),  # Smoove Start positions
    'SMOOVE_END': pe_sr_df['Smoove_End'].str.strip().tolist(),  # Smoove End positions
    'SVTYPE': pe_sr_df['SVTYPE'].str.strip().tolist() if 'SVTYPE' in pe_sr_df.columns else [''] * len(pe_sr_df),  # SVTYPE column (if present)
    'PE_SMOOVE': pe_sr_df['PE_Smoove'].str.strip().tolist(),  # Paired-End support values
    'SR_SMOOVE': pe_sr_df['SR_Smoove'].str.strip().tolist(),  # Split-Read support values
    'SMOOVE_COVERAGE': pe_sr_df['Smoove_Coverage'].str.strip().tolist()  # Smoove Coverage values
}

print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(key + ": " + str(values[:5]) + "...")

match_count = 0  # Counter for the number of matching rows found
updated_rows = []  # List to store rows with added percentages, coverage data, and SVTYPE

# Open the benchmarking BED file for reading
with open("benchmarking_HG002_header.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")  # Read the file as a tab-delimited file
    
    header_row = next(reader)  # Read the header row
    
    # Extend the header with new columns for PE, SR, Coverage, and SVTYPE
    header_row.extend(["PE_SMOOVE", "SR_SMOOVE", "SMOOVE_COVERAGE"])
    
    for row in reader:
        # Ensure each row has enough columns; fill missing ones with empty strings
        row = row + [''] * (len(header_row) - len(row) - 3)
        
        # Extract Chromosome, Start, and End positions from the current row
        chrom = row[0].strip() if len(row) > 0 else ""  # Chromosome
        start = int(row[1].strip()) if len(row) > 1 else 0  # Start position
        end = int(row[2].strip()) if len(row) > 2 else 0  # End position
        sv_type = row[3].strip() if len(row) > 3 else "" 
        smoove_binary = int(row[5].strip()) if len(row) > 5 else 0  # Smoove binary value (1 indicates active)

        # Initialize default values for PR, SR, coverage, and SVTYPE
        pe_percentage = '.'
        sr_percentage = '.'
        smoove_coverage = '.'

        if smoove_binary == 1:
            # Filter PR/SR data for rows matching the Chromosome value
            matching_rows = pe_sr_df[pe_sr_df['Chromosome'].str.strip() == chrom]
            
            # Iterate through matching rows and compare start and end positions
            for _, pe_row in matching_rows.iterrows():
                smoove_start = int(pe_row['Smoove_Start'].strip())  # Extract Smoove Start position
                smoove_end = int(pe_row['Smoove_End'].strip())  # Extract Smoove End position
                bed_start = int(pe_row['BED_Start'].strip())
                bed_end = int(pe_row['BED_End'].strip())
            
      
                if (start - 100 <= smoove_start <= end) or (start <= smoove_end <= end + 100):
                    match_count += 1  # Increment match count

                    # Assign values from the matching PR/SR row
                    pe_percentage = float(pe_row['PE_Smoove'].strip()) if pe_row['PE_Smoove'].strip().isdigit() else 0.0
                    sr_percentage = float(pe_row['SR_Smoove'].strip()) if pe_row['SR_Smoove'].strip().isdigit() else 0.0
                    smoove_coverage = float(pe_row['Smoove_Coverage'].strip()) if pe_row['Smoove_Coverage'].strip().isdigit() else 0.0
                    break  # Break after the first match
                
                elif sv_type == pe_row['SVTYPE'].strip():
                    if (start == bed_start) and (end == bed_end):
                        pe_percentage = pe_row['PE_Smoove'].strip()  # PE_Delly
                        sr_percentage = pe_row['SR_Smoove'].strip() if pd.notna(pe_row['SR_Smoove']) else 0.0
                        delly_coverage = pe_row['Smoove_Coverage'].strip()  # Delly_Coverage
                        break 

        # Add the PR, SR, Coverage, and SVTYPE values to the current row
        row_with_percentages = row + [pe_percentage, sr_percentage, smoove_coverage]
        updated_rows.append(row_with_percentages)  # Add updated row to the list

# Write the updated data to a new BED file
output_file = "benchmarking_HG002_header_smoove.bed"
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
