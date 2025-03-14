import pandas as pd
import csv

# Load the PR/SR file
pr_sr_file = "HG002_pr_sr_percentage.tsv"
pr_sr_df = pd.read_csv(pr_sr_file, sep="\t", dtype=str)

# Create a dictionary to store only the relevant columns (Chromosome, Manta_Start, Manta_End, PR_Percentage, SR_Percentage)
header_data = {}

# Store only the desired columns in the dictionary
header_data['Chromosome'] = pr_sr_df['Chromosome'].str.strip().tolist()  # Strip spaces
header_data['Manta_Start'] = pr_sr_df['Manta_Start'].str.strip().tolist()  # Strip spaces
header_data['Manta_End'] = pr_sr_df['Manta_End'].str.strip().tolist()  # Strip spaces
header_data['PR_Percentage'] = pr_sr_df['PR_Percentage'].str.strip().tolist()  # Strip spaces
header_data['SR_Percentage'] = pr_sr_df['SR_Percentage'].str.strip().tolist()  # Strip spaces

# Print the dictionary with the first 5 entries for each column
print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(f"{key}: {values[:5]}...")  # Display only the first 5 values for brevity

# Initialize a counter for the matches
match_count = 0

# List to store the updated rows for the benchmark file
updated_rows = []

# Open the benchmarking BED file and compare each row with the PR/SR file
with open("benchmarking_HG002_header.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    
    # Skip the header row
    header_row = next(reader)  # Skipping the first line, which is the header
    
    # Add new columns for PR_Percentage and SR_Percentage to the header row
    header_row.extend(["PR_Percentage", "SR_Percentage"])
    
    # Process the rows from the benchmarking file
    for row in reader:
        chrom = row[0].strip()  # Strip spaces from Chromosome column
        start = int(row[1].strip())  # Strip spaces and convert to int
        end = int(row[2].strip())  # Strip spaces and convert to int
        
        # Initialize PR and SR percentages as 0
        pr_percentage = 0.0
        sr_percentage = 0.0
        
        # Now, filter only the rows in the PR/SR file where Chromosome matches
        matching_rows = pr_sr_df[pr_sr_df['Chromosome'].str.strip() == chrom]
        
        # Compare the Manta_Start and Manta_End positions for the matching rows
        for _, pr_row in matching_rows.iterrows():
            manta_start = int(pr_row['Manta_Start'].strip())  # Strip spaces and convert to int
            manta_end = int(pr_row['Manta_End'].strip())  # Strip spaces and convert to int
            
            manta_binary = int(row[5].strip())  # Get the MANTA value from the benchmark file (column 5)
            
            # Only perform the comparison if the MANTA value is 1
            if manta_binary == 1:
                # Check if Manta_Start and Manta_End fall within the Benchmark's start and end range
                if manta_start >= start and manta_end <= end:
                    match_count += 1
                
                    # Get PR percentage
                    pr_percentage = float(pr_row['PR_Percentage'].strip())  # Convert to float

                    # Check if SR_Percentage is not NaN and convert it to float, otherwise set it to 0
                    if pd.isna(pr_row['SR_Percentage']):
                        sr_percentage = 0.0
                    else:
                        sr_percentage = float(pr_row['SR_Percentage'].strip())  # Convert to float
        
        # Add the PR and SR percentages to the row
        row_with_percentages = row + [pr_percentage, sr_percentage]
        
        # Append the updated row to the list
        updated_rows.append(row_with_percentages)

# Write the updated rows with PR and SR percentages to a new file
output_file = "lala.bed"
with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")
    
    # Write the updated header row
    writer.writerow(header_row)
    
    # Write the data rows with new PR and SR percentage columns
    writer.writerows(updated_rows)

print(f"\nTotal number of matches: {match_count}")
print(f"Updated BED file saved to: {output_file}")

