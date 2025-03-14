import pandas as pd
import csv

pr_sr_file = "HG002_pr_sr_percentage.tsv"
pr_sr_df = pd.read_csv(pr_sr_file, sep="\t", dtype=str)

header_data = {}

header_data['Chromosome'] = pr_sr_df['Chromosome'].str.strip().tolist()  
header_data['Manta_Start'] = pr_sr_df['Manta_Start'].str.strip().tolist()  
header_data['Manta_End'] = pr_sr_df['Manta_End'].str.strip().tolist()  
header_data['PR_Percentage'] = pr_sr_df['PR_Percentage'].str.strip().tolist()  
header_data['SR_Percentage'] = pr_sr_df['SR_Percentage'].str.strip().tolist() 

print("Header Data Dictionary (first 5 entries):")
for key, values in header_data.items():
    print(key + ": " + str(values[:5]) + "...") 
    
match_count = 0

updated_rows = []

#Comparing each row with the PR/SR file
with open("benchmarking_HG002_header.bed", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    
    header_row = next(reader)  
    
    header_row.extend(["PR_Percentage", "SR_Percentage", "Manta_Start", "Manta_End"])
    
    #Processing rows in BED file
    for row in reader:
        chrom = row[0].strip()  
        start = int(row[1].strip()) 
        end = int(row[2].strip())
        
        pr_percentage = 0.0
        sr_percentage = 0.0
        manta_start = ""
        manta_end = ""
        
        #Filtering rows in the PR/SR file where Chromosome matches
        matching_rows = pr_sr_df[pr_sr_df['Chromosome'].str.strip() == chrom]
        
        #Comparing  Manta_Start and Manta_End positions for the matching rows:
        for _, pr_row in matching_rows.iterrows():
            pr_manta_start = int(pr_row['Manta_Start'].strip())  
            pr_manta_end = int(pr_row['Manta_End'].strip())  
            
            manta_binary = int(row[5].strip())  #Getting binary MANTA value
            
            #Comparison only if manta == 1
            if manta_binary == 1:
                # Checking if Manta_Start and Manta_End fall within the Benchmark's start and end range
                if pr_manta_start >= start and pr_manta_end <= end:
                    match_count += 1
                
                    pr_percentage = float(pr_row['PR_Percentage'].strip())  # Convert to float
                    sr_percentage = float(pr_row['SR_Percentage'].strip()) if pd.notna(pr_row['SR_Percentage']) else 0.0
                    
                    manta_start = pr_manta_start
                    manta_end = pr_manta_end
                    break  
        
        row_with_percentages = row + [pr_percentage, sr_percentage, manta_start, manta_end]
        
        updated_rows.append(row_with_percentages)

output_file = "benchmarking_updated_with_percentages_coverage.bed"
with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")
    
    writer.writerow(header_row)
    
    writer.writerows(updated_rows)

print("Total number of matches:", match_count)
print("Updated BED file saved to:", output_file)













