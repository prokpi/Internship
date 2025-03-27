import csv

input_file = "HG002_pr_sr_mantaonly.tsv"
output_file = "HG002_pr_sr_percentage.tsv"

def compute_percentage(partial, total_sum):
    if total_sum == 0:
        return 0
    return round((partial / total_sum) * 100, 2)  

def clean_tuple(value):
    if value and value != "(None)" and value.strip(): 
        try:
            cleaned_value = value.strip("()").replace(" ", "") 
            cleaned_value = cleaned_value.rstrip(')').lstrip('(')
            return tuple(map(int, cleaned_value.split(",")))  
        except ValueError as e:
            print("Error parsing value:", value, "-", e)
            return ()  
    return ()  

def process_file(input_file, output_file):
    results = []

    with open(input_file, 'r') as file:
        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            chrom = row["Chromosome"]
            bed_start = row["BED_Start"]
            bed_end = row["BED_End"]
            manta_start = row["Manta_Start"]
            manta_end = row["Manta_End"]
            svtype = row.get("SVTYPE", "Unknown")  # Include SVTYPE column

            # Extracting and cleaning PR (Paired-Read) data
            pr_tuple = clean_tuple(row["PR"])
            if pr_tuple:
                pr_percentage = compute_percentage(pr_tuple[1], sum(pr_tuple))
            else:
                pr_percentage = ""  

            # Extracting and cleaning SR (Split-Read) data
            sr_tuple = clean_tuple(row["SR"])
            if sr_tuple:
                sr_percentage = compute_percentage(sr_tuple[1], sum(sr_tuple))
            else:
                sr_percentage = ""  

            # Get Manta_Coverage directly from the input file (as it's already calculated)
            manta_coverage = row.get("Manta_Coverage", "")

            # Append all information to results, including SVTYPE
            results.append([chrom, bed_start, bed_end, manta_start, manta_end, svtype, pr_percentage, sr_percentage, manta_coverage])

    # Write to output file
    with open(output_file, 'w', newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["Chromosome", "BED_Start", "BED_End", "Manta_Start", "Manta_End", "SVTYPE", "PR_Percentage", "SR_Percentage", "Manta_Coverage"])
        writer.writerows(results)

    print("File processed and saved to", output_file)

# Call the function to process the file
process_file(input_file, output_file)


