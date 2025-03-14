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
            manta_start = row["Manta_Start"]
            manta_end = row["Manta_End"]

            pr_tuple = clean_tuple(row["PR"])
            if pr_tuple:
                pr_percentage = compute_percentage(pr_tuple[1], sum(pr_tuple))
            else:
                pr_percentage = ""  

            sr_tuple = clean_tuple(row["SR"])
            if sr_tuple:
                sr_percentage = compute_percentage(sr_tuple[1], sum(sr_tuple))
            else:
                sr_percentage = ""  
            
            results.append([chrom, manta_start, manta_end, row["PR"], row["SR"], pr_percentage, sr_percentage])

    with open(output_file, 'w', newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["Chromosome", "Manta_Start", "Manta_End", "PR", "SR", "PR_Percentage", "SR_Percentage"])
        writer.writerows(results)

    print("File processed and saved to", output_file)

process_file(input_file, output_file)



