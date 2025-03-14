import pysam
import csv

manta_vcf_file = "HG002.diploidSV.vcf.gz"  
coordinates_bed_file = "HG002_merged_mantaonly.bed"
output_pr_sv_tsv_file = "HG002_pr_sr_mantaonly.tsv"


def extract_pr_sr(manta_vcf, coordinates_bed, output_file, margin=0):
    results = []  

    with open(coordinates_bed, 'r') as bed_file:
        bed_entries = [line.strip().split("\t") for line in bed_file]

    with pysam.VariantFile(manta_vcf) as vcf:
        for bed_entry in bed_entries:
            chrom, start, end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])  
            
            start_with_margin = max(0, start - margin)
            end_with_margin = end + margin

            for record in vcf.fetch(chrom, start_with_margin, end_with_margin):
                manta_start = record.pos         
                manta_end = record.stop           
                
                pr = record.samples.values()[0].get("PR", None)
                sr = record.samples.values()[0].get("SR", None)
                
                pr_str = "(" + str(pr[0]) + ", " + str(pr[1]) + ")" if pr else ""
                sr_str = "(" + str(sr[0]) + ", " + str(sr[1]) + ")" if sr else ""
                
                results.append([chrom, manta_start, manta_end, pr_str, sr_str])

    with open(output_file, "w", newline="") as output_tsv:
        tsv_writer = csv.writer(output_tsv, delimiter="\t")
        tsv_writer.writerow(["Chromosome", "Manta_Start", "Manta_End", "PR", "SR"])  
        
        for result in results:
            tsv_writer.writerow([str(result[0]).ljust(12), str(result[1]).ljust(12), str(result[2]).ljust(12),
                     str(result[3]).ljust(15), str(result[4]).ljust(15)])

    
    print("PR and SR data with coordinates saved to", output_file)

def main():
    extract_pr_sr(manta_vcf_file, coordinates_bed_file, output_pr_sv_tsv_file)

if __name__ == "__main__":
    main()
