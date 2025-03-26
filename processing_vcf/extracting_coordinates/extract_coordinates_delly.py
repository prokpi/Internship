import pysam 
import csv 

survivor_vcf = "HG002_merged.vcf"
output_bed = "HG002_merged_dellyonly.bed"


def extract_coordinates(vcf_path, bed_path):
    
    with pysam.VariantFile(vcf_path) as vcf:
        with open(bed_path, "w", newline ="") as bed_file:
            bed_writer = csv.writer(bed_file, delimiter = '\t') 
            
            for record in vcf:
                supp_vec = record.info.get("SUPP_VEC") 
                if supp_vec.endswith("1"): 
                    chrom = record.chrom           
                    start = record.pos - 1        
                    end = record.stop
                
                    bed_writer.writerow([chrom, start, end])
                
                    print("Chromosome: ", chrom, "Start: ", start, "End:", end)

extract_coordinates(survivor_vcf, output_bed)
