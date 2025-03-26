import pysam
import csv

survivor_vcf = "HG002_merged.vcf"
output_bed = "HG002_merged_mantaonly_filtered.bed"

def extract_coordinates(vcf_path, bed_path):
    with pysam.VariantFile(vcf_path) as vcf:
        with open(bed_path, "w", newline="") as bed_file:
            bed_writer = csv.writer(bed_file, delimiter='\t')
            
            for record in vcf:
                supp_vec = record.info.get("SUPP_VEC")
                
                # Check if SUPP_VEC matches any of the desired patterns
                if supp_vec in ['010', '110', '011', '111']:  
                    chrom = record.chrom           
                    start = record.pos - 1  # Convert VCF 1-based to BED 0-based
                    end = record.stop
                
                    # Write to the BED file
                    bed_writer.writerow([chrom, start, end])
                
                    print("Chromosome:", chrom, "Start:", start, "End:", end)

extract_coordinates(survivor_vcf, output_bed)
print("Filtered coordinates extracted and saved to:", output_bed)
