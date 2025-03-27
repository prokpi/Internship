import pysam
import csv
import os

# Input and output files
delly_vcf_file = "HG002_HG002-delly-recode.vcf-filtered.vcf.gz"  # Delly output VCF file
coordinates_bed_file = "HG002_merged_dellyonly.bed"  # Coordinates from BED file
output_pr_sv_tsv_file = "HG002_pe_sr_dellyonly.tsv"  # Output TSV file

def extract_pr_sr(delly_vcf, coordinates_bed, output_file, margin=100):
    results = []  # To collect Delly coordinates, PE, SR, and SVTYPE data

    # Read BED file
    with open(coordinates_bed, 'r') as bed_file:
        bed_entries = [line.strip().split("\t") for line in bed_file]

    # Open VCF file with pysam
    with pysam.VariantFile(delly_vcf) as vcf:
        for bed_entry in bed_entries:
            chrom, start, end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])

            # Adjust start and end with the margin
            start_with_margin = max(0, start - margin)
            end_with_margin = end + margin

            # Fetch records from the VCF that overlap with the coordinates
            for record in vcf.fetch(chrom, start_with_margin, end_with_margin):
                delly_start = record.pos
                delly_end = record.stop

                # Extract PE, SR, and SVTYPE from the INFO field
                pe = record.info.get("PE", None)  # Paired-end support
                sr = record.info.get("SR", None)  # Split-read support
                svtype = record.info.get("SVTYPE", "Unknown")  # Structural variant type, default "Unknown"

                # Calculate coverage
                pe_coverage = pe if pe else 0
                sr_coverage = sr if sr else 0
                total_coverage = pe_coverage + sr_coverage

                # Format PE, SR, and SVTYPE as strings
                pe_str = str(pe) if pe is not None else ""
                sr_str = str(sr) if sr is not None else ""
                svtype_str = str(svtype)

                # Add results
                results.append([chrom, start, end, delly_start, delly_end, pe_str, sr_str, total_coverage, svtype_str])
                
    # Write results to output file
    with open(output_file, "w", newline="") as output_tsv:
        tsv_writer = csv.writer(output_tsv, delimiter="\t")
        tsv_writer.writerow(["Chromosome", "BED_Start", "BED_End", "Delly_Start", "Delly_End", "PE_Delly", "SR_Delly", "Delly_Coverage", "SVTYPE"])  # Header

        for result in results:
            tsv_writer.writerow([
                str(result[0]).ljust(12),
                str(result[1]).ljust(12),
                str(result[2]).ljust(12),
                str(result[3]).ljust(15),
                str(result[4]).ljust(15),
                str(result[5]).ljust(10),  # PE_Delly with ljust
                str(result[6]).ljust(10),  # SR_Delly with ljust
                str(result[7]).ljust(15),  # Delly_Coverage
                str(result[8]).ljust(10)   # SVTYPE
            ])

    print("Delly PE, SR, and SVTYPE data with coordinates saved to", output_file)

# Main function
def main():
    extract_pr_sr(delly_vcf_file, coordinates_bed_file, output_pr_sv_tsv_file)

if __name__ == "__main__":
    main()