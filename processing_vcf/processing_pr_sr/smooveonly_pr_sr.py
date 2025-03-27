# before running the programme do these: 
# bgzip HG002_HG002-smoove.vcf.gz-filtered.vcf
# tabix -p vcf HG002_HG002-delly-recode.vcf-filtered.vcf.gz 
# for smoove VCF output, the fields PR and SR are specifically in INFO and FORMAT sections
# INFO fields:
    # PE - number of paired-end reads supporting the variant (PE=1 indicates 1 read supporting the variant)
    # SR - number of split-reads supporting the variant (e.g. SR=12 indicates 10 split-reads supporting the variant)
# FORMAT fields:
    # GT - Genotype
    # SU - total number of supporting reads (combined evidence from PE and SR)
    # PE - paired-end read count (evidence from discordant pairs)
    # SR - split-read count (evidence from split alignment)
# Interpretation would be: 1 paired-end read and 12 split-reads support this deletion. The total support is 13 SU

# PE (Paired-End Reads) and SR (Split-Reads) provide overall support counts.
# INFO    PE=1;SR=12
# FORMAT field Provides detailed support counts for individual samples
# e.g.: GT:SU:PE:SR  ./.:13:1:12; PE: 1 paired-end read; SR: 12 split-reads.


#bgzip
#bcftools sort HG002_HG002-smoove.vcf.filtered.vcf.gz -o HG002_HG002-smoove.sorted.vcf.gz
#bgzip -c HG002_HG002-smoove.sorted.vcf.gz > HG002_HG002-smoove.sorted.vcf.gz
#tabix -p vcf HG002_HG002-smoove.sorted.vcf.gz




import csv

delly_vcf_file = "HG002_HG002-smoove.vcf.gz-filtered.vcf"  # Smoove VCF file 
coordinates_bed_file = "HG002_merged_smooveonly.bed"  # BED file with regions of interest
output_pr_sv_tsv_file = "HG002_pe_sr_smooveonly.tsv"  # Output file

def extract_pr_sr(smoove_vcf, coordinates_bed, output_file, margin=100):
    results = []

    # Read BED file
    with open(coordinates_bed, 'r') as bed_file:
        bed_entries = [line.strip().split("\t") for line in bed_file]

    # Iterate over each bed entry
    for bed_entry in bed_entries:
        chrom, start, end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])

        # Define the region to search in the VCF with margin
        start_with_margin = max(0, start - margin)
        end_with_margin = end + margin

        # Open VCF file fresh for each BED region
        with open(smoove_vcf, 'r') as vcf_file:
            vcf_lines = vcf_file.readlines()

        # Process VCF file
        for line in vcf_lines:
            # Skip header lines
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            vcf_chrom = columns[0]
            vcf_start = int(columns[1])
            vcf_end = None
            svtype = None                               #
            info_field = columns[7]

            # Extract END value from INFO field and SVTYPE
            for info in info_field.split(';'):
                if info.startswith('END='):
                    vcf_end = int(info.split('=')[1])
                if info.startswith('SVTYPE='):          #
                    svtype = info.split('=')[1]         #

            # Ensure that both start and end of the variant are within the region (with margin)
            if vcf_chrom == chrom and (start_with_margin <= vcf_start <= end_with_margin or
                                        (vcf_end and start_with_margin <= vcf_end <= end_with_margin)):
                # If variant start and end are in the region of interest, process it
                smoove_start = vcf_start
                smoove_end = vcf_end if vcf_end else vcf_start  # Use start if no end is available

                # Extract sample data (assuming only one sample column)
                format_field = columns[8]
                sample_data = columns[9]
                format_values = sample_data.split(':')

                # Extract PE and SR from sample data
                pe = None
                sr = None
                if 'PE' in format_field:
                    pe_index = format_field.split(':').index('PE')
                    pe = format_values[pe_index] if len(format_values) > pe_index else None
                
                if 'SR' in format_field:
                    sr_index = format_field.split(':').index('SR')
                    sr = format_values[sr_index] if len(format_values) > sr_index else None
                
                # Format PE and SR as single values (they may be None)
                pe_str = str(pe) if pe else ""
                sr_str = str(sr) if sr else ""

                # Only calculate coverage once PE and SR are assigned
                pe_coverage = int(pe) if pe else 0
                sr_coverage = int(sr) if sr else 0
                total_coverage = pe_coverage + sr_coverage

                # Append results
                results.append([vcf_chrom, start, end, smoove_start, smoove_end, pe_str, sr_str, total_coverage, svtype])

    # Write the results to the output file
    with open(output_file, "w", newline="") as output_tsv:
        tsv_writer = csv.writer(output_tsv, delimiter="\t")
        tsv_writer.writerow(["Chromosome", "BED_Start", "BED_End", "Smoove_Start", "Smoove_End", "PE_Smoove", "SR_Smoove", "Smoove_Coverage", "SVTYPE"])

        for result in results:
            tsv_writer.writerow([str(result[0]).ljust(12), 
                             str(result[1]).ljust(12), 
                             str(result[2]).ljust(12),
                             str(result[3]).ljust(15), 
                             str(result[4]).ljust(15), 
                             str(result[5]).ljust(10),   # PE_Smoove 
                             str(result[6]).ljust(10),   # SR_Smoove 
                             str(result[7]).ljust(15),
                             str(result[8]).ljust(10)]) 

    print("Smoove PE and SR data with coordinates saved to", output_file)

def main():
    extract_pr_sr(delly_vcf_file, coordinates_bed_file, output_pr_sv_tsv_file)

if __name__ == "__main__":
    main()
