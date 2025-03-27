# Internship

Directory processing_vcf contains scripts for:

extracting_coordinates/
1. Extracting coordinates: from HG002_merged.vcf in order to extract indipendent coordinates for SMOOVE MANTA and DELLY structural variant (SV) calls
   
processing_pr_sr/
2. Obtaining PR and SR values: usign extracted cordinates pr (pe for SMOOVE and DELLY) and sr values were extracted from individual vcf files of all 3 tools
3. Converting PR and SR to percentage: For manta PR and SR was converted to a percentange, while SMOOVE and DELLY have single values for PE and SR

benchmarking_update/
4. Updating BED file: benchmarking file contains extracted features for machine learning model training for structural variance (SV) detection.
