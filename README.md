# Tempus Bioinformatics Challenge
**Will Peterson** <br />
Email: wcp7cp@virginia.edu <br />
LinkedIn: https://www.linkedin.com/in/willpeterson76/ <br />
Repo: https://github.com/your_username/repo_name <br />

## Instructions:
Tempus Bioinformatics Challenge: For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to annotate each variant in the file.

Each variant must be annotated with the following pieces of information:
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Using the VEP hgvs API, get the gene of the variant, type of variation (substitution,
insertion, CNV, etc.) and their effect (missense, silent, intergenic, etc.). The API
documentation is available here: https://rest.ensembl.org/#VEP
5. The minor allele frequency of the variant if available.
6. Any additional annotations that you feel might be relevant.

## Overview:
The variant annotation tool was developed and organized into several components. 
* The `utils/` directory stores all necessary functions used in annotating and writing results
  * `parser.py`: Function definitions for processing VCF files and annotating variants with additional information obtained from the Ensembl Variant Effect Predictor (VEP) API
  * `vcf_writer.py`: Functions used to write results to CSV and VCF files
* `requirements.txt` store all necessary dependencies for running the tool, these dependencies are also provided below:
  * `aiohttp==3.8.4`
  * `numpy==1.21.6`
  * `pandas==1.3.5`
  * `tqdm==4.64.0`
* `main.py` is the entry point for reading, extracting, and annotating inputted Variant Call Format (VCF) files.
  * Filename must be specified when executing main.py
  * Optional `--filter` tag may be included to only consider variants who pass the quality check
  * Example: `python main.py -f FILENAME [--filter]`
* `results/` contains the CSV and VCF files obtained from running the variant annotation tool with the provided test file ('data/test_vcf_data.txt')


## Results:
The depth of sequence coverage, number of reads supporting variant, and percentage of reads supporting the variant vs. reference could all be collected from the input VCF file. 
  1. **DepthCoverage**: Depth of sequence coverage was computed from the 'TC' INFO field.
  2. **NumVariantReads**: Number of reads supporting was found from the 'TR' INFO field. 
  3. **PercVariant**: Percentage of reads supporting the variant: TR / TC
  4. **PercReference**: Percentage of reads supporting reference: 1 - (TR / TC)
     
The remainder of required variant annotations were collected using the Ensembl VEP API: `https://grch37.rest.ensembl.org/vep/human/hgvs/`.

**Note:** that because the FASTA file (`/data/ref_genome/human_g1k_v37.fasta`) provided in `test_vcf_data.txt` referenced 'v37,' it was assumed the assembly was 'GRCh37.' Therefore, the Ensembl API Url was corrected to use the GRCh37 assembly in each request. With the Ensembl API, the following information was collected:
* **MAF**: Minor Allele Frequency
* **VariantID**: Variant ID
* **GeneSymbols**: Genes affected by the variant
* **Biotypes**: Gene and/or transcript classification
* **CQ**: All potential consequences arising from variant
* **Impact**:Degree of impact arising from variant
* **SevereCQ**: Most severe consequences arising from variant
* **AlleleString**: Shorthand allele notation
* **VariantClass**: Classification of variant

All results from running the variant annotation tool against the provided test file are stored under `results/`







