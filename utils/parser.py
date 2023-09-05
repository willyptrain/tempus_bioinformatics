import pandas as pd
import re
from aiohttp import ClientSession
import numpy as np
import logging
import json


def create_df(fp, filter_variants=False):
    """ 
    Read the input file into a DataFrame; optionally consider only variants passing filter with filter_variants=True 
    
    :param fp: File path of input VCF
    :param filter_variants: If true, consider only variants passing quality check.
    :return: VCF data lines stored as Pandas DataFrame
    """

    print("\nReading", fp, '\n')

    # Read lines of VCF into a list
    raw_vcf = open(fp,'r').readlines()
    # Find index of row containing VCF header (i.e. line containing #CHROM, ...)
    header_idx = get_header(raw_vcf)

    # Using header index, collect and clean up header lines
    header = raw_vcf[header_idx].strip('#\n').split('\t')

    # Read all lines following header row into dataframe / table
    vcf = pd.read_csv(fp, comment='#', delim_whitespace=True, names=header, skiprows=header_idx)
    
    # Number of variants found in input VCF
    print(len(vcf), 'variants found.', '\n')
    
    # If true, consider only variants which pass the quality check.
    if filter_variants:
        old_len = len(vcf)
        vcf = vcf[vcf['FILTER'] == 'PASS']
        print('Filtering...',old_len - len(vcf), 'variants dropped.', '\n')
        
    return vcf


def concat_ensembl_results(vcf: pd.DataFrame, results: list):
    """ 
    Parse all responses from Ensembl REST API and concatenate into DataFrame
    
    :param vcf: DataFrame containing data lines from VCF input file
    :param results: Byte string responses from Ensembl request
    :return: VCF data lines with necessary info fields concattenated to input vcf
    """


    vep_dict = {
        'MAF':[], # Minor Allele Frequency
        'VariantID':[], # Variant ID
        'GeneSymbols':[], # Symbols of Genes Affected by Variant
        'Biotypes':[], # Gene and/or transcript classification
        'CQ':[], # All potential consequences arising from variant
        'Impact':[], # Degree of impact arising from variant
        'SevereCQ':[], # Most severe consequences arising from variant.
        'AlleleString':[], # Shorthand allele notatio
        'VariantClass':[] # Classification of variant
    }

    # Iterate over list of JSON-formatted responses
    for resp in results:
        # Parse byte string as JSON
        json_resp = json.loads(resp)[0]
        
        # Read Ensembl results and store all relevant data in 'vep_dict'
        parsed_info = parse_vep_dict(json_resp)
        for key in vep_dict.keys():
            vep_dict[key].append(parsed_info[key])

    # Create a DataFrame from the parsed Ensembl results and concatenate with the original VCF DataFrame
    vep_df = pd.DataFrame(vep_dict)
    vcf = pd.concat((vcf.reset_index(drop=True),vep_df.reset_index(drop=True)),axis=1) 
    
    # Drop helper columns no longer needed
    vcf = vcf.drop('alt_split',axis='columns')
    vcf = vcf.drop('variant_prefix',axis='columns')
    
    return vcf


def explode_multiallele(vcf: pd.DataFrame):
    """ 
    Expand multi-allelic variants into multiple rows, s.t. each row contains only two alleles
    
    :param vcf: DataFrame containing data lines from VCF input file
    :return: VCF data lines expanded, with each row containing strictly two alleles
    """

    # Insert helper columns to aid in expanding dataframe
    vcf.insert(5, 'alt_split', vcf.apply(lambda row: row['ALT'].split(','),axis=1))
    vcf.insert(6, 'multi_allelic',vcf.apply(lambda row: True if ',' in row['ALT'] else False,axis=1))
    
    # For all ALT alleles containing ',' (implying multi-allelic), expand into multiple rows for each ALT allele
    vcf_biallelic = vcf.explode('alt_split', ',')

    # For each allele pair (REF and ALT), predict class of variant by comparing reference vs. alternative string. 
    vcf_biallelic['variant_prefix'] = vcf_biallelic.apply(lambda x: get_variant(x['REF'], x['ALT']),axis=1)

    return vcf_biallelic


def get_variant_statistics(vcf: pd.DataFrame):
    """ 
    Collect depth of sequence coverage, and degree of variant coverage (#,% supporting reads) using information in VCF INFO column
    
    :param vcf: DataFrame containing data lines from VCF input file
    :return: VCF data lines with new information added under new columns
    """

    # Depth of sequence coverage at the site of variation, stored in 'TC' field: Total coverage at this locus
    vcf['DepthCoverage'] = vcf.apply(lambda row: float(get_info_for_key('TC', row['INFO'])),axis=1)
    # Number of reads supporting the variant, stored in 'TR' field: Total number of reads containing this variant
    vcf['NumVariantReads'] = vcf.apply(lambda row: float(get_info_for_key('TR', row['INFO'])),axis=1)
    # Percentage reads supporting variant = Reads supporting variant / Total coverage = TR / TC
    vcf['PercVariant'] = vcf.apply(lambda row: row['NumVariantReads'] / row['DepthCoverage'],axis=1)
    # Percentage reads supporting reference = 1 - (TR / TC)
    vcf['PercReference'] = 1 - vcf['PercVariant']
    
    return vcf


def get_header(raw_vcf):
    """ 
    Return index of line where header information is stored
    
    :param raw_vcf: list containing each line in VCF file
    :return: row number / index where header information is located
    """
    for i, line in enumerate(raw_vcf):
        # '#CHROM' is signature of line containing header information
        if '#CHROM' in line:
            return i
        
    # If no header is found, return None
    return None


def get_info_for_key(key, info_str):
    """ 
    Return value assigned to given info field, specified by key
    
    :param key: Info field to be read (i.e. 'TC')
    :param info_str: INFO item to be parsed (represented as 'key=value')
    :return: value collected from INFO column for given key
    """
    found_info = re.findall('{0}=(\w.*?);'.format(key), info_str)[0]
    # If multi-allelic, return first found instance
    if ',' in found_info:
        return float(found_info[0])
    return float(found_info)


def get_reverse_complement(seq):
    """ 
    Compute reverse complement string for given sequence

    :param seq: Sequence from which to compute reverse complement
    :return: reverse complement string of provided sequence
    """
    
    # Dictionary mapping for each base to its complement
    reverse_comp = {'A':'T','C':'G','G':'C','T':'A'}
    # Reverse sequence
    reverse_seq = seq[::-1]

    return "".join([reverse_comp[i] for i in reverse_seq])


def get_variant(ref, alt):
    """ 
    Compute the type of genetic variant (snp, ins, del, etc.) using reference and alternative allele strings
    This function is used to aid in determining HGVS notation for each variant, therefore variant classes are restricted to {snp, ins, indel, del} to ensure they are recognized by Ensembl
    

    :param ref: reference allele
    :param alt: alternative allele
    :return: predicted variant classification
    """

    # Note: "Tandem duplications may also be represented as insertions" 
    # Source: "https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/NormalizingSmallTandemDuplications_fDG_dtSW.htm"
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1:
        # If reference length & alternative length == 1, bases must be changed -> SNP
        if alt_len == 1:
            return 'snp'
        # If reference length < alternative length, bases must be inserted 
        return 'ins' if ref[0] == alt[0] else 'indel'
    else:
        # If reference length > 1 & alternative is reverse complement: inversion
        if alt == get_reverse_complement(ref):
            return 'inv'
        # If reference length < alternative length & reference length > 1, bases must be added (opt. additional deletions)
        if ref_len < alt_len:
            # If reference in alternative allele, no deletions occured to reference, strictly insertion (else indel / insertion+deletion)
            return 'ins' if ref in alt else 'indel'
        
        # If reference length > altertnative length, bases must be removed: deletion
        elif ref_len > alt_len:
            return 'del'
        
        # If reference length == alternative length, bases are inserted and deleted but length is preserved
        else:
            return 'indel'    
    
        


async def make_ensembl_call(row, session: ClientSession):
    """ 
    Compute HGVS notation and send request to Ensembl Variant Effect Predictor API

    :param row: variant row from VCF dataframe
    :param session: ClientSession
    :return: response from Ensembl
    """

    # Specifying assembly to GRCh37 because Reference comes from v37 of fasta file, as specified in test_csv_data.txt: "/data/ref_genome/human_g1k_v37.fasta"
    vep_base = 'https://grch37.rest.ensembl.org/vep/human/hgvs/'
    var_type = row['variant_prefix']

    # Default variant type is SNP
    hgvs_notation = "{0}:g.{1}{2}{3}>{4}".format(row['CHROM'], str(row['POS']), row['REF'], var_type, row['alt_split'])
    
    # Formatting rules for HGVS nomenclature: https://varnomen.hgvs.org/recommendations/DNA/
    if var_type == 'ins':
        hgvs_notation = "{0}:g.{1}_{2}{3}ins{4}?variant_class=True".format(row['CHROM'], row['POS'], row['POS']+1, row['REF'], row['alt_split'])
    elif var_type == 'snp':
        hgvs_notation = "{0}:g.{1}{2}>{3}?variant_class=True".format(row['CHROM'], row['POS'], row['REF'], row['alt_split'])
    elif var_type == 'del':
        hgvs_notation = "{0}:g.{1}_{2}del{3}?variant_class=True".format(row['CHROM'], row['POS'], row['POS']+len(row['REF'])-1, row['alt_split'])
    elif var_type == 'indel':
        hgvs_notation = "{0}:g.{1}_{2}delins{3}?variant_class=True".format(row['CHROM'], row['POS'], row['POS']+len(row['REF'])-1, row['alt_split'])
    elif var_type == 'inv':
        hgvs_notation = "{0}:g.{1}_{2}inv?variant_class=True".format(row['CHROM'], row['POS'], row['POS']+len(row['REF'])-1)
        
    try:
        # Send request to Ensembl, specify response format as JSON
        async with session.get(url=vep_base + hgvs_notation, headers = {'content-type': 'application/json'}) as response:
            resp = await response.read()
            return resp
    except Exception as e:
        logging.error("Failed fetching data from Ensembl:" + str(e))
    

def parse_vep_dict(res):
    """ 
    Parse Ensembl VEP JSON response

    :param res: single Ensembl VEP response from which to extract relevant information
    :return: parsed information as dictionary
    """


    # Fields to collect from Ensembl responses
    minor_allele_freq = '.'
    variant_id = '.'
    gene_set = set()
    biotype_set = set()
    consequence_set = set()
    impact_set = set()
    most_severe_consequence = '.'
    allele_string = '.'
    type_variant = '.' #stored as variant_class

    try: 
        allele_string = res['allele_string']
        start = res['start']
        type_variant = res['variant_class']
        most_severe_consequence = res['most_severe_consequence']
        
        # Minor allele frequency and variant ID is stored under colocated_variants
        if 'colocated_variants' in res:
            for var in res['colocated_variants']:
                # If colocated variant matches location of variant
                if start == var['start'] and 'minor_allele_freq' in var:
                    minor_allele_freq = var['minor_allele_freq']
                    variant_id = var['id']
                    break
                
        # Gene symbols, biotypes, effect, and consequence terms fall under 'transcript_consequences' field in dictionary
        if 'transcript_consequences' in res:
            # Iterate over all transcript consequences
            for tx_cons in res['transcript_consequences']:
                gene_set.add(tx_cons['gene_symbol'])
                biotype_set.add(tx_cons['biotype'])
                impact_set.add(tx_cons['impact'])
                # All consequences are stored as list
                for term in tx_cons['consequence_terms']:
                    consequence_set.add(term)


    # Throw error if accessing key not found in dictionary
    except KeyError as err:
        print("KeyError while accessing Ensembl JSON response dict:", err)

    # Convert sets to comma-separated string
    gene_str = ",".join(gene_set) if len(gene_set) > 0 else "."
    biotype_str = ",".join(biotype_set) if len(biotype_set) > 0 else "."
    consequence_str = ",".join(consequence_set) if len(consequence_set) > 0 else "."
    impact_str = ",".join(impact_set) if len(impact_set) > 0 else "."

    return {
        'MAF':minor_allele_freq,
        'VariantID':variant_id,
        'GeneSymbols':gene_str,
        'Biotypes':biotype_str,
        'CQ':consequence_str,
        'Impact':impact_str,
        'SevereCQ':most_severe_consequence,
        'AlleleString':allele_string,
        'VariantClass':type_variant
    }
    
