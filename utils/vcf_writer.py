
import pandas as pd


# All new INFO fields to be added to VCF meta-data
info_strings = [
    '##INFO=<ID=DepthCoverage,Number=.,Type=Float,Description="Depth of sequence coverage at the site of variation.">',
    '##INFO=<ID=NumVariantReads,Number=.,Type=Float,Description="Number of reads supporting the variant.">',
    '##INFO=<ID=PercVariant,Number=.,Type=Float,Description="Percentage of reads supporting the variant">',
    '##INFO=<ID=PercReference,Number=.,Type=Float,Description="Percentage of reads supporting the reference">',
    '##INFO=<ID=MAF,Number=.,Type=Float,Description="Minor Allele Frequency">',
    '##INFO=<ID=VariantID,Number=.,Type=Str,Description="Variant ID">',
    '##INFO=<ID=GeneSymbols,Number=.,Type=Str,Description="Genes affected by the variant">',
    '##INFO=<ID=Biotypes,Number=.,Type=Str,Description="Gene and/or transcript classification">',
    '##INFO=<ID=CQ,Number=.,Type=Str,Description="All potential consequences arising from variant">',
    '##INFO=<ID=Impact,Number=.,Type=Str,Description="Degree of impact arising from variant">',
    '##INFO=<ID=SevereCQ,Number=.,Type=Str,Description="Most severe consequences arising from variant.">',
    '##INFO=<ID=AlleleString,Number=.,Type=Str,Description="Shorthand allele notation">',
    '##INFO=<ID=VariantClass,Number=.,Type=Str,Description="Classification of variant">',
]

# IDs of each new INFO field
new_info_fields = ['DepthCoverage', 'NumVariantReads', 'PercVariant', 'PercReference','MAF', 'VariantID','GeneSymbols', 'Biotypes', 'CQ','Impact', 'SevereCQ', 'AlleleString', 'VariantClass']


def append_info_items(row):
    """ 
    Add new information fields to existing INFO column

    :param row: variant row from VCF dataframe
    :return: info string with new fields appended
    """
    info_string = row['INFO']
    for col in new_info_fields:
        info_string += row[col] + ';'

    return info_string


def write_df_to_vcf(source_fp: str, df: pd.DataFrame):
    """ 
    Write DataFrame / VCF table to VCF format

    :param source_fp: Original VCF file from which to append new data onto
    :param df: Dataframe to be written to VCF
    """

    # Read all lines from original VCF file
    vcf_lines = open(source_fp,'r').readlines()

    # Append new info fields to existing INFO column
    df['INFO'] = df.apply(lambda row: append_info_items(row),axis=1)

    data_cols = ""

    # Boolean set if new information fields have been written to header, used to ensure only written once
    new_fields_appended = False

    # Open new VCF file for writing
    with open('./results/annotations.vcf', 'w') as f:
        # Iterate over existing lines
        for line in vcf_lines:
            # Write new INFO fields to header once
            if not new_fields_appended and '##INFO' in line:
                f.write("\n".join(info_strings))
                f.write("\n")
                new_fields_appended = True
            f.write(line)
            
            # Stop writing if header is reached, will write remaining rows from VCF dataframe
            if '#CHROM' in line:
                data_cols = line.strip('#\n').split('\t')
                break
                
        # Write all data lines using rows from VCF table
        for row in df[data_cols].iterrows():
            f.write("\t".join(row[1].values))

    # Close new VCF file
    f.close()


def write_df_to_csv(vcf,fp='./results/annotations_table.csv'):
    """ 
    Store VCF Dataframe / Table as CSV

    :param vcf: VCF data as DataFrame
    :param fp: File path to save CSV file
    """
    
    vcf.to_csv(fp, index=False)