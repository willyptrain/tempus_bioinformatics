import pandas as pd
import asyncio
from aiohttp import ClientSession
from utils.parser import make_ensembl_call, get_variant_statistics, explode_multiallele, create_df, concat_ensembl_results
from utils.vcf_writer import write_df_to_vcf, write_df_to_csv
import argparse
import tqdm
import time


async def main(fp: str, filter=False):
    """ 
    Main function to process VCF data, fetch information from Ensembl, and aggregate results.

    :param fp: Filepath to the VCF file.
    :param filter: Optional flag to filter variants based on specified criteria.
    :return: Processed VCF DataFrame.
    """

    # Create dataframe from VCF file, optionally filter
    vcf = create_df(fp,filter)

    # Expand multi-allelic variants into separate rows and compute variant statistics
    vcf_biallele = explode_multiallele(vcf)
    vcf_biallele = get_variant_statistics(vcf_biallele)

    # Extract rows from the DataFrame
    vcf_rows = [i[1] for i in vcf_biallele.iterrows()]    
    async_res = []
    results = []
    chunk_size = 10

    # Create an asynchronous HTTP session 
    async with ClientSession() as session:
        # Fetch Ensembl data asynchronously for a chunk of VCF row
        # Data must be chunked so as to ensure Ensembl's limit of 15 requests / second is not reached
        for idx in tqdm.tqdm(range(0, len(vcf_biallele), chunk_size), desc='Collecting data from Ensembl API'):
            async_res = await asyncio.gather(*[make_ensembl_call(row, session) for row in vcf_rows[idx:min(len(vcf_biallele),idx+chunk_size)]])
            results += async_res
            # Sleep for 0.5 seconds to ensure limit of calls is not passed
            await asyncio.sleep(0.5)

    # Append Ensembl VEP results to VCF dataframe  
    vcf_biallele = concat_ensembl_results(vcf_biallele, results)
    # Group and aggregate s.t. multi-allelic sites are no longer on separate rows
    vcf = vcf_biallele.groupby(['CHROM','POS', 'multi_allelic']).agg({col: lambda x: ",".join(set(x.astype(str))) for col in vcf_biallele.columns})
    # Drop the 'multi_allelic' column and reset the index
    vcf = vcf.drop('multi_allelic',axis='columns')
    vcf = vcf.reset_index(drop=True)

    return vcf
    


if __name__ == '__main__':

    # Start measuring time
    start = time.time()


    # Parse command-line arguments using argparse
    # -f, --filename: Must specify name of file to be annotated (required)
    # --filter: If true, only variants that pass the filter will be considered.
    parser = argparse.ArgumentParser(prog='main.py')
    parser.add_argument('-f', '--filename', type=str, required=True, help='Must specify name of file to be annotated.')
    parser.add_argument('--filter', required=False, default=False, action="store_true", help='If true, only variants that pass the filter will be considered.')
    args = parser.parse_args()
    
    fp = args.filename
    vcf = asyncio.run(main(fp, args.filter))

    # Save results as new VCF and CSV files
    write_df_to_vcf(fp, vcf)
    write_df_to_csv(vcf)
    
    # Calculate and display the time taken for processing
    diff = time.time() - start
    mins = int(diff / 60)
    seconds = round(diff % 60,2)

    if diff >= 60: 
        print()
        print('Annotations completed in', mins, 'minute(s),', seconds, 'seconds')
        print()
    else:
        print()
        print('Annotations completed in', seconds, 'seconds')
        print()
    

    
    


