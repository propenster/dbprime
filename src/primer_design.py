
import os
import argparse
import json
import yaml
import csv
import logging
from pathlib import Path
from collections import defaultdict
import align as align 
import utils as utils
from concurrent.futures import ThreadPoolExecutor, as_completed


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
 

try:
    import primer3
except ModuleNotFoundError:
    raise ModuleNotFoundError("primer3-py is required for this script. Please install it via 'pip install primer3-py' or 'conda install -c bioconda primer3-py'.")


def design_primer(left_flank, right_flank, product_size_range=(100, 300), num_primers=5):
    """Design primers using primer3 for given flanking sequences."""
    template = left_flank + right_flank
    product_start = len(left_flank) - 10 if len(left_flank) > 10 else 0
    product_end = len(left_flank) + 10 + product_size_range[1]
    product_length = product_end - product_start
    primer3_result = primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': template,
            'SEQUENCE_TARGET': [product_start, product_length],
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [product_start, product_end - product_start]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_NUM_RETURN': num_primers,
            'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range]
        }
    )
    if primer3_result['PRIMER_PAIR_NUM_RETURNED'] == 0:
        return None
    return {
        'left_primer': primer3_result['PRIMER_LEFT_0_SEQUENCE'],
        'right_primer': primer3_result['PRIMER_RIGHT_0_SEQUENCE'],
        'left_tm': primer3_result['PRIMER_LEFT_0_TM'],
        'right_tm': primer3_result['PRIMER_RIGHT_0_TM'],
        'product_size': primer3_result['PRIMER_PAIR_0_PRODUCT_SIZE']
    }



def run_pairwise_primer_design(args):
    """ Run pairwise primer design workflow. """
    outdir = args.output
    logging.info("Running pairwise primer design...")
    Path(outdir).mkdir(parents=True, exist_ok=True)
    sam_path = os.path.join(outdir, "alignment.sam")
    align.align_sequences(args.query, args.reference, sam_path, threads=args.threads)
    logging.info(f"Alignment completed. See output SAM file at {sam_path}")
    logging.info("Loading sequences")
    ref_seqs = utils.load_fasta(args.reference)
    query_seqs = utils.load_fasta(args.query)
    logging.info("Parsing SAM file for indels")
    candidates_df = utils.parse_samfile(sam_path, query_seqs, ref_seqs, min_indel_length=args.min_indel_length, flank=args.flank)
    logging.info(f"Found {len(candidates_df)} indel candidates.")
    if args.design_primers:
        logging.info("Designing primers for indel candidates")
        for c in candidates_df.itertuples():
            left_flank = c.left_flank
            right_flank = c.right_flank
            primer_info = design_primer(left_flank, right_flank, product_size_range=args.product_size_range, num_primers=args.num_primers)
            if primer_info:
                candidates_df["Primers"] = str(primer_info)
                logging.info(f"Primer design successful for candidate {c.Index}: {primer_info}")
                
                
    output_file = os.path.join(outdir, "indel_candidates_with_primers.tsv")
    utils.export_primer_results_to_file(candidates_df, output_file, outdir)
    logging.info(f"Primer design results saved to {output_file}")

def run_primer_design(fasta_path, outdir, reference, design_primers=True, min_indel_length=1, flank=150, product_size_range=(100,300), num_primers=5, threads=4):
    logging.info(f"Processing sample: {fasta_path}")
    compressed = str(fasta_path).endswith(".gz")
    sample_output = os.path.join(Path(outdir), Path(fasta_path).stem)
    sample_output.mkdir(parents=True, exist_ok=True)
    logging.info("Query vs Reference alignment")
    sam_path = os.path.join(sample_output, "alignment.sam")

    align.align_sequences(fasta_path, reference, sam_path, threads=threads)

    logging.info(f"Alignment completed. See output SAM file at {sam_path}")
    logging.info("Loading sequences")
    ref_seqs = utils.load_fasta(reference, COMPRESSED=True if reference.endswith(".gz") else False)
    query_seqs = utils.load_fasta(fasta_path, COMPRESSED=compressed)
    logging.info("Parsing SAM file for indels")
    candidates_df = utils.parse_samfile(sam_path, query_seqs, ref_seqs, min_indel_length=min_indel_length, flank=flank)
    logging.info(f"Found {len(candidates_df)} indel candidates.")   
    
    if design_primers:
        logging.info("Designing primers for indel candidates")
        for c in candidates_df.itertuples():
            left_flank = c.left_flank
            right_flank = c.right_flank
            primer_info = design_primer(left_flank, right_flank, product_size_range=product_size_range, num_primers=num_primers)
            if primer_info:
                candidates_df["Primers"] = str(primer_info)
                logging.info(f"Primer design successful for candidate {c.Index}: {primer_info}")



        output_file = os.path.join(sample_output, "indel_candidates_with_primers.tsv")
        tsv_output_path, json_output_path =utils.export_primer_results_to_file(candidates_df, output_file, sample_output)
        logging.info(f"Primer design results for sample saved to {output_file}")


    return fasta_path.stem, tsv_output_path, json_output_path, candidates_df

    
def run_multi_primer_design(args):
    """ Run multi sample primer design workflow. """
    outdir = args.output
    run_results = {}
    logging.info("Running multi sample primer design...")
    from tqdm import tqdm
    
    sample_files = [os.path.join(args.queries, f) for f in os.listdir(args.queries) if f.endswith(('.fasta', '.fa', '.fastq', '.fq', '.fa.gz', '.fasta.gz', '.fq.gz', '.fastq.gz'))]
    fasta_files = list(set(sample_files))
    if not fasta_files:
        logging.error(f"No FASTA files found in the queries directory '{args.queries}'.")
        raise RuntimeError("No FASTA files found in the queries directory.")
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(run_primer_design, fasta_path, outdir, args.reference, design_primers=args.design_primers, min_indel_length=args.min_indel_length, flank=args.flank, product_size_range=args.product_size_range, num_primers=args.num_primers, threads=args.threads) for fasta_path in sorted(fasta_files)
        }
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Processing fasta samples",
            unit="sample"
        ):
            
            sample_name, tsv_output_path, json_output_path, candidates_df = future.result()
            run_results[sample_name] = {
                "tsv": tsv_output_path,
                "json": json_output_path,
                "num_indel_candidates": len(candidates_df)
            }
    logging.info("Multi-sample primer design completed.")
    
    summary_path = os.path.join(outdir, "multi_sample_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(run_results, f, indent=4)
    logging.info(f"Multi-sample primer design summary saved to {summary_path}")
