import os
import logging
import pandas as pd
from pathlib import Path


try:
    from Bio import SeqIO    
except ModuleNotFoundError:
    raise ModuleNotFoundError("Biopython is required for this script. Please install it via 'pip install biopython' or 'conda install bioconda::biopython'.")
   

try:
    import pysam
except ModuleNotFoundError:
    raise ModuleNotFoundError("pysam is required for this script. Please install it via 'pip install pysam' or 'conda install -c bioconda pysam'.")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)

def load_fasta(path, COMPRESSED=False):
    """It loads fasta into a SeqRecord and dumps into a dict."""
    seqs = {}    
    if COMPRESSED:
        import gzip
        with gzip.open(path, 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seqs[record.id] = record
        return seqs
        
    
    for record in SeqIO.parse(path, "fasta"):
        seqs[record.id] = record
    return seqs



def export_primer_results_to_file(candidates_df, output_file, outdir):
    """Export primer design results to a CSV file. Also dump json to outdir..."""
    logging.info(f"Exporting primer design results to output directory '{outdir}'")
    Path(outdir).mkdir(parents=True, exist_ok=True)
    
    #make the json whatever output file is but json extension
    tsv_output_path = output_file
    json_output_path = output_file.replace(".tsv", ".json")

    logging.info(f"Exporting primer design results to '{tsv_output_path}' and '{json_output_path}'")

    candidates_df.to_csv(tsv_output_path, sep='\t', index=False)
    candidates_df.to_json(json_output_path, orient='records', lines=True)
    print(f"Primer design results exported to '{tsv_output_path}' and '{json_output_path}'")
    
    return tsv_output_path, json_output_path    

def parse_samfile(sam_file, query_seqs, ref_seqs, min_indel_length=1, flank=150, max_indel_size=10000):
    """Parse SAM file, return dataframe of indel candidates."""

    candidates = []
    index = 0
    with pysam.AlignmentFile(sam_file, "r") as samfile:
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            
            ref_id = samfile.get_reference_name(read.reference_id)
            query_id = read.query_name
            ref_start = read.reference_start
            ref_pos = ref_start
            query_pos = 0
            cigar_tuples = read.cigartuples
            for(cigar_op_code, length) in cigar_tuples:
                #skip matches and mismatches
                if cigar_op_code == 0 or cigar_op_code == 7 or cigar_op_code == 8:
                    ref_pos += length
                    query_pos += length
                if cigar_op_code not in [1, 2]:
                    if cigar_op_code in [3, 6]:
                        ref_pos += length
                    elif cigar_op_code == 4 or cigar_op_code == 5:
                        query_pos += length
                    else:
                        ref_pos += length
                        query_pos += length               
                
                if cigar_op_code == 1:
                    #we're about to parse an insertion
                    if length >= min_indel_length:
                        left_flank_ref = ref_seqs[ref_id].seq[max(0, ref_pos - flank):ref_pos]
                        right_flank_ref = ref_seqs[ref_id].seq[ref_pos:min(len(ref_seqs[ref_id]), ref_pos + flank)]
                        ref_seq = str(left_flank_ref) + str(right_flank_ref)
                        query_seq = query_seqs[query_id].seq[max(0, query_pos - flank):min(len(query_seqs[query_id]), query_pos + length + flank)]
                        inserted_seq = query_seqs[query_id].seq[query_pos:query_pos + length]   
                        candidates.append( (index, 'INS', ref_id, ref_pos, length, str(ref_seq), str(query_seq), str(inserted_seq), query_id, str(left_flank_ref), str(right_flank_ref)) )
                        index += 1
                    query_pos += length
                else:
                    if length >= min_indel_length:
                        #we're about to parse a deletction           
                        left_flank_ref = ref_seqs[ref_id].seq[max(0, ref_pos - flank):ref_pos]
                        right_flank_ref = ref_seqs[ref_id].seq[ref_pos + length:min(len(ref_seqs[ref_id]), ref_pos + length + flank)]
                        ref_seq = str(left_flank_ref) + str(right_flank_ref)
                        query_seq = query_seqs[query_id].seq[max(0, query_pos - flank):min(len(query_seqs[query_id]), query_pos + flank)]
                        deleted_seq = ref_seqs[ref_id].seq[ref_pos:ref_pos + length]
                        candidates.append( (index, 'DEL', ref_id, ref_pos, length, str(ref_seq), str(query_seq), str(deleted_seq), query_id, str(left_flank_ref), str(right_flank_ref)) )
                        index += 1
                    ref_pos += length
                    
                    
            if len(candidates) > max_indel_size:
                break
           
    df = pd.DataFrame(candidates, columns=['Index', 'Type', 'Ref_ID', 'Position', 'Length', 'Ref_Seq', 'Query_Seq', 'Indel_Seq', 'Query_ID', 'Left_Flank', 'Right_Flank']) 
    return df
                
                    
                