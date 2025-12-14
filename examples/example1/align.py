import subprocess

def align_sequences(query, reference, output, threads=4):
    """
    Run minimap2 to align query sequences to a reference genome.

    Parameters:
    - query: Path to the query FASTQ/FASTA file.
    - reference: Path to the reference genome FASTA file.
    - output: Path to the output SAM file.
    - threads: Number of threads to use (default is 4).
    """
    cmd = [
        "minimap2",
        "-t", str(threads),
        "-a", 
        reference,
        query
    ]

    with open(output, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)