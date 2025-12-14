import argparse
import os
import primer_design as primer_design



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="dbprime: A tool for designing primers flanking indels between two sequences.")
    subcommand = parser.add_subparsers(dest="command", help="Sub-commands")

    pairwise = subcommand.add_parser("pair", help="Sub-commands for two samples")

    pairwise.add_argument("-q", "--query", required=True, type=os.path.abspath, help="Path to the query FASTA/FASTQ file.")
    pairwise.add_argument("-r", "--reference", required=True, help="Path to the reference FASTA file.")
    pairwise.add_argument("-o", "--output", required=True, type=os.path.abspath, help="Path to the output directory.")
    pairwise.add_argument("-l", "--min_indel_length", type=int, default=1, help="Minimum indel length to consider. Default is 1.")
    pairwise.add_argument("-f", "--flank", type=int, default=150, help="Length of flanking sequence to extract. Default is 150.")
    pairwise.add_argument("-s", "--product_size_range", type=int, nargs=2, default=(100, 300), help="Desired product size range for primer design. Default is (100, 300).")
    pairwise.add_argument("-p", "--num_primers", type=int, default=5, help="Number of primer pairs to design per indel. Default is 5.")
    pairwise.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for alignment. Default is 4.")
    pairwise.add_argument("-d", "--design_primers", action='store_true', help="Flag to indicate whether to design primers for the detected indels.")


    multiple = subcommand.add_parser("multi", help="Sub-commands for multiple samples")
    multiple.add_argument("-q", "--queries", required=True, type=os.path.abspath, help="Path to the directory containing multiple query FASTA/FASTQ files.")
    multiple.add_argument("-r", "--reference", required=True, type=os.path.abspath, help="Path to the reference FASTA file.")
    multiple.add_argument("-o", "--output", required=True, type=os.path.abspath, help="Path to the output directory.")
    multiple.add_argument("-l", "--min_indel_length", type=int, default=1, help="Minimum indel length to consider. Default is 1.")
    multiple.add_argument("-f", "--flank", type=int, default=150, help="Length of flanking sequence to extract. Default is 150.")
    multiple.add_argument("-s", "--product_size_range", type=int, nargs=2, default=(100, 300), help="Desired product size range for primer design. Default is (100, 300).")
    multiple.add_argument("-p", "--num_primers", type=int, default=5, help="Number of primer pairs to design per indel. Default is 5.")
    multiple.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for alignment. Default is 4.")   
    multiple.add_argument("-d", "--design_primers", action='store_true', help="Flag to indicate whether to design primers for the detected indels.")

    args = parser.parse_args()
    if args.command not in ["pair", "multi"]:
        parser.print_help()
        exit(1)
    if args.command == "pair":
        primer_design.run_pairwise_primer_design(args)
    else:
        primer_design.run_multi_primer_design(args)
        