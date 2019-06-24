# bc4q
script for enhancing clusters from cd-hit, used for ITS database

dependencies: docopt, pandas, biopython, ete3, subprocess, matplotlib, tqdm

## Usage:
    bc4q.py -c cluster_file -e entrez2qiime_file -f fasta_file -o output_basename [--relaxed]

## Options:
    -h --help           show this help
    -c --clust FILE     cluster file from cd-hit
    -e --e2q FILE       entrez2qiime file
    -f --fasta FILE     fasta file with sequences
    -o --out NAME       base name for output files
    --relaxed           accept family level identity (optional)

