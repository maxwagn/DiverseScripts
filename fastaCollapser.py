# !/usr/bin/env python

## scripts that collapses fasta sequences and returns list of identical sequences;
## modified after: https://bioinformatics.stackexchange.com/questions/2817/how-do-i-find-identical-sequences-in-a-fasta-file

import sys
from Bio import SeqIO
from collections import defaultdict

### usage: python3 fasta.Collapser.py {input.fasta} {out-prefix}

fasta_input = sys.argv[1]
output_prefix = sys.argv[2]

dedup_records = defaultdict(list)

for record in SeqIO.parse(str(fasta_input), "fasta"):
    # Use the sequence as the key and then have a list of id's as the value
    dedup_records[str(record.seq)].append(record.id)

with open("{}_allIDs.fasta".format(output_prefix), 'w') as output1,open("{}_uniqueIDs.fasta".format(output_prefix), 'w') as output2, open("{}_haplotype_ids.txt".format(output_prefix), 'w') as id_list:
    id_list.write("### list of identical sequence ids:\n")
    for seq, ids in dedup_records.items():
        output1.write(">{}\n".format('|'.join(ids)))
        output1.write(seq + "\n")
        output2.write(">{}\n".format(ids[0]))
        output2.write(seq + "\n")
        id_list.write('{}\n'.format(','.join(ids)))
