#!/usr/bin/env python3

import os
import io
import dnaio

infile = 'HG03486_hgsvc_pbsq2-ccs_1000_scV12-pereg.fasta.backup'
outfile = 'HG03486_hgsvc_pbsq2-ccs_1000_scV12-pereg.fasta'

target_cluster = 'cluster1'

head_buffer = io.StringIO()
tail_buffer = io.StringIO()

current_buffer = head_buffer

cluster_buffer = io.StringIO()

with open(infile, 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            if line.strip() == '>{}'.format(target_cluster):
                current_buffer = cluster_buffer
                continue
            elif cluster_buffer.tell() > 0:
                current_buffer = tail_buffer
            else:
                current_buffer = head_buffer
        current_buffer.write(line)

splitter = 'N' * 100

cluster_seq = cluster_buffer.getvalue().replace('\n', '').split(splitter)
print('num contigs: {}'.format(len(cluster_seq)))
seq_sizes = [len(s) for s in cluster_seq]
print(seq_sizes)

print(head_buffer.tell())
print(tail_buffer.tell())

cluster_buffer = []

suffices = ['A', 'B', 'C']
suffix_idx = 0

with dnaio.FastaWriter(outfile, line_length=80) as fasta:
    current_block = 0
    for seq_size, seq in zip(seq_sizes, cluster_seq):
        current_block += seq_size
        cluster_buffer.append(seq)
        if current_block > 150e6:
            print('Writing block size: {}'.format(current_block))
            block_name = target_cluster + suffices[suffix_idx]
            fasta.write(block_name, splitter.join(cluster_buffer))
            cluster_buffer = []
            current_block = 0
            suffix_idx += 1
        

    print('Writing block size: {}'.format(current_block))
    block_name = target_cluster + suffices[suffix_idx]
    fasta.write(block_name, splitter.join(cluster_buffer))
    
with open(outfile, 'a') as fasta:
    _ = fasta.write(tail_buffer.getvalue())