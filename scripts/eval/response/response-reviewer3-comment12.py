#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

# plotting script for response letter
# created by Tobias Marschall

d = pd.read_csv('pete-zenodo/variants_freeze3_sv_insdel.tsv.gz', sep='\t')

d = d.assign(PAV=lambda df: True)
d = d.assign(PANGENIE_STRICT=lambda df: df.PG_CONF==4)
d = d.assign(PANGENIE_LENIENT=lambda df: df.PG_CONF>0)
d = d.assign(ILLUMINA=lambda df: df['1KGHC_OVERLAP']=='OVR')

counts = d.groupby(by=['PAV', 'PANGENIE_STRICT','PANGENIE_LENIENT','ILLUMINA']).size()
upsetplot.plot(counts, sort_by='cardinality')
plt.savefig('response-reviewer3-comment12.pdf')

print(counts)


