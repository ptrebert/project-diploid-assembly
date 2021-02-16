#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

# plotting script for response letter
# created by Tobias Marschall

d = pd.read_csv('pete-zenodo/variants_freeze3_sv_insdel.tsv.gz', sep='\t')

d = d.assign(PAV=lambda df: True)
d = d.assign(IN_AUDANO=lambda df: ~df.AUDANO2019.isna())
d = d.assign(IN_CHAISSON=lambda df: ~df.HGSVC1.isna())

print(d.groupby(by=['PAV', 'IN_AUDANO']).size())
print(d.groupby(by=['PAV', 'IN_CHAISSON']).size())

counts = d.groupby(by=['PAV', 'IN_AUDANO','IN_CHAISSON']).size()
print(counts)
upsetplot.plot(counts, sort_by='cardinality')
plt.savefig('response-reviewer3-comment2.pdf')


