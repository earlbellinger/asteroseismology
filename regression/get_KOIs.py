from astroquery.simbad import Simbad
import re 
import numpy as np
import csv 
import os
import time
import pandas as pd

KICs = list(set([str.split(fname, '-')[0] for fname in os.listdir('data/feh')]))

queries = []
for KIC in KICs:
    print(KIC)
    queries += [Simbad.query_objectids('KIC ' + KIC)]
    time.sleep(1) 

queries_cpy = queries

stars = []
pnums = []
refs = []
with open('planets-reference.csv', 'r') as f:
    reader = csv.reader(f) 
    for row in reader: 
        if row[0][0] == '#': 
            continue 
        _, pl_hostname, _, pl_pnum, _, _, _, _, _, pl_disc_reflink = row
        stars += [pl_hostname]
        pnums += [pl_pnum]
        refs += [pl_disc_reflink]

KICs_sorted = sorted([int(i) for i in KICs])

unique_refs = [] 

for ii in range(len(KICs_sorted)):
    KIC_ii = KICs.index(str(KICs_sorted[ii]))
    KIC = KICs[KIC_ii]
    query = queries[KIC_ii]
    lines = str.split(str(query), '\n')
    for line in lines:
        name = re.sub('\*+', '', line)
        name = name.lstrip()
        if name in stars:
            star_idx = stars.index(name)
            indices = [i for i, x in enumerate(stars) if x == name]
            #print(indices)
            ref = list(set(refs[idx] for idx in indices))
            ads = ''
            for href_i in range(len(ref)):
                href = ref[href_i]
                first = str.split(href, ' target')[0]
                second = str.split(first, '/')[-1]
                third = str.split('='+second, '=')[-1]
                unique_refs += [third]
                if href_i > 0:
                    ads += ','
                ads += third
            print(KIC, pnums[star_idx], ads)

len(unique_refs)
len(set(unique_refs))
