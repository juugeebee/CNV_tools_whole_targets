# Author: Julie BOGOIN

import os
import pandas as pd
from collections import OrderedDict
import gzip

###########
#FUNCTIONS#
###########

VCF_HEADER = ['contig', 'start', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'end', 'FORMAT', 'CN']

def _count_comments(filename):
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                comments += 1
            else:
                break
    return comments

def dataframe(filename, large=True):
    if large:
        # Set the proper argument if the file is compressed.
        comp = 'gzip' if filename.endswith('.gz') else None
        # Count how many comment lines should be skipped.
        comments = _count_comments(filename)
        # Return a simple DataFrame without splitting the INFO column.
        return pd.read_table(filename, compression=comp, skiprows=comments,
                             names=VCF_HEADER, usecols=range(10))

    # Each column is a list stored as a value in this dict. The keys for this
    # dict are the VCF column names and the keys in the INFO column.
    result = OrderedDict()
    # Parse each line in the VCF file into a dict.
    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i
        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)



###########
#PRINCIPAL#
###########

print("\nGATK CNV results program openning.\n")

path = '.'
files = os.listdir(path)
li = []

for name in files:

    if "genotyped-segments" in name and ".vcf.gz.tbi" not in name:
        df = dataframe(name, large=True)
        df.dropna(how='all')
        sample_name = name.split('.')
        df['sample'] = sample_name[1]
        li.append(df)

concat = pd.concat(li, axis=0, ignore_index=True)

df_sex = pd.read_csv('../samples.txt', header = [0], sep="\t", index_col=None)

frame = pd.merge(concat, df_sex, left_on='sample', right_on='sample')

total = frame.shape[0] - 12

frame['end'] = frame['end'].str[4:]
info = frame['CN'].str.split(':', expand=True)
frame['CN'] = info[1]
frame['CN'] = frame['CN'].astype('str').astype('int')

del frame['FILTER']
del frame['ID']
del frame['REF']
del frame['QUAL']
del frame['FORMAT']
del frame['ALT']

frame.loc[frame.CN<1, 'effect'] = "deletion"
frame.loc[frame.CN>2, 'effect'] = "duplication"
frame.query('CN>2 or CN<1', inplace=True)

cols = ['sample', 'sex', 'contig', 'start', 'end', 'CN', 'effect']
frame = frame[cols]

print(frame)

print('{0} CNV lines filtred among {1} lines found by gatk.'.format(frame.shape[0], total))

if os.path.isfile('gatk_results.csv'):
    os.remove('gatk_results.csv')
    print('Previous results file removed.')

frame.to_csv('gatk_results.csv', index=False)
print("gatk_results.csv generated.\n")
print("GATK CNV results job done!\n")
