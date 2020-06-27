# Author: Julie BOGOIN

import os
import pandas

print('\ncnvkit CNV results program openning.\n')

path = '.'
files = os.listdir(path)
li = []

for name in files:
    
    if ".dedup.call.cns" in name:
        df = pandas.read_csv(name, sep='\t',index_col=None, header=[0])
        df.dropna(how='all')
        sample_name = name.split('.')
        df['sample'] = sample_name[0]
        li.append(df)

concat = pandas.concat(li, axis=0, ignore_index=True)

df_sex = pandas.read_csv('../samples.txt', header = [0], sep="\t", index_col=None)

frame = pandas.merge(concat, df_sex, left_on='sample', right_on='sample')

total = frame.shape[0] - 12

frame.rename(columns={'log2': 'log2copy_ratio'}, inplace=True)
frame.rename(columns={'chromosome': 'contig'}, inplace=True)
frame.rename(columns={'cn': 'CN'}, inplace=True)
frame.rename(columns={'probes': 'targets_number'}, inplace=True)

del frame['gene']
del frame['depth']
del frame['p_ttest']
del frame['weight']

frame.query('log2copy_ratio>1 or log2copy_ratio<-1', inplace=True)

frame['cnv_ratio'] = frame['log2copy_ratio']**2

frame['effect']='i'
frame.loc[frame.log2copy_ratio>1, 'effect'] = "duplication"
frame.loc[frame.log2copy_ratio<-1, 'effect'] = "deletion"

cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect', 'targets_number']
frame = frame[cols]

print('{0} CNV lines filtred among {1} lines found by cnvkit.'.format(frame.shape[0], total))

if os.path.isfile('cnvkit_results.csv'):
    os.remove('cnvkit_results.csv')
    print('Previous results file removed.')

frame.to_csv('cnvkit_results.csv', index=False)                                                                                           
print("cnvkit_results.csv generated.\n")
print("cnvkit CNV results job done!\n")
