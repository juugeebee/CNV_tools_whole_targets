# Author: Julie BOGOIN

import os
import pandas

print('\nEXCAVATOR2 CNV results program openning.\n')

path = '.'
folders = os.listdir(path)
li = []

for fold in folders:

    if "w10K_results." in fold:
        fold_name = fold.split('.')
        sample_name = fold_name[1]
        result_path = fold + '/' + 'Results/' + sample_name 

        files =  os.listdir(result_path)
        
        for name in files:
            if "FastCallResults_" in name:
                txt_file = result_path + '/' + name
                df = pandas.read_csv(txt_file, sep='\t',index_col=None, header=[0])
                df.dropna(how='all')
                df['sample'] = fold_name[1]
                li.append(df)

concat = pandas.concat(li, axis=0, ignore_index=True)

df_sex = pandas.read_csv('../samples.txt', header = [0], sep="\t", index_col=None)

frame = pandas.merge(concat, df_sex, left_on='sample', right_on='sample')

total = frame.shape[0] - 12

frame.rename(columns={'Chromosome': 'contig'}, inplace=True)
frame.rename(columns={'Start': 'start'}, inplace=True)
frame.rename(columns={'End': 'end'}, inplace=True)
frame.rename(columns={'Segment': 'log2copy_ratio'}, inplace=True)
frame.rename(columns={'CNF': 'cnv_ratio'}, inplace=True)

frame.query('log2copy_ratio>1 or log2copy_ratio<-1', inplace=True)

frame['effect']='i'
frame.loc[frame.log2copy_ratio>1, 'effect'] = "duplication"
frame.loc[frame.log2copy_ratio<-1, 'effect'] = "deletion"

del frame['Call']
del frame['ProbCall']

cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect']
frame = frame[cols]

print('{0} CNV lines filtred among {1} lines found by excavator2.'.format(frame.shape[0], total))

if os.path.isfile('excavator2_results.csv'):
    os.remove('excavator2_results.csv')
    print('Previous results file removed.')

frame.to_csv('excavator2_results.csv', index=False)                                                                                           
print("excavator2_results.csv generated.\n")
print("EXCAVATOR2 CNV results job done!\n")
