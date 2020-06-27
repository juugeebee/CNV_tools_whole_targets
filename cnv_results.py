# Author: Julie BOGOIN

import os
import pandas

path = "."
dirs = os.listdir(path)
li = []

print("\nMerged CNV results program openning.\n")

for directory in dirs:
    
    if directory == 'gatkcnv_output':
        file_path = './' + directory + '/gatk_results.csv'
        if os.path.isfile(file_path):
            df_gatk = pandas.read_csv(file_path,index_col=None)
            df_gatk['cnv_tool'] = 'GATK'
            li.append(df_gatk)
            print('GATK results added.')
    
    if directory == 'cnvkit_output':
        file_path = './' + directory + '/cnvkit_results.csv'
        if os.path.isfile(file_path):
            df_cnvkit = pandas.read_csv(file_path,index_col=None)
            df_cnvkit['cnv_tool'] = 'cnvkit'
            li.append(df_cnvkit)
            print('cnvkit results added.')
    
    if directory == 'exomedepth_output':
        file_path = './' + directory + '/exomedepth_results.csv'
        if os.path.isfile(file_path):
            df_exomedepth = pandas.read_csv(file_path,index_col=None)
            df_exomedepth['cnv_tool'] = 'Exome_Depth'
            li.append(df_exomedepth)
            print('Exome Depth results added.')
    
    if directory == 'cn.mops_output':
        file_path = './' + directory + '/cn.mops_results.csv'
        if os.path.isfile(file_path):
            df_cnmops = pandas.read_csv(file_path,index_col=None)
            df_cnmops['cnv_tool'] = 'cn.mops'
            li.append(df_cnmops)
            print('cn.mops results added.')
    
    if directory == 'excavator2_output':
        file_path = './' + directory + '/excavator2_results.csv'
        if os.path.isfile(file_path):
            df_excavator2 = pandas.read_csv(file_path,index_col=None)
            df_excavator2['cnv_tool'] = 'EXCAVATOR2'
            li.append(df_excavator2)
            print('Excavator2 results added.')

frame = pandas.concat(li, axis=0, ignore_index=True)
cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect', 'cnv_tool','targets_number']
frame = frame[cols]

print('{} CNV lines found.'.format(frame.shape[0]))

if os.path.isfile('cnv_results.csv'):
    os.remove('cnv_results.csv')
    print('Previous results file removed.')

frame.to_csv('cnv_results.csv', index=False)                                                                                           
print("cnv_results.csv generated.\n")
print("Merged CNV results job done!\n")

    
