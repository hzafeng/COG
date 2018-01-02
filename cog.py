import os
import re
import math
from scipy import stats
import argparse
from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
what_i_do = "a script to annotate A list of Genomes with COG database"
parser = argparse.ArgumentParser(description=what_i_do)
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='input genome file in fasta(required) e.g. <dir>*.fasta', default=None)
parser.add_argument('-db', dest='local_db', type=str, nargs=1,
                required=True, help='local COG database dir(required)', default=None)                
parser.add_argument('-o', dest='output_files', type=str, nargs=1,
                required=True, help='the output dir(required)', default=None)
parser.add_argument('-t', dest='threads_num', type=str, nargs=1,
                required=False, help='the num of threads', default='1')
def check_ids(arguments):
    lists = []
    for file in arguments["input_files"]:
        f=open(file,'r')
        data = f.read().split('\n')
        f.close()
        for i in data:
            if '>' in i:
                if i.split(' ')[0] not in lists:
                    lists.append(i.split(' ')[0])
                else:
                    print "The is a reduplicated id in your Genome",i.split(' ')[0]

args = vars(parser.parse_args())           

'''
args={  'local_db': ['/home/zjs/hhf/test/db'], 
        'output_files': ['./'], 
        'threads_num': 1, 
        'input_files': ['maxbinout_W0P5.011.fasta','maxbinout_W0P5.014.fasta']}
'''

check_ids(args)
def diamond(arguments):
    mkdb = 'makeblastdb -in '+os.path.join(args['local_db'][0],'clean.fasta')+' -dbtype prot '+'-out '+os.path.join(args['output_files'][0],'COG')
    print mkdb
#    os.system(mkdb)
    for file in arguments["input_files"]:
        blast = 'blastp'+' -db '+os.path.join(args['local_db'][0],'COG')+' -query '+file+' -evalue 10e-6 -outfmt 6 -num_alignments 1 -num_threads '+str(args['threads_num'][0])+' -out '+os.path.join(args['output_files'][0],file+'.table')
        print blast
        os.system(blast)
'''
cog = {'COG3010': 'G',
        'COG3011': 'S',
        'COG3012': 'S',
        'COG3013': 'S',......}

'''


def EnrichmentAnalysis(arguments):
    fa = open(os.path.join(args['local_db'][0],'clean.fasta'),'r')
    clean = fa.read().replace('\r','').split('\n')
    fa.close()
    fw = open(os.path.join(args['local_db'][0],'whog.txt'),'r')
    whog = fw.read().split('\n')
    fw.close()
    cog={}
    f_cog = open(os.path.join(args['local_db'][0],'cog_lists'),'r')
    coglists = f_cog.read().split('\n')
    f_cog.close()
    dict_cog={}
    for i in range(65,91):
        j=chr(i)
        each_count = 0
        for k in coglists[:-1]:
            if k.split(' ')[0] == j:
                num = int(k.split(',')[-1])
                each_count +=num
        dict_cog[j] = each_count 

    for i in whog:
        if '[' in i:
            i_ = i.split(' ')
            cog[i_[1]] = i_[0].split('[')[1][:-1]
    clean_dicts={}
    for j in clean:
        if '>' in j:
            clean_dicts[j[1:].split(' ')[0]] = j.split(' ')[1]
    for genome_id in arguments["input_files"]:
        f = open(os.path.join(args['output_files'][0],genome_id+'.table'),'r')
        each_table = f.read().split('\n')[:-1]
        f.close()
        first_line = ''
        each_cog = ''
        each_cogs = []
        for j in each_table:
            if j.split('\t')[0]+'\t'+j.split('\t')[1] == first_line:
                pass
            else:
                first_line=j.split('\t')[0]+'\t'+j.split('\t')[1]
                j_ = j.split('\t')
                cog_id = clean_dicts[j_[1]]
                hit =cog_id.split('#')
                if len(hit) == 1:
                    each_cog += cog[cog_id]
                else:
                    for k in hit:
                        each_cog += cog[k]
        for i in each_cog:
            each_cogs.append(i)
        of1 = open(os.path.join(args['output_files'][0],genome_id+'.table.count'),'w')
        for i in list(set(each_cogs)):
            of1.write(i+','+str(each_cog.count(i))+'\n')
        of1.close()

    for file in arguments["input_files"]:
        genome_id = file.split('/')[-1]
        f = open(os.path.join(args['output_files'][0],genome_id+'.table.count'),'r')
        eachdata = f.read().split('\n')
        f.close()
        each_cog_gene = 0
        for each in eachdata:
            if each != '':
                each_ = each.split(',')
                each_cog_gene += int(each_[1])
        oof = open(os.path.join(args['output_files'][0],genome_id+'.stats'),'w')
        for each in eachdata:
            if each != '':
                each_ = each.split(',')
                ff = stats.hypergeom.pmf(int(each_[1]),83675+int(dict_cog[each_[0]]),int(dict_cog[each_[0]]),each_cog_gene)
                oof.write(each_[0]+','+str(ff)+'\n')
        oof.close()
EnrichmentAnalysis(args)
