#!/usr/bin/env python3
### Import the library to run the script ####
import sys                                  #
import argparse                             #
import os                                   #
import pandas as pd
import sys
#sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))

import rotifer.core.cli as corecli
#############################################

#Implementar argparse
def parser_cli():
    parser = argparse.ArgumentParser(description=
                                     'This is an wrapper script that will run the mmseqs program multiple times with different \
                                      parameters and will create a table with all the results together.                        \
                                      Written by ggnicastro and acpguedes.')
    parser.add_argument("-c", "--coverage", default=[0.8], nargs='+', help= '''
                    Insert list of coverage to use on cluster separated by
                    spaces. [default = 0.8]
                   ''')
    parser.add_argument("-d", "--identity", default=[0.3, 0.4, 0.5, 0.6],
                    nargs='+', help= '''
                    Insert list of identity to use on cluster separated by
                    spaces [default = [0.3, 0.4, 0.5, 0.6]]
                   ''')
    parser.add_argument("-e", "--evalue", default= [0.001], nargs='+', help=
                    '''
                    Insert list of e-values to use on cluster separated by
                    spaces [default = 0.001]
                   ''')
    parser.add_argument("-b", "--basename", default="cluster" , help= '''
                    Insert basename for the files [default = cluster]
                   ''')
    parser.add_argument("file", help= '''
                    Insert sequence file
                    ''', nargs = '*'
                        )
    parser.add_argument("-m", "--clustermode", default=1, help= '''
                    Choose your cluster algorithm, see MMseqs [default = 1]
                    ''')
    parser.add_argument("-a", "--coveragemode", default=0, help= '''
                    Choose your coverage mode, see MMseqs [default = 0]
                   ''')
    parser.add_argument("--cpu", default=20, help= '''
                    Choose number of threads [default = 20]
                   ''')
    args = parser.parse_args()
    return args 

def makedir(base = 'cluster'):
    '''
    Create a base output directory
    '''
    default = base
    lsdir =  os.listdir()
    x = 1
    while True:
        if default in os.listdir():
            default = base+'_'+str(x)
            x +=1
        try:
            os.mkdir(default)
            break
        except:
            pass
    return os.path.join(os.getcwd(), default)

def runMMSeq(coverage = [0.8], identities = [0.3, 0.4, 0.5, 0.6],
             evalues = [0.001], base = 'cluster', ipt_file = '~/Projects/Response_reg/work/20180409/rsearch/rsearch.neighbor.fa',
             clustermode = 1, coveragemode = 0, cpu = 20,
             chdir = os.getcwd()):
    print('a')
    mmseq_folder = makedir('rclust')
    os.chdir(mmseq_folder)
    mmseq_clusters = os.mkdir(os.path.join(mmseq_folder, 'clusters'))
    os.system('mmseqs createdb {0} {1}.database'.format(ipt_file, base))
    os.system('mmseqs prefilter --threads {0} -s 7.5 {1}.database {1}.database {1}.prefilter'.format(
                                                                                        cpu, base))
    # enter cluster
    os.chdir(os.path.join(os.getcwd(), 'clusters'))
    print(os.getcwd())
    df_final = pd.DataFrame()
    for cov in coverage:
        for identity in identities:
            for evalue in evalues:
                number = '{0}c.{1}i.{2}e'.format(str(cov), str(identity), str(evalue))
                print ('A\n'*4 )
                os.system('''mmseqs align --threads {0} --min-seq-id {1} --cov-mode {2} \
                            -c {3} -a -e {4} ../{5}.database ../{5}.database {5}.prefilter MMseqsAln_{6}
                            '''.format(str(cpu), str(identity), str(coveragemode),
                                       str(cov), str(evalue),
                                       str(base), str(number) ))
                print ('A\n'*4 )
                os.system('mmseqs clust ../{0}.database MMSeqAln_{1} MMseqsClust_{1} --cluster-mode {2}'.format(
                                                    str(base), str(number), str(clustermode)))
                print ('A\n'*4 )
                os.system('mmseqs createtsv ../{0}.database ../{0}.database MMseqsClust_{1} MMseqsClust_{1}.tsv'.format(base, str(number))  )
                print ('A\n'*4 )
                os.system('mmseqs createtsv ../{0}.database ../{0}.database MMseqsAln_{1} MMseqs_{1}'.format(str(base),
                                                                                                             str(number)))
#                df = pd.read_csv('MMseqsClust_{0}'.format(str(number)), names = ['Cluster', 'Protein'], sep = '\t' )
#                if list(df_final.columns) == []:
#                    df_final['Protein'] = df['Protein'].values
#                set_cluster = df.groupy('Cluter').count().sort_values(by = ['Protein'],
#                                                                      ascending = False).index.tolist
#                dc = {}
#                for x in range(len(set_cluster)):
#                    dc[set_cluster[x]] = x
#                df_final[name] = df['Protein'].map(dc)
#    df_final = df_final.sort_values(list(df_final.columns),
#                                    ascending = [True]*df.shape[1])
#    df_final.to_csv(base+'.rawtable', index = False)
if __name__ == '__main__':
    runMMSeq()
