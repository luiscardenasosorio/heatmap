
import numpy as np
import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from pandas import DataFrame

def bigdf(input_files):
    dfs = []
    for fn in input_files:
        df = pd.read_csv(fn, sep='\t')
        df['pool'] = fn.split('.')[1]
        dfs.append(df)
    return pd.concat(dfs).sort_values('copies', ascending=False)


df = bigdf(glob.glob("uvm_combined/*.pooled.tsv"))


def plot_gene(df, column):

    df2 =df.pivot_table(index = "pool", columns = column, values= "copies", aggfunc= np.sum).fillna(0)
    mask = df2 == 0
    sns.clustermap(df2, standard_scale= 0, cmap="coolwarm", mask=mask).fig.suptitle("V Usage {}".format(df.iloc[0].subject))
    plt.savefig('uvm_combined/{}_{}.pdf'.format(df.iloc[0].subject, column), bbox_inches= 'tight')


def get_aa_usage(df):
    usage = df['cdr3_aa'].apply(lambda cdr3: pd.Series(
        Counter(cdr3))
    ).fillna(0)
    usage = usage.mul(df['copies'], axis = 0)
    df2 = df[['pool']].join(usage).groupby('pool').agg(np.sum)
    mask = df2 == 0
    sns.clustermap(df2, standard_scale=0, cmap="coolwarm", mask= mask).\
        fig.suptitle("CDR3_AA Usage {}".format(df.iloc[0].subject))
    plt.savefig('uvm_combined/{}_cdr3aa.pdf'.format(df.iloc[0].subject), bbox_inches= 'tight')


df.groupby('subject').apply(get_aa_usage)
df.groupby("subject").apply(lambda x: plot_gene(x,"v_gene"))
df.groupby("subject").apply(lambda x: plot_gene(x,"j_gene"))
