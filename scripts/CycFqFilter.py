import os
import sys
import math
import gzip
import argparse
from matplotlib import pylab
import pylab as plt
import seaborn as sns
#import numpy as np
import pandas as pd
import mappy as mp
import matplotlib.patches as patches
from Bio.SeqUtils import gc_fraction
#from icecream import ic

def stat_plot(df, ax=None, outpre="Cyclone"):
    # quality plot, boxplot
    l = len(df['quality'])+1
    qual_output = outpre + ".quality.png"
    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))

    sns.catplot(x='quality', data=df, hue=None, kind="box")
    ax.set_title('sequence quality')
    plt.savefig(qual_output) #save as png

    # gc plot, histogram
    gc_output = outpre + ".gc.png"
    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    sns.displot(x='gc', data=df, hue=None)
    ax.set_title('GC content',size=15)
    plt.savefig(gc_output)

    # length plot, histogram
    length_output = outpre + ".length.png"
    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    sns.displot(data=df, x="length", kde=True)
    ax.set_title('Read Length',size=15)
    plt.savefig(length_output)


def main(args):
    ids = []
    quality_for_plot = []
    gc_for_plot = []
    len_for_plot = []
    count = 0
    filtered = 0
    if args.plot_only == False:
        outfile = args.outpre + ".clean.fq.gz"
        if os.path.exists(outfile) == True:
            sys.exit(f"{outfile} has existed, please check!")
        fo = gzip.open(args.outpre + ".clean.fq.gz", 'wt')
    for read in mp.fastx_read(args.fastx, read_comment=False):
        qual = read[0].split("_")[-1]
        qual = float(qual)
        rlen = len(read[1])
        gc = gc_fraction(read[1])
        count += 1
        if count < args.plot_limit:
            ids.append(read[0])
            quality_for_plot.append(qual)
            len_for_plot.append(rlen)
            gc_for_plot.append(gc)
        if args.plot_only == False:
            if qual >= args.quality_cutoff and rlen >= args.length_cutoff:
                print(f"@{read[0]}\n{read[1]}\n+\n{read[2]}", file=fo)
            else:
                filtered += 1

    if args.plot_only == False:
        fo.close()
        print(f"filtering done!")
        print(f"   Total Reads: {count}")
        print(f"Filtered Reads: {filtered}")
    # list 2 dict
    data = {
            'id': ids,
            'quality': quality_for_plot,
            'gc': gc_for_plot,
            'length': len_for_plot,
    }

    # dict 2 Pandas DataFrame
    df = pd.DataFrame(data)
    return df


if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage = f"""
Usage: python3 {sys.argv[0]}  -q 7 -l 1000 -o output *.fq.gz # filter and plot
       python3 {sys.argv[0]} --plot_only -q 7 -l 1000 -o output *.fq.gz # only plot
        """
        sys.exit(usage)
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("-o", "--outpre", dest="outpre", metavar="<STR>", type=str,
                            help="prefix for outputs", required=True)
        parser.add_argument("--plot_only", dest="plot_only", action="store_true",
                            help="plot distribution for gc content, quality, read length")
        parser.add_argument("-lim", "--plot_limit", dest="plot_limit", metavar="<INT>", type=int,
                            default=50000, help="item limitation for ploting")
        parser.add_argument("-l", "--len", dest="length_cutoff", metavar="<INT>", type=int,
                            default=1000, help="filtering cutoff for read length")
        parser.add_argument("-q", "--qual", dest="quality_cutoff", metavar="<FLOAT>", type=float,
                            default=7.0, help="filtering cutoff for read quality")
        parser.add_argument("fastx", type=str,
                            help="input file, fasta or fastq")
        args = parser.parse_args()

        df = main(args)
        if args.plot_only:
            print(f"ploting...")
            stat_plot(df, ax=None,  outpre=args.outpre)
            print("ploting done!")
