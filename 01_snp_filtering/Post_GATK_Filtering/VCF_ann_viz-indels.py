#!/usr/bin/env python3
"""A script to plot distributions of variant annotations in VCF."""

import sys
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import cyvcf2
from cyvcf2 import VCF

# User provided command line arguments
vcf_fp = sys.argv[1]
out_dir = sys.argv[2]
out_prefix = sys.argv[3]

def plot_density_hist(ann_list, plot_title):
    """Plot density/histogram of annotation"""
    num_bins = int(len(ann_list)/5)
    sns.distplot(ann_list, hist=True, kde=True,
         bins=num_bins, color = 'darkblue',
         hist_kws={'edgecolor':'black'},
         kde_kws={'shade': True, 'linewidth': 1}).set(title=plot_title)

# Load vcf and prepare annotation values
dp = []
qd = []
fs = []
sor = []
mq = []
mqranksum = []
rpranksum = []
for variant in VCF(vcf_fp):
    dp.append(variant.INFO.get('DP'))
    qd.append(variant.INFO.get('QD'))
    fs.append(variant.INFO.get('FS'))
    sor.append(variant.INFO.get('SOR'))
    mq.append(variant.INFO.get('MQ'))
    mqranksum.append(variant.INFO.get('MQRankSum'))
    rpranksum.append(variant.INFO.get('ReadPosRankSum'))

# DP
plot_density_hist(dp, 'DP')
dp_out_fp = out_dir + '/' + out_prefix + '_DP.png'
plt.savefig(dp_out_fp)

# QD
plot_density_hist(qd, 'QualByDepth (QD)')
qd_out_fp = out_dir + '/' + out_prefix + '_QD.png'
plt.savefig(qd_out_fp)

# FS
plot_density_hist(fs, 'FisherStrand (FS)')
fs_out_fp = out_dir + '/' + out_prefix + '_FS.png'
plt.savefig(fs_out_fp)

# SOR
plot_density_hist(sor, 'StrandOddsRatio (SOR)')
sor_out_fp = out_dir + '/' + out_prefix + '_SOR.png'
plt.savefig(sor_out_fp)

# MQ
plot_density_hist(mq, 'RMSMappingQuality (MQ)')
mq_out_fp = out_dir + '/' + out_prefix + '_MQ.png'
plt.savefig(mq_out_fp)

# MQRankSum
plot_density_hist(mqranksum, 'MappingQualityRankSumTest (MQRankSum)')
mqrs_out_fp = out_dir + '/' + out_prefix + '_MQRankSum.png'
plt.savefig(mqrs_out_fp)

# ReadPosRankSum
plot_density_hist(rpranksum, 'ReadPosRankSumTest (ReadPosRankSum)')
rprs_out_fp = out_dir + '/' + out_prefix + '_ReadPosRankSum.png'
plt.savefig(rprs_out_fp)
