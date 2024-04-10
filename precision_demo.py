#!/usr/bin/
# -*- coding: utf-8 -*-

__author__ = 'Liu Yachen'

# 运行前修改r_package/bpca_v2.R
# 运行前修改deps_lib/methods_comp.py
# 运行前修改deps_lib/methods_lib.py

import click
import copy
import os
import datetime
import shutil
import random
import configparser

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib_venn import venn2,venn2_circles
from deps_lib import utils, raw_data, stable_pairs, methods_comp, penda_pro, methods_lib, similarity_lib, type1_error_lib, robustness_lib, kegg_lib, survival_lib

import rpy2
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
r_stats = importr('stats')

plt.rcParams['font.sans-serif'] = ['SimHei']  
sns.set(font_scale=1.5)  
plt.rcParams['axes.unicode_minus']=False 

''' Compare the precision of DEAs with the number of DEPs '''

########################## parameters input ##############################

# workdir
workdir = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/workdir'
mc_workdir = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/workdir/mc_workdir'


# script
r_bpca_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/r_package/bpca_v2.R'
r_penda_fdr_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/r_package/penda_deps_v2.R'
r_penda_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/r_package/penda_deps_v3.R'
reoa_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/deps_lib/reoa/bin/reoa'


# preprocess parameters
LOG_LABEL = True
IMPUT_LABEL = True
NA_LABEL = ' '
NA_RATE = 0.3
NORMALIZATION = 'z-score'

INDEX_OLD_CHAR = ['-', ' ']
INDEX_NEW_CHAR = '.'
COLUMNS_OLD_CHAR = ['-', ' ']
COLUMNS_NEW_CHAR = '.'

# method parameters
CYCLE_RANKC = 128
FDR = 0.05
MAX_EXCEPTION_RANKC = 0.05

# methods comparison parameter
data_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/data/test_run_data/methodsComp/data.csv'
normal_cohort_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/data/test_run_data/methodsComp/normal.txt'
tumor_cohort_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/data/test_run_data/methodsComp/tumor.txt'


paired_data_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/data/test_run_data/methodsComp/paired_data.csv'
paired_samples_path = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD/data/test_run_data/methodsComp/paired_samples.txt'

METHODS_LIST = ['RankComp', 'Penda', 'T-test', 'Wilcoxon', 'Quantile']

# read data
data = utils.read_data(data_path)
normal_cohort = utils.read_normal_cohort(normal_cohort_path)
tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)
paired_data = utils.read_paired_data(paired_data_path)
paired_samples = utils.read_paired_samples(paired_samples_path)

rd = raw_data.InputData(data=data, 
                        normal_cohort=normal_cohort, 
                        tumor_cohort=tumor_cohort, 
                        specific_protein=None, 
                        HAS_SPECIFIC_PROTEIN=False,
                        paired_data=paired_data, 
                        paired_samples=paired_samples, 
                        HAS_PAIRED_DATA=True,  
                        INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                        INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                        COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                        COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                        NORMALIZATION=NORMALIZATION, 
                        LOG_LABEL=LOG_LABEL, 
                        IMPUT_LABEL=IMPUT_LABEL,
                        NA_LABEL=NA_LABEL, 
                        NA_RATE=NA_RATE)

preprocess_workdir = mc_workdir + '/preprocess'

rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                    r_bpca_path=r_bpca_path)

mc = methods_comp.methodsComp(rd = rd, 
                                r_penda_path = r_penda_path, 
                                r_penda_fdr_path = r_penda_fdr_path,
                                reoa_path = reoa_path,
                                CYCLE_RANKC = CYCLE_RANKC, 
                                FDR = FDR, 
                                MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                PLOT_METHODS_COMPARE_RESULT = False,
                                METHODS_LIST = METHODS_LIST)

precision_mess, positive_mess = mc.run_methodsComp(mc_workdir = mc_workdir)
