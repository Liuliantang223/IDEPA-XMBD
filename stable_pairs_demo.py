#!/usr/bin/
# -*- coding: utf-8 -*-

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
from matplotlib_venn import venn2, venn2_circles
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
plt.rcParams['axes.unicode_minus'] = False

# base directory
idea_dir = '/cluster2/huanglab/jiaao/individualized_analysis/IDEPA-XMBD'

# workdir
workdir = os.path.join(idea_dir, 'workdir')
sp_workdir = os.path.join(workdir, 'sp_workdir')

# script
r_bpca_path = os.path.join(idea_dir, 'r_package/bpca_v2.R')
r_penda_fdr_path = os.path.join(idea_dir, 'r_package/penda_deps_v2.R')
r_penda_path = os.path.join(idea_dir, 'r_package/penda_deps_v3.R')
reoa_path = os.path.join(idea_dir, 'deps_lib/reoa/bin/reoa')

# preprocess parameters
LOG_LABEL = True
IMPUT_LABEL = True
NA_LABEL = 'NA'
NA_RATE = 0.5
NORMALIZATION = 'z-score'

INDEX_OLD_CHAR = ['-', ' ']
INDEX_NEW_CHAR = '.'
COLUMNS_OLD_CHAR = ['-', ' ']
COLUMNS_NEW_CHAR = '.'

dataset_name = 'test_run_data'

# stable pairs parameter
data_path = os.path.join(idea_dir, f'data/{dataset_name}/methodsComp/data.csv')
normal_cohort_path = os.path.join(idea_dir, f'data/{dataset_name}/methodsComp/normal.txt')
tumor_cohort_path = os.path.join(idea_dir, f'data/{dataset_name}/methodsComp/tumor.txt')

HAS_SPECIFIC_PROTEIN = False
if HAS_SPECIFIC_PROTEIN:
    specific_protein_path = os.path.join(idea_dir, 'data/specific_protein/prognosis_unfavorable.txt')
else:
    specific_protein_path = None

# Read and preprocess data
data = utils.read_data(data_path)
normal_cohort = utils.read_normal_cohort(normal_cohort_path)
tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)
specific_protein = utils.read_specific_protein(specific_protein_path) if HAS_SPECIFIC_PROTEIN else None

rd = raw_data.InputData(data=data, 
                        normal_cohort=normal_cohort, 
                        tumor_cohort=tumor_cohort, 
                        specific_protein=specific_protein, 
                        HAS_SPECIFIC_PROTEIN=HAS_SPECIFIC_PROTEIN,
                        paired_data=None, 
                        paired_samples=None, 
                        HAS_PAIRED_DATA=None,  
                        INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                        INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                        COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                        COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                        NORMALIZATION=NORMALIZATION, 
                        LOG_LABEL=LOG_LABEL, 
                        IMPUT_LABEL=IMPUT_LABEL,
                        NA_LABEL=NA_LABEL, 
                        NA_RATE=NA_RATE)

preprocess_workdir = os.path.join(sp_workdir, 'preprocess')

# Data preprocess
rd.data_preprocess(preprocess_workdir=preprocess_workdir, r_bpca_path=r_bpca_path)

sp = stable_pairs.StablePairs(rd=rd, 
                              sp_workdir=sp_workdir, 
                              SP_THRES=5, 
                              RANDOM_VISUAL=True, 
                              NUM_VISUAL=5, 
                              CYCLE_RANKC=128, 
                              FDR=0.05, 
                              MAX_EXCEPTION_RANKC=0.05, 
                              fig_save_path=sp_workdir,
                              reoa_path=reoa_path)

sp.run_stablePairs()
