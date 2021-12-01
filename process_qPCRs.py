"""
===========================
process_qPCRs
===========================
:Author: Jack Riley
:Date: 01/12/2021

Overview
========
Pipeline used in the analysis of qPCR data for RNA stability experiments

Usage
=====
To execute this pipeline in an interactive ShARC session or from a personal
device (not recommended) use the following command:

"python <path_to_pipeline_folder>/process_qPCRs.py make full -v5"

To execute this pipeline to the ShARC cluster, use the following:

"submit_pipeline <path_to_pipeline_folder>/process_qPCRs.py make full -v5"


Configuration
=============
Coming soon...


Input files
===========
csv files placed into inputs.dir folder. Files should be named as follows "$GENENAME_$ISOFORM.csv".


Requirements
============
Packages required:
    - Python3
    - cgat/cgatcore
    - R packages:
        - dplyr
        - tidyr
        - stringr
        - data.table
        - tibble
        - ggplot2

All reqs are installed in the conda environment "stem_utrons". 


Pipeline output
===============
Outlying techrep/biorep-filtered .csv file with raw CTs per biorep (sem = sem of techreps) --> "$GENENAME/$GENENAME_$ISOFORM_CTs.csv
Merged 2^-ddCT values --> "$GENENAME/$GENENAME_ddCTs.csv
RNA stability plot comparing retaining and splicing isoform (normalized to spliced isoform) --> "$GENENAME/$GENENAME_stability_plot.png"


Code
====
"""

###################
##### imports #####
###################

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
import csv
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile

############################################
##### Create a directory for each gene #####
############################################

@follows(mkdir("outputs.dir"))
@transform("inputs.dir/*.csv", regex("(.+)/(.+)_(.+).csv"), output=r"outputs.dir/\2/\2_\3.csv")
def create_dir_structure(infile, outfile):
    folder_name = os.path.dirname(outfile)
    statement = """mkdir -p %(folder_name)s &&
                cp %(infile)s %(outfile)s"""
    to_cluster = False
    P.run(statement)

#########################
##### RUN R SCRIPTS #####
#########################

@follows(create_dir_structure)
@transform("outputs.dir/*/*.csv", regex("(.+)/(.+)/(.+)_(.+).csv"), output=r"\1/\2/outlying_techreps_removed/otr_\3_\4.csv")
def remove_outlying_techreps(infile, outfile):
    folder_name = os.path.dirname(outfile)
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    
    statement = """mkdir -p %(folder_name)s &&
                    Rscript %(pipeline_directory)s/remove_outlying_techreps.R %(infile)s %(outfile)s"""
    P.run(statement)

@follows(remove_outlying_techreps)
@transform("outputs.dir/*/outlying_techreps_removed/otr_*", regex("(.+)/(.+)/(.+)/otr_(.+)_(.+).csv"), output=r"\1/\2/merged_techreps/merged_\4_\5.csv")
def merge_techreps(infile, outfile):
    folder_name = os.path.dirname(outfile)
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)

    statement = """mkdir -p %(folder_name)s &&
                    Rscript %(pipeline_directory)s/summarize_CTs.R %(infile)s %(outfile)s"""
    P.run(statement)

@follows(merge_techreps)
@transform("outputs.dir/*/merged_techreps/*.csv", regex("(.+)/(.+)/(.+)/merged_(.+)_(.+).csv"), output=r"\1/\2/outlying_bioreps_removed/obr_\4_\5.csv")
def remove_outlying_bioreps(infile, outfile):
    folder_name = os.path.dirname(outfile)
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)

    statement = """mkdir -p %(folder_name)s &&
                    Rscript %(pipeline_directory)s/remove_outlying_bioreps.R %(infile)s %(outfile)s"""
    P.run(statement)    

@follows(remove_outlying_bioreps)
@collate("outputs.dir/*/outlying_bioreps_removed/obr*", regex("(.+)/(.+)/(.+)/obr_(.+)_(.+).csv"), output=r"\1/\2/ddCT/\4_ddCTs.csv")
def calc_ddCTs(infile, outfile):
    folder_name = os.path.dirname(outfile)
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    file1, file2 = infile

    statement = """mkdir -p %(folder_name)s &&
                    Rscript %(pipeline_directory)s/calc_ddCT.R %(file1)s %(file2)s %(outfile)s"""
    to_cluster=False
    P.run(statement)

@follows(calc_ddCTs)
@transform("outputs.dir/*/ddCT/*", regex("(.+)/(.+)/(.+)/(.+)_ddCTs.csv"), output=r"\1/\2/\4_stability_plot.png")
def draw_stability_plot(infile, outfile):
    folder_name = os.path.dirname(outfile)
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)

    statement = """Rscript %(pipeline_directory)s/draw_stability_plots.R %(infile)s %(outfile)s"""
    to_cluster=False
    P.run(statement)


###################
##### utility #####
###################

@follows(create_dir_structure, remove_outlying_techreps, merge_techreps, remove_outlying_bioreps, calc_ddCTs, draw_stability_plot)
def full():
    pass

##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))