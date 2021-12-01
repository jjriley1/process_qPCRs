"""
===========================
pipeline_rMATS
===========================

:Author: Jack Riley

Overview
========

Pipeline to detect differential alternative splicing from RNA-sequencing data
by using the rMATS software (https://github.com/Xinglab/rmats-turbo/blob/v4.1.1). 

Usage
=====

To execute this pipeline in an interactive ShARC session or from a personal
device (not recommended) use the following command:

"python <path_to_pipeline_folder>/pipeline_rMATS.py make full -v5"

To execute this pipeline to the ShARC cluster, use the following:

"submit_pipeline <path_to_pipeline_folder>/pipeline_rMATS.py make full -v5"


Configuration
=============

In order for pipeline_rMATS to know what comparison to make, 2 important files are required
(without which the pipeline will not work): 

    - design.tsv = tsv file where each line a comparison of 2 named variables 
        - example design.tsv can be found at pipeline_rMATS/pipeline_rMATS/design.tsv
    - file_naming.tsv = tsv file which links the named variables to regex to the files of interest by
        - example file_naming.tsv can be found at pipeline_rMATS/pipeline_rMATS/file_naming.tsv


Input files
===========

Files are given as .bam format. Usually following completion of pipeline_utrons. (.fastq files can be
used with rMATS however this has not yet been configured). In order to compare against fixed event sets, 
files from utron_beds (output from pipeline_utrons) will also be needed.


Requirements
============

Packages required:
    - cgat/cgatcore
    - rmats-turbo (v4.1.1 was used during creation of pipeline)
        - a compatible conda environment has been created and
        stored in sudlab shared drive ('rmats-env')

These packages are installed in the conda environment "qapa-env".
R packages for final analysis/reports are installed in "qapa-env-R". 


Pipeline output
===============

rMATS produces tabular output for each type of alternative splicing event. These will be used by a custom
Rmarkdown to produce a human-readable summary of the rMATS results. Subsequent data manulation can then be
conducted manually.

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


############################
#####       main       #####
############################

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

RnaSeq.PARAMS = PARAMS


#############################
##### create prep files #####
#############################
# Only do this for the files included in design.tsv

@follows(mkdir("input_in_design"))
def filter_in_design_matrix():
    design_matrix = open("design.tsv")
    design = csv.reader(design_matrix, delimiter="\t")
    for row in design:
        comp1 = row[0]
        comp2 = row[1]
        print(comp1)
        print(comp2)
        file_naming = open("file_naming.tsv")
        file_names = csv.reader(file_naming, delimiter="\t")
        for file_row in file_names:
            for item in file_row:
                if comp1 in item: 
                    comp1_files = item.strip(comp1).strip()
                    os.system("cp -P input/" + comp1_files + " input_in_design/")
                if comp2 in item:
                    comp2_files = item.strip(comp2).strip()
                    os.system("cp -P input/" + comp2_files + " input_in_design/")
        file_naming.close()
    design_matrix.close()

@follows(filter_in_design_matrix)
def checkpaired():
    #add param lookup to check if paired analysis is being done
    paired = PARAMS["rmats_paired"]
    if paired == True :
        os.system("""for i in input_in_design/*; do
                        base=${i%-*}
                        base=${base%-*}
                        match=`ls ${base}-* | wc -l`
                        if [ $match -ne 2 ]; then
                            rm $i
                        fi
                    done""")
    
@follows(checkpaired, mkdir("prep"))
@transform("input_in_design/*.bam", regex("(.+)/(.+)-(.+)-R0.star.bam"), output=r"prep/\2-\3-prep.txt")
def create_prep_files(infile, outfile):
    to_cluster=False
    statement = "echo %(infile)s > %(outfile)s"
    P.run(statement)

@follows(create_prep_files, mkdir("prep/outputs.dir"), mkdir("post"))
@transform("prep/*.txt", regex("(.+)/(.+)-(.+)-prep.txt"), output=r"prep/outputs.dir/\2-\3-output_temp")
def run_prep(infile, outfile):
    job_threads=4
    job_memory="16G"
    job_condaenv="rmats-env"
    gtf_loc=PARAMS["gtf_path"]
    paired_reads = PARAMS["reads_paired"]
    read_length = PARAMS["reads_length"]
    final_outfile = "post/"
    star_index = PARAMS["star_index"]

    if paired_reads == True:
        statement = """rmats.py --b1 %(infile)s 
                            --gtf %(gtf_loc)s
                            -t paired
                            --readLength %(read_length)s
                            --nthread %(job_threads)s
                            --od %(final_outfile)s
                            --tmp %(outfile)s
                            --task prep"""
    if paired_reads == False: 
         statement = """rmats.py --b1 %(infile)s 
                            --gtf %(gtf_loc)s
                            -t single
                            --readLength %(read_length)s                            
                            --nthread %(job_threads)s
                            --od %(final_outfile)s
                            --tmp %(outfile)s
                            --task prep"""       
                        
    P.run(statement, job_memory=job_memory, job_threads=job_threads, job_condaenv=job_condaenv)

@follows(run_prep)
def create_post_files():
    #this works best if it is a paired analysis, but if it is not then this function will
    #still split the post step into 2 files, selecting every other file from the list 
    #where this is paired, doing so will split the conditions into each post_file.
        os.system("""ITER=1
                        for i in input_in_design/*; do
                            if [ $((ITER%2)) -ne 0 ]
                            then
                                if [ $ITER -eq 1 ]
                                then 
                                    echo "$i" >> post/condition1_post.txt      
                                else 
                                    echo ",$i" >> post/condition1_post.txt  
                                fi                        
                            else 
                                if [ $ITER -eq 2 ]
                                then 
                                    echo "$i" >> post/condition2_post.txt      
                                else 
                                    echo ",$i" >> post/condition2_post.txt
                                fi
                            fi 
                            ((ITER++))
                        done
                        echo $(tr -d '\n' < post/condition1_post.txt) > post/condition1_post.txt 
                        echo $(tr -d '\n' < post/condition2_post.txt) > post/condition2_post.txt""")

@follows(create_post_files, mkdir("post/tmp/"))
@transform("prep/outputs.dir/*/*.rmats", regex("prep/outputs.dir/(.+)-(.+)-output_temp/(.+).rmats"), output=r"post/tmp/\1-\2_\3.rmats")
def copy_rmats_prep(infile, outfile):
    statement="cp %(infile)s %(outfile)s"
    to_cluster=False
    P.run(statement)

@follows(copy_rmats_prep)
@merge(["post/condition1_post.txt", "post/condition2_post.txt"], "post/post_executed.txt")
def run_post(infiles, outfile):
    condition1, condition2 = infiles
    job_threads=16
    job_memory="8G"
    job_condaenv="rmats-env"
    gtf_loc=PARAMS["gtf_path"]
    paired_reads = PARAMS["reads_paired"]
    read_length = PARAMS["reads_length"]
    paired = PARAMS["rmats_paired"]
    novel = PARAMS["rmats_detect_novel_ss"]
    star_index = PARAMS["star_index"]

    #TO DO , ADD NOVEL SS INTEGRATION

    if paired_reads == True:
        if paired == True:
            statement = """rmats.py --b1 %(condition1)s
                            --b2 %(condition2)s 
                            --gtf %(gtf_loc)s
                            -t paired
                            --readLength %(read_length)s
                            --nthread %(job_threads)s
                            --od post/
                            --tmp post/tmp/
                            --paired-stats
                            --task post &&
                            echo "complete" >> post/post_executed.txt"""
        if paired == False:
            statement = """rmats.py --b1 %(condition1)s
                            --b2 %(condition2)s 
                            --gtf %(gtf_loc)s
                            -t paired
                            --readLength %(read_length)s
                            --nthread %(job_threads)s
                            --od post/
                            --tmp post/tmp/
                            --task post &&
                            echo "complete" >> post/post_executed.txt"""
    if paired_reads == False: 
        if paired == True:
            statement = """rmats.py --b1 %(condition1)s
                            --b2 %(condition2)s 
                            --gtf %(gtf_loc)s
                            -t single
                            --readLength %(read_length)s
                            --nthread %(job_threads)s
                            --od post/
                            --tmp post/tmp/
                            --paired-stats
                            --task post &&
                            echo "complete" >> post/post_executed.txt"""
        if paired == False:
            statement = """rmats.py --b1 %(condition1)s
                            --b2 %(condition2)s 
                            --gtf %(gtf_loc)s
                            -t single
                            --readLength %(read_length)s
                            --nthread %(job_threads)s
                            --od post/
                            --tmp post/tmp/
                            --task post &&
                            echo "complete" >> post/post_executed.txt"""

    P.run(statement, job_condaenv=job_condaenv, job_memory=job_memory, job_threads=job_threads)

###################
##### utility #####
###################

#@follows()
#def full():
#    pass


##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))