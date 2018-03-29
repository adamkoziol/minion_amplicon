#!/usr/bin/env python3
from argparse import ArgumentParser
from glob import glob
from time import time
import pytest
import shutil
import sys
import os

testpath = os.path.abspath(os.path.dirname(__file__))
scriptpath = os.path.join(testpath, '..')
sys.path.append(scriptpath)
from minion_amplicon_pipeline import MinionPipeline

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = ArgumentParser()
    v.path = os.path.join(testpath, 'testdata')
    v.targetpath = os.path.join(v.path, 'targets')
    v.targetfile = os.path.join(v.targetpath, 'genesofinterest.fasta')
    v.wheat = os.path.join(v.targetpath, 'genesofinterest_wheat.fasta')
    v.db = os.path.splitext(v.targetfile)[0]
    v.wheat_db = os.path.splitext(v.wheat)[0]
    v.fastq = os.path.join(v.path, 'testdata.fastq.gz')
    v.length = 50
    v.kmersize = 27
    v.hdist = 2
    v.readlength = 400
    v.startingtime = time()
    v.scriptpath = scriptpath
    v.overhang = 1
    v.assemble = True
    return v


@pytest.fixture()
def method_init(variables):
    method = MinionPipeline(variables)
    return method


method = method_init(variables())


def variable_update():
    global method
    method = method_init(variables())


def test_bait(variables):
    method.bait()
    assert os.path.isfile(os.path.join(variables.path, 'bait', 'filteredfastq.fasta'))


def test_make_blastdb(variables):
    method.make_blastdb(variables.targetfile, variables.db)
    assert os.path.isfile(os.path.join(variables.targetpath, 'genesofinterest.nhr'))


def test_blast(variables):
    method.blast(variables.db)
    assert os.path.isfile(os.path.join(variables.path, 'blast', 'blast_report.csv'))


def test_blast_parser():
    method.blast_parser()
    assert method.resultdict['CTP2']


def test_populate_unique():
    method.populate_unique()
    assert method.read_dict['173eaa32-64e2-40ad-9f10-849e9e3b8a2a'] == 2


def test_find_unique():
    method.find_unique()
    assert len(method.unique_dict['e35S']) == 1


def test_create_lists(variables):
    method.create_lists()
    assert os.path.isfile(os.path.join(variables.path, 'outputs', 'CTP2.txt'))


def test_bin_fastq(variables):
    method.bin_fastq()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2.fastq.gz'))


def test_assemble(variables):
    method.assemble()
    assert os.path.isdir(os.path.join(variables.path, 'assemblies', 'CTP2'))


def test_target_creation(variables):
    method.target_creation(variables.targetfile)
    assert os.path.isfile(os.path.join(variables.path, 'targets', 'CTP2.tfa'))


def test_bowtie_build(variables):
    method.bowtie_build()
    assert os.path.isfile(os.path.join(variables.targetpath, 'CTP2.1.bt2'))


def test_bowtie_run(variables):
    method.bowtie_run()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2', 'CTP2_sorted.bam'))


def test_samtools_index(variables):
    method.samtools_index()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2', 'CTP2_sorted.bam.bai'))


def test_extract_overhangs():
    method.extract_overhangs()
    assert len(method.rightrecords['CTP2']) == 4


def test_overhang_makedb(variables):
    method.make_blastdb(variables.wheat, variables.wheat_db)
    assert os.path.isfile(os.path.join(variables.targetpath, 'genesofinterest_wheat.nhr'))


def test_blast_overhangs(variables):
    method.blast_overhangs(variables.wheat_db)
    assert os.path.isfile(os.path.join(variables.path, 'overhangs', 'CTP2', 'CTP2_left_blast_results.csv'))


def test_parse_overhang_blast():
    method.parse_overhang_blast()
    assert '46c520b6-7af3-4df3-9688-d543b270bac8' in method.overhang_blast_results['CTP2']['left']['e35S']


def test_output_overhang_results(variables):
    method.output_overhang_results()
    assert os.path.isfile(os.path.join(variables.path, 'overhang_bins', 'CTP2', 'CTP2_e35S_left.fasta'))


def test_overhang_bowtie_run(variables):
    method.overhang_bowtie_run()
    assert os.path.isfile(os.path.join(variables.path, 'overhangs', 'CTP2', 'CTP2_e35S_left_sorted.bam'))


def test_clean_folder(variables):
    directories = ['assemblies', 'bait', 'blast', 'outputs', 'overhangs', 'overhang_bins', 'sequences']
    for directory in directories:
        shutil.rmtree(os.path.join(variables.path, directory))


def test_clean_indices(variables):
    indices = glob(os.path.join(variables.targetpath, '*.bt2'))
    for index in indices:
        os.remove(index)


def test_clean_targets(variables):
    target_files = glob(os.path.join(variables.targetpath, '*.tfa'))
    for target_file in target_files:
        os.remove(target_file)


def test_clean_fai(variables):
    fai_files = glob(os.path.join(variables.targetpath, '*.fai'))
    for fai in fai_files:
        os.remove(fai)


def test_clean_blast(variables):
    database_files = glob(os.path.join(variables.targetpath, '*.n*'))
    for db_file in database_files:
        os.remove(db_file)
