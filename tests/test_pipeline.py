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
from minion_amplicon_pipeline import ParseBlast

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = ArgumentParser()
    v.path = os.path.join(testpath, 'testdata')
    v.targetpath = os.path.join(v.path, 'targets')
    v.fastqfile = os.path.join(v.path, 'testdata.fastq.gz')
    v.length = 50
    v.kmersize = 27
    v.hdist = 2
    v.readlength = 400
    v.startingtime = time()
    v.scriptpath = scriptpath
    return v


@pytest.fixture()
def method_init(variables):
    method = ParseBlast(variables)
    return method


method = method_init(variables())


def variable_update():
    global method
    method = method_init(variables())


def test_bait(variables):
    method.bait()
    assert os.path.isfile(os.path.join(variables.path, 'bait', 'filteredfastq.fastq'))


def test_fastq_to_fasta(variables):
    method.fastq_to_fasta()
    assert os.path.isfile(os.path.join(variables.path, 'bait', 'filteredfastq.fasta'))


def test_make_blastdb(variables):
    method.make_blastdb()
    assert os.path.isfile(os.path.join(variables.targetpath, 'combined', 'genesofinterest.nhr'))


def test_blast(variables):
    method.blast()
    assert os.path.isfile(os.path.join(variables.path, 'blast', 'blast_report.csv'))


def test_blast_parser():
    method.blast_parser()
    assert method.resultdict['CTP2']


def test_populate_unique():
    method.populate_unique()
    assert method.read_dict['173eaa32-64e2-40ad-9f10-849e9e3b8a2a'] == 1


def test_find_unique():
    method.find_unique()
    assert len(method.unique_dict['e35S']) == 2


def test_create_lists(variables):
    method.create_lists()
    assert os.path.isfile(os.path.join(variables.path, 'outputs', 'CTP2.txt'))


def test_bin_fastq(variables):
    method.bin_fastq()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2.fastq.gz'))


def test_assemble(variables):
    method.assemble()
    assert os.path.isdir(os.path.join(variables.path, 'assemblies', 'CTP2'))


def test_bowtie_build(variables):
    method.bowtie_build()
    assert os.path.isfile(os.path.join(variables.targetpath, 'CTP2.1.bt2'))


def test_bowtie_run(variables):
    method.bowtie_run()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2', 'CTP2_sorted.bam'))


def test_samtools_index(variables):
    method.samtools_index()
    assert os.path.isfile(os.path.join(variables.path, 'sequences', 'CTP2', 'CTP2_sorted.bam.bai'))


def test_clean_folder(variables):
    directories = ['assemblies', 'bait', 'blast', 'outputs', 'sequences']
    for directory in directories:
        shutil.rmtree(os.path.join(variables.path, directory))


def test_clean_indices(variables):
    indices = glob(os.path.join(variables.targetpath, '*.bt2'))
    for index in indices:
        os.remove(index)


def test_clean_blast(variables):
    database_files = glob(os.path.join(variables.targetpath, 'combined', '*.n*'))
    for db_file in database_files:
        os.remove(db_file)


def test_clean_fai(variables):
    fai_files = glob(os.path.join(variables.targetpath, '*.fai'))
    for fai in fai_files:
        os.remove(fai)
