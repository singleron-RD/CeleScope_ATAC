import glob
import gzip
import importlib
import logging
import os
import re
import subprocess
import time
import unittest
import json
import sys
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps

import pandas as pd
import pysam

from celescope.__init__ import ROOT_PATH


def add_log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stderr)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper

def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj


def read_one_col(file):
    """
    Read file with one column. Strip each line.
    Returns col_list, line number
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    col1 = [item.strip() for item in col1]
    num = len(col1)
    return col1, num


def read_fasta(fasta_file, equal=False):
    """
    Args:
        equal: if True, seq in fasta must have equal length
    Returns:
        {seq_id: seq} dict
    """
    fa_dict = {}
    length = None
    with pysam.FastxFile(fasta_file) as infile:
        for index, record in enumerate(infile):
            seq = record.sequence
            if index == 0:
                length = len(seq)
            if equal:
                if length != len(seq):
                    raise Exception(f"{fasta_file} have different seq length")
            fa_dict[record.name] = seq
    return fa_dict, length


def hamming_correct(string1, string2):
    threshold = len(string1) / 10 + 1
    if hamming_distance(string1, string2) < threshold:
        return True
    return False


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def format_number(number: int) -> str:
    return format(number, ",")


def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))


class MultipleFileFoundError(Exception):
    pass

def glob_file(pattern_list: list):
    """
    glob file among pattern list
    Returns:
        PosixPath object
    Raises:
        FileNotFoundError: if no file found
        MultipleFileFound: if more than one file is found
    """
    if not isinstance(pattern_list, list):
        raise TypeError('pattern_list must be a list')

    match_list = []
    for pattern in pattern_list:
        files = glob.glob(pattern)
        if files:
            for f in files:
                match_list.append(f)
    
    if len(match_list) == 0:
        raise FileNotFoundError(f'No file found for {pattern_list}')
    
    if len(match_list) > 1:
        raise MultipleFileFoundError(
            f'More than one file found for pattern: {pattern_list}\n'
            f'File found: {match_list}'
        )
    
    return match_list[0]


def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'


def fasta_line(name, seq):
    return f'>{name}\n{seq}\n'


def find_assay_init(assay):
    init_module = importlib.import_module(f"celescope.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    file_path_dict = {
        'assay': f'{ROOT_PATH}/{assay}/{step}.py',
        'tools': f'{ROOT_PATH}/tools/{step}.py',
    }

    init_module = find_assay_init(assay)
    if os.path.exists(file_path_dict['assay']):
        step_module = importlib.import_module(f"celescope.{assay}.{step}")
    elif hasattr(init_module, 'IMPORT_DICT') and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict['tools']):
        step_module = importlib.import_module(f"celescope.tools.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module


def find_step_module_with_folder(assay, step):
    init_module = find_assay_init(assay)
    folder = ""
    try:
        step_module = importlib.import_module(f"celescope.{assay}.{step}")
        folder = assay
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"celescope.tools.{step}")
            folder = 'tools'
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")
            folder = module_path.split('.')[1]

    return step_module, folder


def sort_bam(input_bam, output_bam, threads=1, by='pos'):
    cmd = (
        f'samtools sort {input_bam} '
        f'-o {output_bam} '
        f'--threads {threads} '
        '2>&1 '
    )
    if by == "name":
        cmd += " -n"
    subprocess.check_call(cmd, shell=True)


def index_bam(input_bam):
    cmd = f"samtools index {input_bam} 2>&1 "
    subprocess.check_call(cmd, shell=True)


def add_tag(seg, id_name, correct_dict):
    """
    Args:
        seg: pysam bam segment
        id_name: {gene_id: gene_name}
        correct_dict: {low_seq: high_seq}

    Returns:
        seg with tag added

    """
    attr = seg.query_name.split('_')
    barcode = attr[0]
    umi = attr[1]
    seg.set_tag(tag='CB', value=barcode, value_type='Z')
    if umi in correct_dict:
        umi = correct_dict[umi]
    seg.set_tag(tag='UB', value=umi, value_type='Z')
    # assign to some gene
    if seg.has_tag('XT'):
        gene_id = seg.get_tag('XT')
        # if multi-mapping reads are included in original bam,
        # there are multiple gene_ids
        if ',' in gene_id:
            gene_name = [id_name[i] for i in gene_id.split(',')]
            gene_name = ','.join(gene_name)
        else:
            gene_name = id_name[gene_id]
        seg.set_tag(tag='GN', value=gene_name, value_type='Z')
        seg.set_tag(tag='GX', value=gene_id, value_type='Z')

    return seg


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


def get_assay_text(assay):
    """
    Deprecated
    add sinlge cell prefix
    deprecated
    """
    return 'Single-cell ' + assay


def check_arg_not_none(args, arg_name):
    """
    check if args.arg_name is not None
    Args:
        args: argparser args
        arg_name: argparser arg name
    Return:
        bool
    """
    arg_value = getattr(args, arg_name, None)
    if arg_value and arg_value.strip() != 'None':
        return True
    else:
        return False
    

def get_fastx_read_number(fastx_file):
    """
    get read number using pysam
    """
    n = 0
    with pysam.FastxFile(fastx_file) as f:
        for _ in f:
            n += 1
    return n

@add_log
def dump_dict_to_json(d, json_file):
    with open(json_file, 'w') as f:
        json.dump(d, f, indent=4)

@add_log
def barcode_list_stamp(barcode_list, cut=500):
    bc_list, bc_num = read_one_col(barcode_list)
    n, m = 0, 0
    stamp = defaultdict(list)
    for i in bc_list:
        m += 1
        if m <= cut:
            stamp[n].append(i)
        else:
            n += 1
            m = 1
            stamp[n].append(i)
    return stamp, bc_num