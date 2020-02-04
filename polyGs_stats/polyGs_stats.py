# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import gzip
import pandas as pd
import subprocess
import multiprocessing


def chunks(input_files: list, chunk_number: int) -> list:
    # Adapted from:
    # https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    if len(input_files) <= chunk_number:
        files_chunks = [[x] for x in input_files]
    else:
        n = len(input_files) // chunk_number
        m = len(input_files) % chunk_number
        if n > 1:
            files_chunks = [input_files[i:i + n] for i in range(0, len(input_files), n)]
        else:
            files_chunks = [input_files[:m+1]] + [[x] for x in input_files[m+1:]]

    return files_chunks


def get_input_files(i_folders: list) -> list:
    """
    Collect the files paths of the fastqs.

    :param i_folders: path to folder containing the fastq files.
    :return: list of actual fastq / fastq.gz files.
    """
    input_files = []
    for i_folder in i_folders:
        for root, dirs, files in os.walk(i_folder):
            for fil in files:
                p = root + '/' + fil
                if fil.endswith('.fastq'):
                    typ = 'fastq'
                elif fil.endswith('.fastq.gz'):
                    typ = 'gz'
                else:
                    continue
                input_files.append([p, typ])
    return input_files


def write_fil_table(fils_tabs: list, out: str) -> None:
    fils_tab = pd.concat(fils_tabs)
    fils_tab.to_csv(out, index=False, sep='\t')


def run_counting(out, files, p_bases, p_motif_len):
    fils_tabs = []
    for (fil, typ) in files:
        if typ == 'fastq':
            tab = parse_fastq(fil, p_bases, p_motif_len)
        elif typ == 'gz':
            tab = parse_gz(fil, p_bases, p_motif_len)
        else:
            print('Why is the file not ending with fastq or fastq.gz?')
            sys.exit(1)
        fils_tabs.append(tab)
    write_fil_table(fils_tabs, out)
    return fils_tabs


def make_fil_table(fil: str, all_reads: int, poly_gs: dict) -> pd.DataFrame:
    tab_ls = [[fil, 'reads', all_reads]]
    for key, value in poly_gs.items():
        tab_ls.append([fil, key, value])
    tab = pd.DataFrame(tab_ls, columns = ['path', 'count_variable', 'count_value'])
    return tab


def get_counts(cur_seq: str, p_bases: list, poly_gs: dict, p_motif_len: int) -> dict:
    for p_base in p_bases:
        cur_seq_poly_g_strip = cur_seq.rstrip(p_base)
        poly_g_len = len(cur_seq) - len(cur_seq_poly_g_strip)
        poly_g_key = '%s_k%s' % (p_base, poly_g_len)
        if poly_g_key in poly_gs:
            poly_gs[poly_g_key] += 1
        else:
            poly_gs[poly_g_key] = 1
        if poly_g_len and len(cur_seq_poly_g_strip) >= p_motif_len:
            cur_seq_poly_g_strip_motif = cur_seq_poly_g_strip[-p_motif_len:]
            motif_key = '%s_%s' % (p_base, cur_seq_poly_g_strip_motif)
            if motif_key in poly_gs:
                poly_gs[motif_key] += 1
            else:
                poly_gs[motif_key] = 1
    return poly_gs


def parse_fastq(fil, p_bases, p_motif_len):
    poly_gs = {}
    all_reads = 0
    try:
        with open(fil) as f:
            idx = 0
            for ldx, line in enumerate(f):
                idx += 1
                if idx == 2:
                    all_reads += 1
                    poly_gs = get_counts(line.strip(), p_bases, poly_gs, p_motif_len)
                if idx == 4:
                    idx = 0
    except OSError:
        print('OSError with:', fil)
        sys.exit(1)

    tab = make_fil_table(fil, all_reads, poly_gs)
    return tab


def parse_gz(fil, p_bases, p_motif_len):
    poly_gs = {}
    all_reads = 0
    try:
        with gzip.open(fil, 'rb') as f:
            idx = 0
            for ldx, line in enumerate(f):
                idx += 1
                if idx == 2:
                    all_reads += 1
                    poly_gs = get_counts(str(line.decode().strip()), p_bases, poly_gs, p_motif_len)
                if idx == 4:
                    idx = 0
    except OSError:
        print('OSError with:', fil)
        sys.exit(1)

    tab = make_fil_table(fil, all_reads, poly_gs)
    return tab


def polyGs_stats(
        i_folders,
        o_table,
        p_chunks,
        p_bases,
        p_motif_len
):
    input_files = get_input_files(i_folders)
    files_chunks = chunks(input_files, p_chunks)

    jobs = []
    all_outs = []
    for n, chunk in enumerate(files_chunks):
        cur_out = '%s_%s.txt' % (o_table, p_chunks)
        all_outs.append(cur_out)
        p = multiprocessing.Process(target=run_counting, args=(cur_out, chunk, p_bases, p_motif_len,))
        p.start()
        jobs.append(p)

    for j in jobs:
        j.join()

    all_outs_pd = pd.concat([pd.read_csv(x, header=0, sep='\t') for x in all_outs])
    all_outs_pd.to_csv(o_table, index=False, sep='\t')
    for all_out in all_outs:
        subprocess.call(['rm', all_out])


