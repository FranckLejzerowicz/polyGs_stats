# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd
from os.path import dirname, abspath
from polyGs_stats.polyGs_stats import get_input_files, chunks, make_fil_table, get_counts, parse_fastq, parse_gz

from pandas.testing import assert_frame_equal


class TestPolyGsStats(unittest.TestCase):

    def setUp(self):
        self.root = dirname(abspath(__file__))
        self.fastq = '%s/test.R1.fastq' % self.root
        self.gz = '%s/test.R1.fastq.gz' % self.root
        self.poly_g = {'G_k7': 1, 'G_AAAAAA': 1}
        self.output_pd = pd.DataFrame({
            'path': ['filename', 'filename', 'filename'],
            'count_variable': ['reads', 'G_k0', 'G_k1'],
            'count_value': [10, 8, 2]
        })
        self.fastq_pd = pd.DataFrame({
            'path': [self.fastq] * 6,
            'count_variable': ['G_AAAA', 'G_CCCC', 'G_k0', 'G_k4', 'G_k6', 'reads'],
            'count_value': [1, 1, 1, 1, 1, 3]
        })
        self.gz_pd = pd.DataFrame({
            'path': [self.gz] * 6,
            'count_variable': ['G_AAAA', 'G_CCCC', 'G_k0', 'G_k4', 'G_k6', 'reads'],
            'count_value': [1, 1, 1, 1, 1, 3]
        })

    def test_get_input_files(self):
        fastq_paths = sorted(get_input_files([self.root]))
        self.assertEqual(fastq_paths, [['%s/test.R1.fastq' % self.root, 'fastq'],
                                       ['%s/test.R1.fastq.gz' % self.root, 'gz']])

    def test_chunks(self):
        chunks_list = chunks([1,2,3,4,5,6], 8)
        self.assertEqual(chunks_list, [[1],[2],[3],[4],[5],[6]])
        chunks_list = chunks([1,2,3,4,5,6,7], 4)
        self.assertEqual(chunks_list, [[1,2,3,4],[5],[6],[7]])
        chunks_list = chunks([1,2,3,4,5,6,7,8], 4)
        self.assertEqual(chunks_list, [[1,2],[3,4],[5,6],[7,8]])

    def test_make_fil_table(self):
        make_fil_table_pd = make_fil_table('filename', 10, {'G_k0': 8, 'G_k1': 2})
        assert_frame_equal(make_fil_table_pd, self.output_pd)

    def test_get_counts(self):
        counts = get_counts('ACTGTGCAAAAAAGGGGGGG', ['G'], {}, 6)
        self.assertEqual(counts, {'G_k7': 1, 'G_AAAAAA': 1})
        counts = get_counts('ACTGTGCAAAAAAGGGGGGG', ['G'], self.poly_g, 6)
        self.assertEqual(counts, {'G_k7': 2, 'G_AAAAAA': 2})
        counts = get_counts('ACTGTGCAAAAAAGGGGGGG', ['G'], self.poly_g, 3)
        self.assertEqual(counts, {'G_k7': 3, 'G_AAA': 1, 'G_AAAAAA': 2})

    def test_parse_fastq(self):
        parse_fastq_pd = parse_gz(self.gz, 'G', 4).sort_values(['path', 'count_variable'])
        parse_fastq_pd.index = range(parse_fastq_pd.shape[0])
        assert_frame_equal(parse_fastq_pd, self.gz_pd)
        parse_fastq_pd = parse_fastq(self.fastq, 'G', 4).sort_values(['path', 'count_variable'])
        parse_fastq_pd.index = range(parse_fastq_pd.shape[0])
        assert_frame_equal(parse_fastq_pd, self.fastq_pd)


if __name__ == '__main__':
    unittest.main()
