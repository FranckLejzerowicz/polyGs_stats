# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from polyGs_stats.polyGs_stats import polyGs_stats
from polyGs_stats import __version__

@click.command()
@click.option(
    "-i", "--i-folders", required=True, multiple=True,
    help="Folder of parse for fastq and/or faqst.gz files."
)
@click.option(
    "-o", "--o-table", required=True, type=str,
    help="Count table output."
)
@click.option(
    "-n", "--p-chunks", required=False, type=int, default=8,
    help="Number of files chunks (default = 8)."
)
@click.option(
    "-b", "--p-bases", required=False, default='G',
    type=click.Choice(['A', 'C', 'G', 'T']),
    help="Bases to count (default = G)"
)
@click.option(
    "-m", "--p-motif-len", required=False, type=int, default=5,
    help="Number of bases before the starts of polyG."
)
@click.version_option(__version__, prog_name="polyGs_stats")


def run_polyGs_stats(
        i_folders,
        o_table,
        p_chunks,
        p_bases,
        p_motif_len
):
    polyGs_stats(
        i_folders,
        o_table,
        p_chunks,
        p_bases,
        p_motif_len
    )


if __name__ == "__main__":
    run_polyGs_stats()
