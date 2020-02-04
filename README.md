# polyGs_stats

Count reads and polyGs and get the motifs before the Gs
     
 # Install

```
git clone https://github.com/FranckLejzerowicz/polyGs_stats.git
cd polyGs_stats
python setup.py build_ext --inplace --force install
```
*_Note that python should be python3_

## Input

One input is passed to `-i` or `--i-folder`: the path to a folder containing one or more .fastq of .fastq.gz files

## Outputs

A data frame written in one file:


path | count_variable | count_value
:---:|:---:|:---:
path/to/file.fastq.gz | reads | 10
path/to/file.fastq.gz | G_AAAA | 1
path/to/file.fastq.gz | G_CCC | 2
path/to/file.fastq.gz | G_k0 | 8
path/to/file.fastq.gz | G_k1 | 1 
... | ... | ... 
  
- "reads" is the total number of sequences in the file. 
- "G_AAAA" is the number of sequences having a "AAAA" motif right before the first G of the trailing Gs. 
- "G_CCC" is the number of sequences having a "CCC" motif right before the first G of the trailing Gs.
- "G_k0" is the number of sequences having no trailing Gs.
- "G_k1" is the number of sequences having a trailing Gs series of length 1.
 
## Usage

```
python3 ./polyGs_stats/scripts/run_polyGs_stats.py -i <input_path> -o <output_path> [OPTIONS]
```

By default, only the "Gs" will be searched, as this is developed to look at problems related to NovaSeq data.

*It's possible that you first need to* `chmod 755 ./polyGs_stats/script/run_polyGs_stats.py` *or to add the tool to your* `$PATH`

### Optional arguments

``` 
Usage: run_polyGs_stats.py [OPTIONS]

Options:
  -i, --i-folders TEXT       Folder of parse for fastq and/or fastq.gz files.
                             [required]
  -o, --o-table TEXT         Count table output.  [required]
  -n, --p-chunks INTEGER     Number of files chunks (default = 8).
  -b, --p-bases [A|C|G|T]    Bases to count (default = G)
  -m, --p-motif-len INTEGER  Number of bases before the starts of polyG.
  --version                  Show the version and exit.
  --help                     Show this message and exit.
```


### Bug Reports

contact `flejzerowicz@health.ucsd.edu`
