# COG

## Usage 
cog.py [-h] -i INPUT_FILES [INPUT_FILES ...] -db LOCAL_DB -o
              OUTPUT_FILES [-t THREADS_NUM]

A script to annotate A list of Genomes with COG database

optional arguments:
  -h, --help            show this help message and exit
  
  -i INPUT_FILES [INPUT_FILES ...]      input genome file in fasta(required) e.g. <dir>*.fasta
  
  -db LOCAL_DB          local COG database dir(required)
  
  -o OUTPUT_FILES       the output dir(required)
  
  -t THREADS_NUM        the num of threads

## Dependency

### 1.Blast(best Diamond)
### 2.Python2
### 3.SciPy
### 4.Seaborn
