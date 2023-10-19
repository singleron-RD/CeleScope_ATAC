# User guide

## Installation
1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope_ATAC.git
```

2. Create conda environment and install conda packages. 
It is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda):
```
conda install mamba
cd CeleScope_ATAC
mamba create -n celescope_atac -y --file conda_pkgs.txt
```


3. Install celescope_atac

Make sure you have activated the `celescope_atac` conda environment before running `pip install celescope_atac`. 
```
conda activate celescope_atac
pip install celescope_atac
```

## Usage
1. Make a atac reference

### Homo sapiens

```
mkdir -p /genome/scATAC
cd /genome/scATAC
wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz
tar -xvzf Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz
wget http://cistrome.org/~galib/MAESTRO/references/giggle.all.tar.gz
tar -xvzf giggle.all.tar.gz
cd Refdata_scATAC_MAESTRO_GRCh38_1.1.0
chromap -i -r GRCh38_genome.fa -o GRCh38_chromap.index
```

### Mus musculus

```
mkdir -p /genome/scATAC
cd /genome/scATAC
wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz
tar -xvzf Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz
wget http://cistrome.org/~galib/MAESTRO/references/giggle.all.tar.gz
tar -xvzf giggle.all.tar.gz
cd Refdata_scATAC_MAESTRO_GRCm38_1.1.0
chromap -i -r GRCm38_genome.fa -o GRCm38_chromap.index
```

2. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
multi_atac \
        --mapfile ./20230427.mapfile \
        --thread 8 \
        --chemistry atac \
        --reference /SGRNJ06/randd/USER/cjj/soft/scpipe/ref/latest/hg38.fa \
        --organism hg38 \
        --mod shell
``` 
`--mapfile` Required.  Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1
fastq_prefix2	fastq_dir2	sample1
fastq_prefix3	fastq_dir1	sample2

$ls fastq_dir1
fastq_prefix1_1.fq.gz   fastq_prefix1_2.fq.gz	fastq_prefix1_3.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz	fastq_prefix3_3.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz	fastq_prefix2_3.fq.gz
```

`--refernece` Required. The path of the genome reference directory after downloading from ucsc.

`--organism` Required. hg38 or mm10.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Main output
- `{sample}/02.atac/{sample}_clonetypes.csv` The output directory.