# User guide

## Installation
1. Clone repo
```
git clone -b customized_ref https://github.com/singleron-RD/CeleScope_ATAC.git
```

2. Create conda environment and install conda packages. 
It is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda):
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install mamba -c conda-forge
cd CeleScope_ATAC
mamba create -n celescope_atac -y --file conda_pkgs.txt
```

3. Install celescope

Make sure you have activated the conda environment before running `pip install .`.
```
conda activate celescope_atac
pip install .
```

## Usage
1. Make atac referenceDir.

### Homo sapiens
```
mkdir hs_ensembl_99
cd hs_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

conda activate celescope_atac
celescope atac mkref \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.gtf.gz \
```

### Mus musculus
```
mkdir mmu_ensembl_99
cd mmu_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz

conda activate celescope_atac
celescope atac mkref \
 --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm38.99.gtf \
```

### Customized species
```
mkdir customized_species
cd customized_species

wget customized_species.fa.gz
wget customized_species.gtf.gz

gunzip customized_species.fa.gz
gunzip customized_species.gtf.gz

conda activate celescope_atac
celescope atac mkref \
 --fasta customized_species.fa \
 --gtf customized_species.gtf \
```

2. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
multi_atac \
        --mapfile mapfile \
        --thread 8 \
        --reference ref_path \
        --genomesize "2.7e+9"
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

`--refernece` Required. The path of the genome reference directory that includes genome.fa, genome.index, promoter.bed file.

`--genomesize` Required. genome size. Refer to www.genomesize.com. for example, 2.7e+9 for hs, 1.87e+9 for mm, 1.2e+8 for fruitfly.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Main output
- `{sample}/outs/` The output directory includes Filtered peak barcode matrix in hdf5 format, Barcoded and aligned fragment file, Fragment file index, Per-barcode fragment counts & metrics.