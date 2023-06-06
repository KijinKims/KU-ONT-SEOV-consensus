# KU-ONT-SEOV-consensus
program to generate consensus sequence from Nanopore sequencing dataset of SEOV

## installation
```bash
git clone https://github.com/KijinKims/seov_consensus.git
cd seov_consensus
conda env create-f environment.yml # you could use mamba or micromamba for faster installation.
conda activate seov-consensus
pip install medaka-cpu
```

## Set reference sequences path
```bash
cd References
readlink -f SEOV_L.fasta
```

Then, open the `nextflow.config` at home folder.
Paste the command output path between the quotes of `L = ""`.

For example,
```bash
L = "/home/kijnkims/seov_consensus/References/SEOV_L.fasta"
```

Do the same step for M and S segments.
```bash
M = "/home/kijnkims/seov_consensus/References/SEOV_M.fasta"
S = "/home/kijnkims/seov_consensus/References/SEOV_S.fasta"
```
Save and close the file.

## Set custom script path
Likewise, please set the custom script path to `nextflow.config`.

```bash
cd ../
readlink -f filter_indel_with_sr.py
```

Paste the command output path between the quotes of `indel_filter_script = ""`.

## Test
```bash
conda activate seov_consensus
cd Test
nextflow ../consensus.nf --fastq test.fastq --prefix test --outdir test_output
cd test_output
```

If the installation is done successfully, you could see the consensus of each segment in the `test_output` folder.

## Funding
This work was supported by the National Research Foundation of Korea (NRF) grant funded by the Korea government (MSIT) (2023R1A2C2006105) and the Institute for Basic Science Korea (IBS-R801-D92022-a03). In addition, this study was funded by Korea Institute of Marine Science & Technology Promotion (KIMST) funded by the Ministry of Oceans and Fisheries, Korea (20210466) and the Basic Research Program through the National Re-search Foundation of Korea (NRF) by the Ministry of Education (NRF-2021R1I1A2049607).
