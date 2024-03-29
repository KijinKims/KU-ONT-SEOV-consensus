# KU-ONT-SEOV-consensus
program to generate consensus sequence from Nanopore sequencing dataset of SEOV

## installation
```bash
git clone https://github.com/KijinKims/KU-ONT-SEOV-consensus.git
cd KU-ONT-SEOV-consensus
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
L = "/home/kijnkims/KU-ONT-SEOV-consensus/References/SEOV_L.fasta"
```

Do the same step for M and S segments.
```bash
M = "/home/kijnkims/KU-ONT-SEOV-consensus/References/SEOV_M.fasta"
S = "/home/kijnkims/KU-ONT-SEOV-consensus/References/SEOV_S.fasta"
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
nextflow ../main.nf --fastq test.fastq --prefix test --outdir test_output
cd test_output
```

If the installation is done successfully, you could see the consensus of each segment in the `test_output` folder.

## Funding
This work was supported by the National Research Foundation of Korea (NRF) grant funded by the Korean government (MSIT) (2023R1A2C2006105), the Institute for Basic Science (IBS), Republic of Korea, under project code IBS-R801-D9-A03, and the Korea University. In addition, this study was funded by the Korea Institute of Marine Science and Technology Promotion (KIMST) fund by the Ministry of Oceans and Fisheries, Korea (20210466), and the Basic Research Program through the National Research Foundation of Korea (NRF) by the Ministry of Education (NRF2021R1I1A2049607).
