# seov_consensus
program to generate consensus sequence from Nanopore sequencing dataset of SEOV

## installation
```bash
git clone https://github.com/KijinKims/seov_consensus.git
cd seov_consensus
conda env create-f environment.yml
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

## Test
```bash
cd Test
nextflow ../consensus.nf --fastq test.fastq --prefix test --outdir test_output
cd test_output
```

If the installation is done successfully, you could see the consensus of each segment in the `test_output` folder.