# README
SCpre-seq, an end-to-end sequence-based tool using a residue level mutation-based 3D structure information to assess single point mutation effects on protein thermodynamic stability and applying to dingle-domain monomeric proteins. More information can be checked at the [tutorial](https://mu3dsp.readthedocs.io/en/latest/) and [webserver](http://101.42.149.84:8080/).
## Clone the package

```
git clone https://github.com/hurraygong/MU3DSP.git
cd MU3DSP
```

## Install dependencies
#### 1. Install HHsuite

Using conda install hhsuite. It is also can be installed by other methods shown on official website,  If you have installed, pass this step.
```
conda install -c conda-forge -c bioconda hhsuite
```
Install HHsuite database.
```
mkdir HHsuitDB
cd HHsuitDB
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
mkdir -p UniRef30_2020_06
tar xfz UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06
rm UniRef30_2020_06_hhsuite.tar.gz
```


#### 2. Install DSSP

Using conda install DSSP programe (https://anaconda.org/salilab/dssp). It is also can be installed by other methods shown on official website,  If you have installed, pass this step.
```
conda install -c salilab dssp
```
Then to find the execute programe `mkdssp`. You can change the `mkdssp`
to `dssp` using `cp` command by yourself.
```
whereis mkdssp
```
#### 3.Install BLAST++

Install [BLAST++]( (https://anaconda.org/bioconda/blast)) using `conda`. To install this package with conda run one of the following:
```
conda install -c bioconda blast
conda install -c bioconda/label/cf201901 blast
```
It is also can be installed by other methods shown on official website,  If you have installed, pass this step.

Install BLAST++ database.
```
mkdir PsiblastDB/
cd PsiblastDB
wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz.md5
tar -zxvf swissprot.tar.gz
rm swissprot.tar.gz
rm swissprot.tar.gz.md5
```


## Usage

If provide MSA files with hhm format from HHblits. It can run with:
```python

  python MU3DSP.py -p residueposition -w wildtyperesidue -m mutationresidue -o -s sequencepath --pdbpath pdbfilepath --dssppath dsspfilepath --dsspbin mkdssp-path --psiblastbin  psiblast-path --hhblitsbin hhblits-path --psiblastout psioutfile-path --psiblastpssm pssmoutfile-path --psiblastdb swissprot-path --hhblitshhm hhmoutfile-path --hhmfile HHblits-hhm-file --printout

```
For example:
```python

    python MU3DSP.py -p 9 -w Q -m H -s ./examples/SEQ/2ocj_A.fasta --pdbpath ./examples/PDBtest --dssppath ./examples/DSSPtest --dsspbin dssp --psiblastbin  psiblast --hhblitsbin hhmake --psiblastout ./examples/psiout --psiblastpssm ./examples/pssmout --psiblastdb ./PsiblastDB/swissprot --hhblitshhm ./examples/hhmout --outfilepath ./examples/2ocj_Q9H.npy --hhmfile /root/MU3DSP/examples/2ocj_A.hhm --printout
```

If provide fasta files and HHblits is accessable.

```python

    python MU3DSP.py -p 9 -w Q -m Y --hhsuite True --outfilepath /root/MU3DSP/tmp/2ocj_A_Q9Y.npy -s /root/MU3DSP/tmp/2ocj_A.fasta --G2s --hhmfile /root/MU3DSP/examples/2ocj_A.hhm --pdbpath ./examples/PDBtest --dssppath ./examples/DSSPtest --dsspbin mkdssp --psiblastbin psiblast --psiblastout ./examples/psiout --psiblastpssm ./examples/pssmout --psiblastdb ../blastDB/swissprot --hhblitsdb ../HHblitDB/UniRef30_2020_06 --hhblitsout ./examples/hhblitout --hhblitshhm ./examples/hhmout -o --printout

```

Because the mutation-based structures are real-time, the DDG will changes with the structures changes and G2s update, we trained MU3DSP-S1676 with S1676 dataset at 12 Sep,2020 and MU3DSP-S5296 on 2 Nov,2022.
## References
