# README
SCpre-seq, an end-to-end sequence-based tool using a residue level mutation-based 3D structure information to assess single point mutation effects on protein thermodynamic stability and applying to dingle-domain monomeric proteins. More information can be checked at the [tutorial](https://scpre-seq.readthedocs.io/en/latest/).
## Install dependencies
#### 1. Install Anaconda3

Download [Anaconda](https://anaconda.org/) and install. If your computer has alread installed, pass this step.

#### 2. Install HHsuite

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

#### 3. Install DSSP

Using conda install DSSP programe. (https://anaconda.org/salilab/dssp)
```
conda install -c salilab dssp
```
Then to find the execute programe `mkdssp`. You can change the `mkdssp`
to `dssp` using `cp` command by yourself.
```
whereis mkdssp
```
#### 4.Install BLAST++

We also recommend to install [BLAST++]( (https://anaconda.org/bioconda/blast)) using `conda`. To install this package with conda run one of the following:
```
conda install -c bioconda blast
conda install -c bioconda/label/cf201901 blast
```

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
## Clone the package

```
git clone https://github.com/hurraygong/MU3DSP.git
cd MU3DSP
```
## Usage


```python


  python MU3DSP.py -p residueposition -w wildtyperesidue -m mutationresidue -o outfilename -s sequencepath --pdbpath pdbfilepath --dssppath dsspfilepath --dsspbin mkdssp-path --psiblastbin  psiblast-path --hhblitsbin hhblits-path --psiblastout psioutfile-path --psiblastpssm pssmoutfile-path --psiblastdb swissprot-path --hhblitsdb UniRef30_2020_06-path --hhblitsout hhblitoutpath --hhblitshhm hhmoutfilepath

```

## References
