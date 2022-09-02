## Installation

### Clone the package

```
git clone https://github.com/hurraygong/MU3DSP.git
cd MU3DSP
```

### Install dependencies

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

Then you can proceed to the next step.
