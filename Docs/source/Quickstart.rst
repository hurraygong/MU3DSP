Quickstart
----------

Quikly Run MU3DSP
^^^^^^^^^^^^^^^^^
This is the parameters to run MU3DSP.

::

    parser = argparse.ArgumentParser(description='Download PDBfiles from G2s & get features')
    # parser.add_argument('-l', '--variant-site', dest='variant_list', type=str, help='A list of variants, one per line in the format "POS WT MUT", a file')

    parser.add_argument('-w', '--variant-wildtype', dest='variant_wildtype', type=str,
                        required=True, help='wild-type residue, for example residue "A"')
    parser.add_argument('-m', '--variant-mutation', dest='variant_mutation', type=str,
                        required=True, help='mutation residue, for example residue "A"')
    parser.add_argument('-p', '--variant-position', dest='variant_position', type=int,
                        required=True, help='variant position in sequence')

    parser.add_argument('-s', '--sequence',  type=str,
                        required=True, help='A protein primary sequence in a file in the format fasta.')
    parser.add_argument('--hhmfile', type=str,help='MSA files with hhm format from HHblits')
    # background mutation Features from G2s

    parser.add_argument('--pdbpath',  type=str,default='/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s1676_seq/PDB/',
                        required=True, help='The path for mutation-based structures')



    parser.add_argument('--dssppath', type=str,default='/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s1676_seq/dssp/',
                        required=True, help='The path for DSSP output files')
    parser.add_argument('--dsspbin', type=str, default='mkdssp',
                        required=True, help='DSSP binary executable.')
    parser.add_argument('--psiblastbin', type=str, default='psiblast',
                        required=True, help='psiblast binary executable.')
    parser.add_argument('--hhblitsbin', type=str, default='hhblits',help='Binary executable for computing hhblits profile, "hhblits" for fasta input file and "hhmake" for A3M,')
    parser.add_argument('--hhsuite', type=bool, default=False,help='HHM file')

    parser.add_argument('--psiblastout',  type=str,default='./Sequence/psiout',
                        required=True, help='psiblast output files, a path')
    parser.add_argument('--psiblastpssm',  type=str,default='./Sequence/pssmout',
                        required=True, help='A path for PSSM files')
    parser.add_argument('--psiblastdb',  type=str,default='/home/gongjianting/tools/PsiblastDB/swissprot',
                        required=True, help='background database for align in psiblast')


    parser.add_argument('--hhblitsdb',  type=str,default='/home/gongjianting/tools/HHsuitDB/UniRef30_2020_06', help='background database for align in tools hhblits')

    parser.add_argument('--hhblitsout',  type=str,default='./Sequence/hhblitout',
                        required=True, help='A path for hhblits output files')
    parser.add_argument('--hhblitshhm',  type=str,default='./Sequence/hhmout',
                        required=True, help='A path for storing HHM files')
    parser.add_argument("-v", "--version", action="version")
    parser.add_argument("-o", "--outfile", default='False', action='store_true',help='Whether save the result or not')
    parser.add_argument("--printout",  default='False', action='store_true', help='Whether print the result or not')
    parser.add_argument("--outfilepath", type=str,default='./input.npy',help='Output file path')
    parser.add_argument("-G","--G2s",  default='False', action='store_true', help='Fast Version, structures are unavailable in Q4')


Take mutation Q to H at position 104 of p53 protein as example. Its position in sequence is 9. So the run command as follows:

If provide MSA files with a3m format from HHblits. It can run with:
::

  python MU3DSP.py -p residueposition -w wildtyperesidue -m mutationresidue -o outfilename -s sequencepath --pdbpath pdbfilepath --dssppath dsspfilepath --dsspbin mkdssp-path --psiblastbin  psiblast-path --hhblitsbin hhblits-path --psiblastout psioutfile-path --psiblastpssm pssmoutfile-path --psiblastdb swissprot-path --hhblitshhm hhmoutfile-path --seqa3m MSA-a3m-file


For example:
::

    python MU3DSP.py -p 9 -w Q -m H -o True -s ./examples/SEQ/2ocj_A.fasta --pdbpath ./examples/PDBtest --dssppath ./examples/DSSPtest --dsspbin dssp --psiblastbin  psiblast --hhblitsbin hhmake --psiblastout ./examples/psiout --psiblastpssm ./examples/pssmout --psiblastdb ./PsiblastDB/swissprot --hhblitshhm ./examples/hhmout --outfilepath ./examples/2ocj_Q9H.npy --seqa3m ./examples/a3m/p53.a3m --printout True


If provide fasta files and HHblits is accessable.

::

  python MU3DSP.py -p residueposition -w wildtyperesidue -m mutationresidue -o outfilename -s sequencepath --pdbpath pdbfilepath --dssppath dsspfilepath --dsspbin mkdssp-path --psiblastbin  psiblast-path --hhblitsbin hhblits-path --psiblastout psioutfile-path --psiblastpssm pssmoutfile-path --psiblastdb swissprot-path --hhblitsdb UniRef30_2020_06-path --hhblitsout hhblitoutpath --hhblitshhm hhmoutfilepath
