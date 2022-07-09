Quickstart
----------

Quikly Run scGNN
^^^^^^^^^^^^^^^^
This is the parameters to run SCpre-seq.

::

    parser.add_argument('-l', '--variant-site', dest='variant_list', type=str,required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')

    parser.add_argument('-w', '--variant-wildtype', dest='variant_wildtype', type=str,
                        required=True, help='wild-type residue, for example residue "A"')
    parser.add_argument('-m', '--variant-mutation', dest='variant_mutation', type=str,
                        required=True, help='mutation residue, for example residue "A"')
    parser.add_argument('-p', '--variant-position', dest='variant_position', type=int,
                        required=True, help='variant position in sequence')

    parser.add_argument('-s', '--sequence',  type=str,
                        required=True, help='A protein primary sequence in a file in the format fasta.')

    # background mutation Features from G2s

    parser.add_argument('--pdbpath',  type=str,default='/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s1676_seq/PDB/',
                        required=True, help='A path for storage matched PDB structures')
    parser.add_argument('--dssppath', type=str,default='/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s1676_seq/dssp/',
                        required=True, help='A path for storage dssp output files')
    parser.add_argument('--dsspbin',  type=str,default='mkdssp',
                        required=True, help='Execute bin path for DSSP')
    parser.add_argument('--psiblastbin',  type=str,default='psiblast',
                        required=True, help='Execute bin path for psiblast')
    parser.add_argument('--hhblitsbin',  type=str,default='hhblits',
                        required=True, help='Execute bin path for hhblits')

    parser.add_argument('--psiblastout',  type=str,default='./Sequence/psiout',
                        required=True, help='A path for storage psiblast out files')
    parser.add_argument('--psiblastpssm',  type=str,default='./Sequence/pssmout',
                        required=True, help='A path for storage psiblast pssm files')
    parser.add_argument('--psiblastdb',  type=str,default='/home/gongjianting/tools/PsiblastDB/swissprot',
                        required=True, help='background database for align in psiblast')


    parser.add_argument('--hhblitsdb',  type=str,default='/home/gongjianting/tools/HHsuitDB/UniRef30_2020_06',
                        required=True, help='background database for align in tools hhblits')

    parser.add_argument('--hhblitsout',  type=str,default='./Sequence/hhblitout',
                        required=True, help='A path for storage hhblits hhr files')
    parser.add_argument('--hhblitshhm',  type=str,default='./Sequence/hhmout',
                        required=True, help='A path for storage hhblits hhm files')
    parser.add_argument("-v", "--version", action="version")
    parser.add_argument("-o", "--outfile",type=bool, default=False,help='Whether save the result or not')
    parser.add_argument("-printout", type=bool, default=True, help='Whether print the result or not')
    parser.add_argument("-outpath", "--outfilepath", type=str,default='./',help='Output file path')


Take mutation Q to H at position 104 of p53 protein as example. Its position in sequence is 9. So the run command as follows:


::

     python StatGetPDBlist.py -p 9 -w Q -m H -o True -s ./examples/SEQ/2ocj_A129D.fasta --pdbpath ./examples/PDBtest --dssppath ./examples/DSSPtest --dsspbin mkdssp --psiblastbin  psiblast --hhblitsbin hhblits --psiblastout ./examples/psiout --psiblastpssm ./examples/pssmout --psiblastdb /home/gongjianting/tools/PsiblastDB/swissprot --hhblitsdb /home/gongjianting/tools/HHsuitDB/UniRef30_2020_06 --hhblitsout ./examples/hhblitout --hhblitshhm ./examples/hhmout -outpath ./examples/