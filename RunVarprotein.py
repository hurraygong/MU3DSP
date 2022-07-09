
import time,os


lines = open('/home/gongjianting/gongjt/SCpre-seq/PDB2Seq/varbench_seq.txt','r').readlines()

for eachline in lines:
    if eachline[0] == '>':
        fastafilename = eachline[1:].split('_')[0]+eachline.split('_')[1]+'_'+eachline.split('_')[2]
        seqpath ='/home/gongjianting/gongjt/SCpre-seq/PDB2Seq/VarbenchFasta/'+fastafilename+'.fasta'
        fastaseq = open(seqpath,'r').readlines()
        seq = ''
        for subseq in fastaseq[1:]:
            seq = seq + subseq.strip('\n').strip()
        i = int(eachline.split('_')[2][1:-1])
        Aacid = eachline.split('_')[2][-1]
        cmd = "python StatGetPDBlist.py -p "+  str(i)+ " -w "+ seq[i-1] + " -m " + Aacid + " -o " + str(i)+Aacid  + " -s "+seqpath+" --pdbpath ./PDB2Seq/VarbenchPDB/ --dssppath ./PDB2Seq/DSSPtest --dsspbin mkdssp --psiblastbin  psiblast --hhblitsbin hhblits --psiblastout ./Sequence/psiout --psiblastpssm ./VarSequence/pssmout --psiblastdb /home/gongjianting/tools/PsiblastDB/swissprot --hhblitsdb /home/gongjianting/tools/HHsuitDB/UniRef30_2020_06 --hhblitsout ./VarSequence/hhblitout --hhblitshhm ./VarSequence/hhmout"
        os.system(cmd)