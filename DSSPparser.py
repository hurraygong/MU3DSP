#coding:utf-8 
#

#提取出PDBID和链的ID,突变位点以及位置

import os
import pickle

#根据DSSP文件，提取出特征，可选则，
{'Wilke': { 
         'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 
         'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 
         'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 
        'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 
         'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0 
     }}
class DSSPparser():
    def file_name(self,file_dir):   
        L=[]   
        for dirpath, dirnames, filenames in os.walk(file_dir):  
            #print(dirpath, dirnames, filenames)
            for file in filenames :  
                if os.path.splitext(file)[1] == '.dssp':  
                    L.append(os.path.join(dirpath, file))  
        return L 
      
    def dssp_paser(self,filepath,withpath=False,filepathinfo = True ):
             
        if filepathinfo == True:
            files = self.file_name(filepath) 
        else:
            files = [filepath]
            
        alldssp = {}
        for filename in files:
            onedsspdict = self.onedssp_parser(filename)
            PDB = os.path.split(filename)[-1].split('.')[0]
            alldssp.update({PDB:onedsspdict})
        return alldssp 
                
    
    def onedssp_parser(self,filename) :
        
        lines = open(filename,'r').readlines()

        onefiledssp = {} 
        featureinfoindex = 0
        feainfo = False
        for lineindex in range(len(lines)):

            eachline = lines[lineindex]
            
            if eachline.split('\n')[0] != '':
                
                if eachline.split()[0] == '#':
                    featureinfoindex = lineindex
                    feainfo = True             
                if feainfo == True:
                    
                    if lineindex > featureinfoindex:                       
                        if eachline[5:10].strip() != '':
                            chain,pos_AA, fea = self.perline(eachline) 
                            if chain not in onefiledssp.keys():
                                onefiledssp.update({chain:{pos_AA:fea}})
                            else:
                                onefiledssp[chain].update({pos_AA:fea})                          
                        
                else:
                    pass   
            else:
                pass
        #print(onefiledssp)
        return onefiledssp
        
    def perline(self,eachline):
        
        chain = eachline[11]
        pos_AA = eachline[5:10].strip() + '_' + eachline[13]

        if eachline[16] == ' ':
            ss = 'C'
        else:
            ss = eachline[16]
        ASA = float(eachline[34:38].strip())
        
        kappa = float(eachline[91:97].strip())      
        alpha = float(eachline[97:103].strip())
        phi  = float(eachline[103:109].strip())
        psi = float(eachline[109:115].strip())
        
        xca = float(eachline[115:122].strip())
        yca = float(eachline[122:129].strip())
        zca = float(eachline[129:136].strip())
        
        fea = [ss,ASA,kappa,alpha,phi,psi,xca,yca,zca]   
        
        return chain,pos_AA, fea
        
if __name__ == '__main__':
    
    #targetpath = 'C:\\Myself\\GNN\\fromwang\\ease-mm\\DSSPfiles' #for train 1676
    targetpath = 'E:\\stability\\Sta_Dssp\\PDBDssp236'#for test 236
    
    dssp = DSSPparser()
    alldsspparsed = dssp.dssp_paser(targetpath,withpath=False)
    print(alldsspparsed)

    pickle.dump(alldsspparsed, open('E:\\stability\\Sta_Dssp\\263stability_dsspfea.pic','wb'), protocol=2)
    #pickle.dump(alldsspparsed, open('C:\\Myself\\GNN\\fromwang\\ease-mm\\DSSPfiles\\263stability_dsspfea.pic','wb'), protocol=2)
    




