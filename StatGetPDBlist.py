# Python 3.6
# A Simple Python template program reading G2S APIs
# download PDBfiles and thier dssp files based on

import requests
import json
import pickle
import wget
import os
import argparse
import time
from argparse import ArgumentParser
from DSSPparser import *
import sys
import numpy as np
from collections import Counter
from HHMwm import ParsePssm,ParseHHm
import pandas as pd
import lightgbm as lgb
from sklearn import metrics
from sklearn.metrics import *
from scipy import stats
# Setup API, e.g. https://g2s.genomenexus.org/api/alignments/uniprot/P02185/residueMapping?positionList=10
# Can change to any valid URL in G2S API, UniprotID or positionist
############################### DSSP #################################
def get_global_features(examples, alldsspparsed,msite,singlePDB=True):

    siglefea = []
    for example in examples:
        fea = get_siglesite(example, alldsspparsed)
        if fea != []:
            siglefea.append(fea)

    if siglefea == []:
        return [],[],[],[],[]
    else:
        if singlePDB!=True:
            SS_type, asa, mean_ang = Fea_parse(siglefea, msite)
        else:
            SS_type, asa, mean_ang = onePDBFea_parse(siglefea, msite)

    if SS_type == []:
        return [],[],[],[],[]
    else:
        return SS_type, asa[0],asa[1],asa[2], mean_ang

def onePDBFea_parse(siglefea, msite):

    maxASA = {
        'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0,
        'C': 148.0, 'Q': 214.0, 'E': 214.0, 'G': 97.0,
        'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
        'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0,
        'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0
    }
    ss = siglefea[0][0]
    rasa = float(siglefea[0][1]) / maxASA[msite]
    angle = siglefea[0][2:-3]
    if rasa > 0.25:
        ex = 1
    else:
        ex = 0
    return ss, [rasa,ex,rasa],angle


def Fea_parse(siglefea, msite):

    maxASA = {
        'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0,
        'C': 148.0, 'Q': 214.0, 'E': 214.0, 'G': 97.0,
        'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
        'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0,
        'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0
    }
    ss = []
    rasa = []
    angle = []
    for each in siglefea:
        #print(each,msite)
        ss.append(each[0])
        rasa.append(float(each[1]) / maxASA[msite])
        angle.append(each[2:-3])

    e_count = 0
    b_count = 0
    e_num = 0
    b_num = 0
    for ea in rasa:
        if ea > 0.25:
            e_count += 1
            e_num += ea
        else:
            b_count += 1
            b_num += ea

    muta_all = (b_num + e_num) / (e_count + b_count)
    if e_count >= b_count:
        b = e_num / e_count
        muta_f = 1
    else:
        b = b_num / b_count
        muta_f = 0

    asa = [muta_all, muta_f, b]
    SS_type = Counter(ss).most_common(1)[0][0]
    mean_ang = np.mean(np.array(angle), axis=0)
    return SS_type, asa, mean_ang


def get_siglesite(example, alldsspparsed):
    pdbid = example[0].split('_')[0]
    pdbchain = example[0].split('_')[1]
    mutasite = example[-1]
    mutasitekey = str(example[1]) + '_' + mutasite
    mutafea = []
    if pdbid in alldsspparsed.keys():
        if pdbchain in alldsspparsed[pdbid].keys():
            if mutasitekey in alldsspparsed[pdbid][pdbchain].keys():
                mutafea = alldsspparsed[pdbid][pdbchain][mutasitekey]
            else:
                mutafea = []
    else:
        pass
    return mutafea

def get_imprecise_match(example, alldsspparsed):
    #print('###',example)
    pdbid = example[0].split('_')[0]
    pdbchain = example[0].split('_')[1]
    mutasite = example[-1]
    mutafea = []
    if pdbid in alldsspparsed.keys():
        if pdbchain in alldsspparsed[pdbid].keys():
            for eachsitekey in alldsspparsed[pdbid][pdbchain].keys():
                if eachsitekey.split('_')[0] == str(example[1]):
                    mutafea= alldsspparsed[pdbid][pdbchain][eachsitekey]
                    break
    #print(mutafea)
    return mutafea
############################### DSSP #################################


############################### AAPS #################################

def global_dist(example,alldsspparsed,dssppath=None,precise_match=True):
    pdb = example[0].split('_')[0]
    pdbchain = example[0].split('_')[1]
    mutasite = example[-1]
    mutasitekey = str(example[1]) + '_' + mutasite
    if  dssppath == None:
        pdbid =pdb
    else:
        pdbid= os.path.join(dssppath,pdb)
    cal_dist = True

    muxca,muyca,muzca = None,None,None
    if pdbid in alldsspparsed.keys():
        if pdbchain in alldsspparsed[pdbid].keys():
            if precise_match==True:
                if mutasitekey in alldsspparsed[pdbid][pdbchain].keys():
                    sitedssp = alldsspparsed[pdbid][pdbchain][mutasitekey]
                    muxca = sitedssp[-3]
                    muyca = sitedssp[-2]
                    muzca = sitedssp[-1]
                else:
                    cal_dist = False

            else:
                for eachsitekey in alldsspparsed[pdbid][pdbchain].keys():
                    if eachsitekey.split('_')[0] == str(example[1]):
                        sitedssp = alldsspparsed[pdbid][pdbchain][eachsitekey]
                        muxca = sitedssp[-3]
                        muyca = sitedssp[-2]
                        muzca = sitedssp[-1]
                        break
        else:
            cal_dist = False

        #check unformat residues in PDB files, such as missing residues
        if (muxca == None and muyca == None and muzca==None):
            resi_exit = False
        else:
            resi_exit = True


        mutadict1 = {}
        if cal_dist == False or resi_exit == False:
            mutad2 = []
            pass
        else:
            for eachsitekey in alldsspparsed[pdbid][pdbchain].keys():
                Xfea = alldsspparsed[pdbid][pdbchain][eachsitekey]
                a = abs(float(Xfea[-3]) - float(muxca))
                b = abs(float(Xfea[-2]) - float(muyca))
                c = abs(float(Xfea[-1]) - float(muzca))
                distance = (a * a + b * b + c * c) ** 0.5
                distance1 = round(distance, 4)
                # list1.append(distance1)
                mutadict1.update({eachsitekey: distance1})
                mutad2 = sorted(mutadict1.items(), key=lambda item: item[1])
    else:
        mutad2 = []

    return mutad2


def get_feamatrix(distencelist, r=10):
    pe = {'A': 0, 'V': 0, 'L': 0, 'I': 0, 'M': 0, 'C': 0, 'F': 0, 'W': 0, 'Y': 0, 'H': 0, 'S': 0, 'T': 0, 'N': 0,
          'Q': 0, 'K': 0, 'R': 0, 'D': 0, 'E': 0, 'G': 0, 'P': 0, 'X': 0}
    aa_count = 0
    aa_distence = []
    for z in range(1, len(distencelist)):
        eachitem = distencelist[z]
        if eachitem[1] <= r:
            aa_count = aa_count + 1
            # print(eachitem[0])
            try:
                pe[eachitem[0].split('_')[1]] = pe[eachitem[0].split('_')[1]] + 1
            except:
                pe['X'] = pe['X'] + 1

    for key in ['A', 'V', 'L', 'I', 'M', 'C', 'F', 'W', 'Y', 'H', 'S', 'T', 'N',
          'Q', 'K', 'R', 'D', 'E', 'G', 'P']:
    #for key in pe.keys():
        aa_distence.append(round(pe[key] / aa_count, 4))

    return aa_distence

def get_str_features(examples, alldsspparsed, dssppath=None,rd=10,singlePDB=True):
    mutafea = []
    for example in examples:
        dis = global_dist(example, alldsspparsed,dssppath)
        if dis == []:
            pass
        else:
            muta = get_feamatrix(dis, r=rd)
            mutafea.append(muta)
        # mutafea.append(muta)
    if mutafea == []:
            meanfea = []
        #pass
    else:
        if singlePDB!=True:
            meanfea = np.nanmean(np.array(mutafea), axis=0)
        else:
            meanfea = np.array(mutafea[0])

    return mutafea, meanfea
############################### AAPS #################################

def file_name(file_dir, fileformat):
    # fileformat : like .pdb or .dssp
    L = []
    for dirpath, dirnames, filenames in os.walk(file_dir):
        # print(dirpath, dirnames, filenames)
        for file in filenames:
            if os.path.splitext(file)[1] == fileformat:
                L.append(os.path.join(dirpath, file))
    return L

def Open_pickle_file(picklefile):
    with open(picklefile, 'rb') as f:
        content = pickle.load(f)
    return content

def Get_wildsitepdb(wild_key, maplist):
    pdblist = []
    for mappdb in maplist:
        if mappdb['pdbaAminoacid'] == wild_key:
            eachpdb = [mappdb['name'], mappdb['pdbPosition'], mappdb['pdbaAminoacid']]
            pdblist.append(eachpdb)
    return pdblist

def G2s_mapsite_feature(wm,mutabasepath = r"./DsspFea/"):
    mutasite380_ss_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_ss_fea.pic'), 'rb'))
    mutasite380_asamean_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_asamean_fea.pic'), 'rb'))
    mutasite380_asaall_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_asaall_fea.pic'), 'rb'))
    mutasite380_asa01_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_asa01_fea.pic'), 'rb'))
    mutasite380_angle_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_angle_fea.pic'), 'rb'))
    return mutasite380_ss_fea[wm],mutasite380_asaall_fea[wm],mutasite380_asa01_fea[wm],mutasite380_asamean_fea[wm],mutasite380_angle_fea[wm]

def G2s_mapsite_strfeature(wm,mutabasepath = r"./StructureFea/",rd=20):
    mutasite380_str_fea = pickle.load(
        open(os.path.join(mutabasepath,'mutasite380_str_disfea'+str(rd)+'.pic'), 'rb'))
    return mutasite380_str_fea[wm]

def MutationBasedFea(ResID):
    wildss, wildasa_all, wildasa_01, wildasa_me, wildangle = G2s_mapsite_feature(ResID)
    wildstr = G2s_mapsite_strfeature(ResID)
    return wildss, wildasa_all, wildasa_01, wildasa_me, wildangle,wildstr


def MuStructureFea(wildsites,mutasites,alldsspparsed,muta,wild,singlePDB=True):
    if wildsites == [] and mutasites == []:
        wildss, wildasa_all, wildasa_01, wildasa_me, wildangle, wildstr = MutationBasedFea(muta + wild)
        mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle, mutastr = MutationBasedFea(wild + muta)
    else:
        # ddg.append(float(site.split()[-3]))
        if wildsites == []:
            wildss, wildasa_all, wildasa_01, wildasa_me, wildangle, wildstr = MutationBasedFea(muta + wild)
        else:
            wildss, wildasa_all, wildasa_01, wildasa_me, wildangle = get_global_features(wildsites, alldsspparsed, wild,
                                                                                         singlePDB=singlePDB)
            if wildss == []:
                # print('**', site)
                wildss, wildasa_all, wildasa_01, wildasa_me, wildangle, wildstr = MutationBasedFea(
                    muta + wild)
            else:
                _, wildstr = get_str_features(wildsites, alldsspparsed, rd=20, singlePDB=singlePDB)

        if mutasites == []:
            mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle, mutastr = MutationBasedFea(wild + muta)
        else:
            mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle = get_global_features(mutasites, alldsspparsed,
                                                                                         muta, singlePDB=singlePDB)

            if mutass == []:
                # print('##', site)
                mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle, mutastr = MutationBasedFea(
                    wild + muta)
            else:
                _, mutastr = get_str_features(mutasites, alldsspparsed, rd=20, singlePDB=singlePDB)
    return wildss, wildasa_all, wildasa_01, wildasa_me, wildangle, wildstr,mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle, mutastr

def SecondStrucType(X,featype = 1):
    ssdict = {0: 0, 'H': 1, 'B': 2, 'E': 3, 'G': 4, 'I': 5, 'T': 6, 'S': 7, 'C': 0}
    #H=HGI, E=EB, C=STC
    ssdict3 = {0: 0, 'H': 1, 'B': 2, 'E': 2, 'G': 1, 'I': 1, 'T': 0, 'S': 0, 'C': 0}
    if featype == 1:
        zzz1 = ssdict[X]
    elif featype == 3:
        zzz1 = np.zeros(3)
        zzz1[ssdict3[X]] = 1
    elif featype == 8:
        zzz1= np.zeros(8)
        zzz1[ssdict[X]] = 1
    # print(ssdict[SS_wild_type])
    else:
        print('error')
    return zzz1

def  MuSecondStructure(wildss,mutass,inte='cancat',featype=1):
    w = SecondStrucType(wildss, featype)
    m = SecondStrucType(mutass, featype)
    if inte=='minus':
        ssfea = m-w
        if featype == 1:
            return [ssfea]
        else:
            return ssfea.tolist()
    elif inte == 'cancat':
        if featype == 1:
            return [w,m]
        else:
            return w.tolist()+m.tolist()

def StrAngle(wildangle,mutaangle,AngleFeaType='PolarC'):
    # angle

    #Polar coordinates = PolarC
    if AngleFeaType=='PolarC':
        angle = []
        angle_x = np.cos(wildangle)
        angle_y = np.sin(wildangle)

        angle_mu_x = np.cos(mutaangle)
        angle_mu_y = np.sin(mutaangle)

        for i in range(0, 4):
            if wildangle[i] <= 0:
                angle.append(
                    -((angle_mu_x[i] - angle_x[i]) ** 2 + (angle_mu_y[i] - angle_y[i]) ** 2) ** 0.5)
            else:
                angle.append(((angle_mu_x[i] - angle_x[i]) ** 2 + (angle_mu_y[i] - angle_y[i]) ** 2) ** 0.5)
    elif AngleFeaType=='minus':
        angle = (np.array(mutaangle) - np.array(wildangle)).tolist()
    elif AngleFeaType=='cancat':
        angle = np.array(mutaangle).tolist() + np.array(wildangle).tolist()
    return angle


def MuRasa(mutaasa_me, wildasa_me, rasatype='minus'):
    if rasatype == 'minus':
        dRasa = [mutaasa_me - wildasa_me]
    elif rasatype == 'cancat':
        dRasa = [mutaasa_me, wildasa_me]
    return dRasa


def MuStr(wildstr, mutastr, strtype='minus'):
    if strtype == 'minus':
        AAfs = wildstr - mutastr
    elif rasatype == 'cancat':
        AAfs = [wildstr.tolist(), mutastr.tolist()]
    return AAfs


def DsspFea(wildss, wildasa_me, wildangle, wildstr, mutass,
            mutaasa_me, mutaangle, mutastr, AngleFeaType='PolarC', inte='minus', rasatype='minus',strtype='minus',featype=1):
    Rsa = StrAngle(wildangle, mutaangle, AngleFeaType=AngleFeaType)
    Rss = MuSecondStructure(wildss, mutass, inte=inte, featype=featype)
    Rrasa = MuRasa(mutaasa_me, wildasa_me, rasatype=rasatype)
    AApsstr = MuStr(wildstr, mutastr, strtype=strtype)
    return Rsa, Rss, Rrasa, AApsstr
##################################################### str fea ##############################

def Get_download_url(ID, position):
    hostName = "https://g2s.genomenexus.org/api/"
    apiName = "alignments"
    id_type = 'uniprot'
    apiType = "residueMapping?positionList="
    url = hostName + apiName + "/" + id_type + '/' + ID + "/" + apiType + position
    return url
    # Set up API URL

# https://g2s.genomenexus.org/api/alignments/residueMapping?sequence=RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVGGGCRAKRNNFKSAEDCMRTCGGA&positionList=33%2C34
def Get_download_from_seq(seq, position):
    hostName = "https://g2s.genomenexus.org/api/"
    apiName = "alignments"
    apiType = "residueMapping?sequence="
    positionlist = "&positionList="
    url = hostName + apiName + "/" + apiType + seq + positionlist + str(position)
    return url
    # Set up API URL


def Get_jData(url):
    # Request API
    myResponse = requests.get(url)
    # print (myResponse.status_code)

    # For successful API call, response code will be 200 (OK)
    if (myResponse.ok):
        # Loads (Load String) takes a Json file and converts into python data structure (dict or list, depending on JSON)
        # In this Example, jData are lists of Residues from Genome Mapping to Protein Structures
        jData = json.loads(myResponse.content)
        # pring output
        return jData
    else:
        # If response code is not ok (200), print the resulting http error code with description
        myResponse.raise_for_status()


def Creat_Unp_PDBmap_path(jsondata):
    PDB_id_list = []
    chain_position_list = []

    for each in jsondata:

        pdbid = each['pdbId']
        PDB_id_list.append(pdbid)
        # test position_list number
        for eachposition in each['residueMapping']:
            onepdb = {'pdbid': pdbid, 'name': each['pdbNo'], 'chain': each['chain'],
                      'pdbPosition': eachposition['pdbPosition'], 'pdbaAminoacid': eachposition['pdbAminoAcid'],
                      "pdbFrom": each["pdbFrom"], "seqFrom": each["seqFrom"]}
            if onepdb not in chain_position_list:
                chain_position_list.append(onepdb)
            else:
                pass
    return set(PDB_id_list), chain_position_list


def setlocalpath(path):
    path = path.strip()
    path = path.rstrip("/")

    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)
        # print (path+' created successful')
    else:
        # print (path+' dir exsit')
        pass


# get UNP ID and Position from datasets file of ease-mm
def Get_UnpID_Pos(filename):
    Unp_Id_list = open(filename).readlines()
    entrylist = []

    for line in Unp_Id_list:
        entry = line.split()
        if len(entry) > 2:
            if line[0] != "#":
                entrylist.append([entry[2], entry[3], entry[-3]])
            else:
                pass
        else:
            pass
    # pickle.dumps(entrylist, protocol=1)
    return entrylist


class InputError(Exception):
    def __init__(self,ErrorInfo):
        super().__init__(self) #初始化父类
        self.errorinfo=ErrorInfo
    def __str__(self):
        return self.errorinfo


def checkargs(args):

    if args.variant_wildtype not in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                       'W', 'Y']:
        print('input wildtype residue is not a standard amino acid, please check the input!')
        raise ValueError
    if args.variant_mutation not in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                       'W', 'Y']:
        print('input mutation residue is not a standard amino acid, please check the input')
        raise ValueError
    if seq[varpos - 1] != wildres:
        print('input variant position is not matched to input protein sequence, please check the input')
        raise ValueError

    return
if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Download PDBfiles from G2s & get features')
    # parser.add_argument('-l', '--variant-site', dest='variant_list', type=str,
    #                     required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')

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
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')



    parser.add_argument('--dssppath', type=str,default='/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s1676_seq/dssp/',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--dsspbin',  type=str,default='mkdssp',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--psiblastbin',  type=str,default='psiblast',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--hhblitsbin',  type=str,default='hhblits',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')

    parser.add_argument('--psiblastout',  type=str,default='./Sequence/psiout',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--psiblastpssm',  type=str,default='./Sequence/pssmout',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--psiblastdb',  type=str,default='/home/gongjianting/tools/PsiblastDB/swissprot',
                        required=True, help='background database for align in psiblast')


    parser.add_argument('--hhblitsdb',  type=str,default='/home/gongjianting/tools/HHsuitDB/UniRef30_2020_06',
                        required=True, help='background database for align in tools hhblits')

    parser.add_argument('--hhblitsout',  type=str,default='./Sequence/hhblitout',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument('--hhblitshhm',  type=str,default='./Sequence/hhmout',
                        required=True, help='A list of variants, one per line in the format "POS WT MUT", a file')
    parser.add_argument("-v", "--version", action="version")
    parser.add_argument("-o", "--outfile",type=bool, default=False,help='Whether save the result or not')
    parser.add_argument("-printout", type=bool, default=True, help='Whether print the result or not')
    parser.add_argument("-outpath", "--outfilepath", type=str,default='./',help='Output file path')

    #https: // xgxm.xueguoxue.com /  # /user/receiveLearnCard?cardId=d7ce3f00264d23

    args = parser.parse_args()

    #check input variant
    wildres = args.variant_wildtype
    mutares = args.variant_mutation
    varpos = args.variant_position

    mutaseq = args.sequence

    #fasta path
    fastaseq = open(mutaseq, 'r').readlines()
    fastaid = os.path.splitext(os.path.split(mutaseq)[-1])[0]

    seq = ''
    for subseq in fastaseq[1:]:
        seq = seq + subseq.strip('\n').strip()

    checkargs(args)
    sitesmapping = {}

    ############# run psiblast #####################
    if not os.path.exists(args.psiblastout):
        os.makedirs(args.psiblastout)
    else:
        pass

    if not os.path.exists(args.psiblastpssm):
        os.makedirs(args.psiblastpssm)
    else:
        pass

    blastout = os.path.join(args.psiblastout,fastaid)
    blastpssmout = os.path.join(args.psiblastpssm,fastaid)
    if not os.path.exists(blastpssmout+'.pssm'):
        blastcmd = args.psiblastbin+' -comp_based_stats 1 -evalue 0.001 -num_iterations 3 -db ' + args.psiblastdb +' -query ' +  mutaseq + ' -out ' + blastout + '.out -out_ascii_pssm ' + blastpssmout + '.pssm -num_threads 18'
        os.system(blastcmd)
    else:
        pass


############# run hhblits #####################

    if not os.path.exists(args.hhblitsout):
        os.makedirs(args.hhblitsout)
    else:
        pass

    if not os.path.exists(args.hhblitshhm):
        os.makedirs(args.hhblitshhm)
    else:
        pass
    hhblits_out = os.path.join(args.hhblitsout,fastaid)
    hhblitshhmout = os.path.join(args.hhblitshhm,fastaid)
    if not os.path.exists(hhblitshhmout+'.hhm'):
        hhblitscmd = args.hhblitsbin + ' -i ' + mutaseq + ' -d ' + args.hhblitsdb +' -o ' + hhblits_out  + '.hhr -ohhm ' + hhblitshhmout + '.hhm -e 1e-3 -n 2 -p 20 -Z 250 -z 1 -b 1 -B 250'
        os.system(hhblitscmd)
    else:
        pass

################ download mapped pdb files ###############
    eachurl = Get_download_from_seq(seq, varpos)

    try:
        eachjdata = Get_jData(eachurl)
    except:
        time.sleep(5)
        eachjdata = Get_jData(eachurl)

    PDB_id_list, chain_position_list = Creat_Unp_PDBmap_path(eachjdata)

    # dssppath
    unp_wpm = fastaid + '_' + wildres + str(varpos) + mutares
    sitesmapping.update({unp_wpm: chain_position_list})
    #print(sitesmapping)

    #checking PDBpath

    if not os.path.exists(args.pdbpath):
        os.makedirs(args.pdbpath)
    else:
        pass
    # checking DSSP path
    if not os.path.exists(args.dssppath):
        os.makedirs(args.dssppath)
    else:
        pass

    #print(args.pdbpath,args.dssppath)

    for eachentry in PDB_id_list:
        localpath = os.path.join(args.pdbpath,eachentry.lower() + '.pdb')
        if not os.path.exists(localpath):
            originpath = os.path.join('https://files.rcsb.org/download/', eachentry.lower() + '.pdb')
            try:
                wget.download(originpath, localpath)
            except:
                print('download error ' + localpath)
            # download PDBfiles
            try:
                onedssppath = os.path.join(args.dssppath,eachentry.lower() + '.dssp')
                cmd = args.dsspbin + ' -i ' + localpath + ' -o ' + onedssppath
                os.system(cmd)
            except:
                print('Dssp error ' + localpath)
        else:
            pass

    #find features and rm the results
    dssp = DSSPparser()
    #print(dssp)
    alldsspparsed = dssp.dssp_paser(args.dssppath)

    # Dssp entrysites
    # SS+RASA
    wildsites = Get_wildsitepdb(wildres, sitesmapping[unp_wpm])
    mutasites = Get_wildsitepdb(mutares, sitesmapping[unp_wpm])

    wildss, wildasa_all, wildasa_01, wildasa_me, wildangle, wildstr, mutass, mutaasa_all, mutaasa_01, mutaasa_me, mutaangle, mutastr = MuStructureFea(
        wildsites, mutasites, alldsspparsed, mutares, wildres, singlePDB=False)

    Rsa, Rss, Rrasa, AApsstr = DsspFea(wildss, wildasa_me, wildangle, wildstr,
                                       mutass, mutaasa_me, mutaangle, mutastr, AngleFeaType='PolarC', inte='minus',
                                       rasatype='cancat', strtype='minus', featype=3)
    # print(Rsa, Rss, Rrasa, AApsstr)
    dsspfea = Rss + Rrasa
    AApsstr = AApsstr.tolist()

    # AAP
    with open(r'./AAP/parameters_list.json', 'r') as pic:
        parameters = json.load(pic)
    DDpara = []
    for i in range(0, len(parameters[mutares])):
        DDpara.append(round(parameters[mutares][i] - parameters[wildres][i], 3))

    ############################# PSSM
    Pssm = []
    PSSMtitle = {'A': 0, 'V': 19, 'L': 10, 'I': 9, 'M': 12, 'C': 4, 'F': 13, 'W': 17, 'Y': 18, 'H': 8, 'S': 15, 'T': 16,
                 'N': 2, 'Q': 5, 'K': 3, 'R': 1, 'D': 3, 'E': 6, 'G': 7, 'P': 14}
    pssm = ParsePssm()
    siteseqPssm = pssm.getpssmfeature(args.psiblastpssm)

    oripssm = [1 / (1 + (2.7182) ** (-float(x))) for x in siteseqPssm[fastaid][int(varpos) - 1]]

    wildpssm = oripssm[PSSMtitle[wildres]]
    mutapssm = oripssm[PSSMtitle[mutares]]
    DDpssm = [mutapssm - wildpssm]

    Pssm = oripssm + DDpssm

    ############################## HHM
    HHM = []
    hhm = ParseHHm()
    siteseqhhm = hhm.gethmmfeature(args.hhblitshhm)
    hhmlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    hhmdict = {each: hhmlist.index(each) for each in hhmlist}

    wildhhm = float(siteseqhhm[fastaid][int(varpos) - 1][hhmdict[wildres]])
    mutahhm = float(siteseqhhm[fastaid][int(varpos) - 1][hhmdict[mutares]])
    orihhm = [float(x) for x in siteseqhhm[fastaid][int(varpos) - 1]]
    DDhhm = [mutahhm - wildhhm]
    HHM = orihhm + DDhhm

    preFeatures = dsspfea + AApsstr + DDpara + Pssm  + HHM
    bst = lgb.Booster(model_file='model_rmse_0.664_0.513_0920.txt')
    preds_online = bst.predict([preFeatures], num_iteration=bst.best_iteration)  # 输出概率
    args.outfilename = fastaid + '_'+ wildres+str(varpos)+mutares
    if args.outfile == True:
        outfile = os.path.join(args.outfilepath, args.outfilename + '.npy')
        np.save(outfile, preds_online[0])
    if args.printout == True:
        print(wildres+str(varpos)+mutares, round(preds_online[0],4))

