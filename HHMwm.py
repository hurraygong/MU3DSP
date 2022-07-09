import os, json
import sys
import numpy as np
import pickle


class ParsePssm():

    def _loadFile(self, file):
        ''' Open the file with the path
            @return: The lines of the file
            @param file: The full path of the file, str
        '''
        try:
            with open(file) as fh:
                filelines = fh.readlines()  # read file and get all lines
            return filelines
        except IOError as err:
            print('File error: ' + err)

    def __parseonefilepssm(self, filelines):

        pssmvalue = []
        for line in filelines:
            if len(line.split()) == 44:
                eachitem = line.split()[2:22]
                pssmvalue.append(eachitem)
        return pssmvalue

    def getpssmfeature(self, filepath):

        fileslist = os.listdir(filepath)
        allentrypssm = {}
        for onefile in fileslist:
            path = filepath + '/' + onefile
            filelines = self._loadFile(path)

            pssmvalue = self.__parseonefilepssm(filelines)

            uniprotID = onefile.split('/')[-1].split('.')[0]

            uniprotpssm = {uniprotID: pssmvalue}
            allentrypssm.update(uniprotpssm)
        return allentrypssm


def cal_pssm():
    PSSMtitle = {'A': 0, 'V': 19, 'L': 10, 'I': 9, 'M': 12, 'C': 4, 'F': 13, 'W': 17, 'Y': 18, 'H': 8, 'S': 15, 'T': 16,
                 'N': 2, 'Q': 5, 'K': 3, 'R': 1, 'D': 3, 'E': 6, 'G': 7, 'P': 14}

    # sites253 = pickle.load(open('/u1/home/jghhd/stability_change/Reimplement/253mutation_sites.pickle','rb'))

    # print(examples)
    # lines = open('/u1/home/jghhd/stability_change/Get_DSSPfile_Code/Sep_parse/motify_dataset.out').readlines()
    # lines = open('/u1/home/jghhd/stability_change/wild_feature/Win_seq/onesite.out').readlines()
    lines = open('/storage/htc/joshilab/jghhd/SC/stability_change/Data_pro/S1676').readlines()

    filepath = "/storage/htc/joshilab/jghhd/SC/stability_change/Data_pro/PSSMout/"
    pssm = ParsePssm()
    siteseqPssm = pssm.getpssmfeature(filepath)
    # print(siteseqPssm)
    PSSMwm = []
    PSSMwm21 = []
    for i in range(1, len(lines)):
        line = lines[i].split()
        unpid = line[2].upper()
        wild = line[3][0]
        muta = line[3][-1]
        position = line[3][1:-1]
        try:
            wildpssm = 1 / (1 + (2.7182) ** (-float(siteseqPssm[unpid][int(position) - 1][PSSMtitle[wild]])))
            mutapssm = 1 / (1 + (2.7182) ** (-float(siteseqPssm[unpid][int(position) - 1][PSSMtitle[muta]])))
            DDpssm = round((mutapssm - wildpssm), 5)
            pssm = [round(1 / (1 + (2.7182) ** (-float(x))), 5) for x in siteseqPssm[unpid][int(position) - 1]]
        except:
            unpid = line[2]
            wildpssm = 1 / (1 + (2.7182) ** (-float(siteseqPssm[unpid][int(position) - 1][PSSMtitle[wild]])))
            mutapssm = 1 / (1 + (2.7182) ** (-float(siteseqPssm[unpid][int(position) - 1][PSSMtitle[muta]])))
            DDpssm = round((mutapssm - wildpssm), 5)
            pssm = [round(1 / (1 + (2.7182) ** (-float(x))), 5) for x in siteseqPssm[unpid][int(position) - 1]]
            # 1 / (1 + (2.7182) ** (-int(data[j + 4])))
        pssm.append(DDpssm)
        # print(wildpssm,DDpssm)
        PSSMwm.append([wildpssm, DDpssm])
        PSSMwm21.append(pssm)

    print(PSSMwm21[0])
    print(len(PSSMwm21[0]))

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/wild_feature_seq/PSSMwm/pssm_wm_1676.json', 'w')
    json.dump(PSSMwm, pic)

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/wild_feature_seq/PSSMwm/pssm_21_1676.json', 'w')
    json.dump(PSSMwm21, pic)


class ParseHHm():

    def _loadFile(self, file):
        ''' Open the file with the path
            @return: The lines of the file
            @param file: The full path of the file, str
        '''
        try:
            with open(file) as fh:
                filelines = fh.readlines()  # read file and get all lines
            return filelines
        except IOError as err:
            print('File error: ' + err)

    def __parseonefilehhm(self, fin_data):

        hhm_begin_line = 0
        hhm_end_line = 0
        for i in range(len(fin_data)):
            if '#' in fin_data[i]:
                hhm_begin_line = i + 5
            elif '//' in fin_data[i]:
                hhm_end_line = i
        #print(hhm_end_line,hhm_begin_line)
        feature = np.zeros([int((hhm_end_line - hhm_begin_line) / 3), 30])
        axis_x = 0
        for i in range(hhm_begin_line, hhm_end_line, 3):
            line1 = fin_data[i].split()[2:-1]
            line2 = fin_data[i + 1].split()
            axis_y = 0
            for j in line1:
                if j == '*':
                    feature[axis_x][axis_y] = 9999 / 10000.0
                else:
                    feature[axis_x][axis_y] = float(j) / 10000.0
                axis_y += 1
            for j in line2:
                if j == '*':
                    feature[axis_x][axis_y] = 9999 / 10000.0
                else:
                    feature[axis_x][axis_y] = float(j) / 10000.0
                axis_y += 1
            axis_x += 1
        feature = (feature - np.min(feature)) / (np.max(feature) - np.min(feature))
        return feature

    def gethmmfeature(self, filepath):

        fileslist = os.listdir(filepath)
        allentryhmm = {}
        for onefile in fileslist:
            path = os.path.join(filepath,onefile)
            filelines = self._loadFile(path)

            pssmvalue = self.__parseonefilehhm(filelines)
            uniprotID = onefile.split('/')[-1].split('.')[0]

            uniprothmm = {uniprotID: pssmvalue}
            allentryhmm.update(uniprothmm)
        return allentryhmm

if __name__ == '__main__':
    PSSMtitle = {'A': 0, 'V': 19, 'L': 10, 'I': 9, 'M': 12, 'C': 4, 'F': 13, 'W': 17, 'Y': 18, 'H': 8, 'S': 15, 'T': 16,
                 'N': 2, 'Q': 5, 'K': 3, 'R': 1, 'D': 3, 'E': 6, 'G': 7, 'P': 14}

    # sites253 = pickle.load(open('/u1/home/jghhd/stability_change/Reimplement/253mutation_sites.pickle','rb'))
    hhm = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    hhmdict = {each:hhm.index(each) for each in hhm}

    # print(examples)
    # lines = open('/u1/home/jghhd/stability_change/Get_DSSPfile_Code/Sep_parse/motify_dataset.out').readlines()
    # lines = open('/u1/home/jghhd/stability_change/wild_feature/Win_seq/onesite.out').readlines()
    lines = open('/storage/htc/joshilab/jghhd/SC/stability_change/Data_pro/S236').readlines()
    fastapath = "/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s236_unp/SEQ/"
    filepath = "/storage/htc/joshilab/jghhd/SC/stability_change1/datasets_s236_unp/HHblitshhm/"
    pssm = ParseHHm()
    siteseqPssm = pssm.gethmmfeature(filepath)
    print(len(siteseqPssm))

    PSSMwm = []
    PSSMwm21 = []
    HHMori21,HHMoriwm = [],[]
    for i in range(1, len(lines)):
        line = lines[i].split()
        unpid = line[2].upper()
        wild = line[3][0]
        muta = line[3][-1]
        position = line[3][1:-1]
        print(hhmdict[wild])
        try:
            wildpssm = float(siteseqPssm[unpid][int(position) - 1][hhmdict[wild]])
            mutapssm = float(siteseqPssm[unpid][int(position) - 1][hhmdict[muta]])
            DDpssm = round((mutapssm - wildpssm), 5)
            pssm = [round(float(x), 5) for x in siteseqPssm[unpid][int(position) - 1]]
            oripssm = [float(x) for x in siteseqPssm[unpid][int(position) - 1]]
        except:
            unpid = line[2]
            print(unpid)
            wildpssm = float(siteseqPssm[unpid][int(position) - 1][hhmdict[wild]])
            mutapssm = float(siteseqPssm[unpid][int(position) - 1][hhmdict[muta]])
            DDpssm = round((mutapssm - wildpssm), 5)
            pssm = [round(float(x), 5) for x in siteseqPssm[unpid][int(position) - 1]]
            oripssm = [float(x) for x in siteseqPssm[unpid][int(position) - 1]]
            # 1 / (1 + (2.7182) ** (-int(data[j + 4])))
        pssm.append(DDpssm)
        # print(wildpssm,DDpssm)
        PSSMwm.append([wildpssm, DDpssm])
        PSSMwm21.append(pssm)
        HHMori21.append(oripssm+[mutapssm - wildpssm])
        HHMoriwm.append([wildpssm,mutapssm - wildpssm])

    print(PSSMwm21[0])
    print(len(PSSMwm21))
    print(len(PSSMwm21[0]))

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/IndependentData236/HHblits/hmm5_wm_236.json', 'w')
    json.dump(PSSMwm, pic)

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/IndependentData236/HHblits/hmm5_21_236.json', 'w')
    json.dump(PSSMwm21, pic)

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/IndependentData236/HHblits/hmmori_wm_236.json', 'w')
    json.dump(HHMoriwm, pic)

    pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/IndependentData236/HHblits/hmmori_21_236.json', 'w')
    json.dump(HHMori21, pic)

    # pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/wild_feature_seq/PSSMwm/pssm_wm_1676.json', 'w')
    # json.dump(PSSMwm, pic)
    #
    # pic = open('/storage/htc/joshilab/jghhd/SC/stability_change/wild_feature_seq/PSSMwm/pssm_21_1676.json', 'w')
    # json.dump(PSSMwm21, pic)

