import csv
import sys
import argparse

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def getAABase(AA):
    AABase=''
    for x in AA:
        if x.isdigit():
            AABase+=x
    return AABase

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_out_file', dest='intab', help='input_out_file')
    parser.add_argument('-o', '--output_out_file', dest='outtab', help='output_out_file')
    args = parser.parse_args()

    csv_file = open(args.intab)
    read_csv = list(csv.reader(csv_file, delimiter="\t"))
    out_file = open(args.outtab,'w')
    if len(read_csv)>1:
        out_file.write('\t'.join(read_csv[0])+'\n')
        index=1
        while index<len(read_csv[:-1]):
            aa=getAABase(read_csv[index][7])
            codone=[]
            riga=''
            non_trovato=0
            ultimo_dafare=True
    
            if 'FRAME_SHIFT' not in read_csv[index][5]:
                if '?' not in read_csv[index][7] and '*' not in read_csv[index][7]:
                    if '?' not in read_csv[index+1][7] and '*' not in read_csv[index+1][7]:
                        if aa==getAABase(read_csv[index+1][7]) and aa!='' and read_csv[index][0]==read_csv[index+1][0]:
                            codone.append(read_csv[index])
                            codone.append(read_csv[index+1])
                            index += 1
                            non_trovato+=1
                            if not index<len(read_csv[:-1]):
                                ultimo_dafare=False
                    if '?' not in read_csv[index + 1][7] and '*' not in read_csv[index + 1][7]:
                        if aa==getAABase(read_csv[index+1][7]) and aa!='' and read_csv[index][0]==read_csv[index+1][0]:
                            codone.append(read_csv[index+1])
                            index += 1
                            non_trovato+=1
                            if not index<len(read_csv[:-1]):
                                ultimo_dafare=False
            if non_trovato==0:
                out_file.write('\t'.join(read_csv[index])+'\n')
            if non_trovato==1:
                riga += codone[0][0] + '\t' + codone[0][1]+ '\t'
                cod=[]
                cod1=''
                cod2=''
                codon=''
                for x in range(len(codone)):
                    cod.append(codone[x][6])
                for x in range(0, 7):
                    if x<3:
                        if cod[0][x].isupper() and cod[1][x].islower():
                            codon+=cod[0][x]
                            cod1+=cod[0][x]
                        if cod[0][x].islower() and cod[1][x].isupper():
                            codon+=cod[1][x]
                            cod1+=cod[1][x]
                        if cod[0][x].islower() and cod[1][x].islower():
                            codon += cod[0][x]
                    if x>3:
                        if cod[0][x].isupper() and cod[1][x].islower():
                            codon+=cod[0][x]
                            cod2+=cod[0][x]
                        if cod[0][x].islower() and cod[1][x].isupper():
                            codon+=cod[1][x]
                            cod2+=cod[1][x]
                        if cod[0][x].islower() and cod[1][x].islower():
                            codon += cod[0][x]
                    if x==3:
                        codon+='/'
                mutations=''
                if gencode.get(codon[0:-4].upper())==gencode.get(codon[4:].upper()):
                    mutations+='SYNONYMOUS_CODING'
                else:
                    mutations+='NON_SYNONYMOUS_CODING'
                riga+=cod1+'\t'+cod2+'\t'+codone[0][0]+'\t'+mutations+'\t'+codon+'\t'+gencode.get(codon[0:-4].upper())+aa+gencode.get(codon[4:].upper())+'\n'
                out_file.write(riga)
            if non_trovato==2:
                newcodon=codone[0][6][4]+codone[1][6][5]+codone[2][6][6]
                riga+=codone[0][0]+'\t'+codone[0][1]+'\t'+read_csv[index][6][0:3].upper()+'\t'+newcodon+'\t'+codone[0][0]+'\t'+codone[0][5]+'\t'+read_csv[index][6][0:3].upper()+'/'+newcodon+'\t'+gencode.get(read_csv[index][6][0:3].upper())+aa+gencode.get(newcodon)+'\n'
                out_file.write(riga)
            index+=1
        if ultimo_dafare:
            out_file.write('\t'.join(read_csv[len(read_csv)-1])+'\n')
    else:
        out_file.write('error: no variants detected\n')
    out_file.close

if __name__ == "__main__":
    __main__()
