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

def finding_errors(positions):
    success='NOT_SUCCESS'
    errors={'FRAME_SHIFT','CHROMOSOME_LARGE_DELETION','CODON_CHANGE','CODON_INSERTION','CODON_CHANGE_PLUS_CODON_INSERTION','CODON_DELETION','CODON_CHANGE_PLUS_CODON_DELETION','CODON_INSERTION','CODON_CHANGE_PLUS_CODON_INSERTION','STOP_GAINED', '?', '*', 'STOP_GAINED'}
    degeneration={"R","D","M","N","S","K","W","H","B","V","Y","N"}
    for error in errors:
        for p in positions:
            if p.find(error)!=-1:
                success = 'SUCCESS'
                return success
    for deg in degeneration:
        for p in positions:
            if p==deg:
                success = 'SUCCESS'
                return success
    return success

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
    parser.add_argument('--minmax', dest='minmax', help='minority or majority variants')
    args = parser.parse_args()

    csv_file = open(args.intab)
    read_csv = list(csv.reader(csv_file, delimiter="\t"))
    out_file = open(args.outtab,'w')
    
    read_csv2=[]
    read_csv2.append(['Gene', 'Position', 'Reference', 'Alternative', 'Mutation type', 'Codon change', 'Amino Acid Effect'])
    for line in read_csv[1:]:
        if len(line)<=8 and line[3]!='N':
            del line[4]
            read_csv2.append(line)

    if len(read_csv2)>1:
        out_file.write('\t'.join(read_csv2[0])+'\n')
        index=1
        while index<len(read_csv2[:-1]):
            if finding_errors(read_csv2[index][3])=='SUCCESS':
                    out_file.write(read_csv2[index][0] + "\t" + read_csv2[index][1]  + "\t" + read_csv2[index][2]  + "\t" + read_csv2[index][3]+' - Degeneration\n')
                    index+=1
            else:
                aa=getAABase(read_csv2[index][6])
                codone=[]
                riga=''
                non_trovato=0
                if read_csv2[index][1]!=read_csv2[index+1][1]:
                    position = [read_csv2[index][3], read_csv2[index][4], read_csv2[index][6], read_csv2[index + 1][3], read_csv2[index + 1][4], read_csv2[index + 1][6]]
                    if finding_errors(position)=='NOT_SUCCESS':
                        if aa==getAABase(read_csv2[index+1][6]) and aa!='' and read_csv2[index][0]==read_csv2[index+1][0]:
                            codone.append(read_csv2[index])
                            codone.append(read_csv2[index+1])
                            index += 1
                            non_trovato+=1
                        if index+1 < len(read_csv2) and read_csv2[index][1]!=read_csv2[index+1][1]:
                            position = [read_csv2[index][3], read_csv2[index][4], read_csv2[index][6], read_csv2[index + 1][3], read_csv2[index + 1][4], read_csv2[index + 1][6]]
                            if finding_errors(position) == 'NOT_SUCCESS':
                                if aa==getAABase(read_csv2[index+1][6]) and aa!='' and read_csv2[index][0]==read_csv2[index+1][0]:
                                    codone.append(read_csv2[index+1])
                                    index += 1
                                    non_trovato+=1
                if non_trovato==0:
                    l=len(read_csv2[index][2])-1
                    if args.minmax == 'max':
                        out_file.write('\t'.join(read_csv2[index])+'\n')
                    elif read_csv2[index][4].find('FRAME_SHIFT')==-1:
                        out_file.write('\t'.join(read_csv2[index]) + '\n')
                    elif l%3==0 and l>=3:
                        out_file.write('\t'.join(read_csv2[index]) + '\n')       
                if non_trovato==1:
                    riga += codone[0][0] + '\t' + codone[0][1]+ '\t'
                    cod=[]
                    cod1=''
                    cod2=''
                    codon=''
                    for x in range(len(codone)):
                        cod.append(codone[x][5])
                    for x in range(0, 7):
                        if x<3:
                            if cod[0][x].isupper() and cod[1][x].islower():
                                codon+=cod[0][x]
                                cod1+=cod[0][x]
                            if cod[0][x].islower() and cod[1][x].isupper():
                                codon+=cod[1][x]
                                cod1+=cod[1][x]
                            if cod[0][x].islower() and cod[1][x].islower():
                                codon+=cod[0][x]
                        if x>3:
                            if cod[0][x].isupper() and cod[1][x].islower():
                                codon+=cod[0][x]
                                cod2+=cod[0][x]
                            if cod[0][x].islower() and cod[1][x].isupper():
                                codon+=cod[1][x]
                                cod2+=cod[1][x]
                            if cod[0][x].islower() and cod[1][x].islower():
                                codon+=cod[0][x]
                        if x==3:
                            codon+='/'
                    mutations=''
                    if gencode.get(codon[0:-4].upper())==gencode.get(codon[4:].upper()):
                        mutations+='SYNONYMOUS_CODING'
                        riga += cod1 + '\t' + cod2 + '\t' + mutations + '\t' + codon + '\t' + gencode.get(codon[0:-4].upper()) + str(aa) + '\n'
                    else:
                        mutations+='NON_SYNONYMOUS_CODING'
                        riga += cod1 + '\t' + cod2 + '\t' + mutations + '\t' + codon + '\t' + gencode.get(codon[0:-4].upper()) + str(aa) + gencode.get(codon[4:].upper()) + '\n'
                    out_file.write(riga)
                if non_trovato==2:
                    mutations=''
                    newcodon=codone[0][5][4]+codone[1][5][5]+codone[2][5][6]
                    if gencode.get(read_csv2[index][5][0:3].upper())==gencode.get(newcodon):
                        mutations += 'SYNONYMOUS_CODING'
                        riga += cod1 + '\t' + cod2 + '\t' + mutations + '\t' + codon + '\t' + gencode.get(codon[0:-4].upper()) + str(aa) + '\n'
                    else:
                        mutations += 'NON_SYNONYMOUS_CODING'
                        riga+=codone[0][0]+'\t'+codone[0][1]+'\t'+read_csv2[index][5][0:3].upper()+'\t'+newcodon+'\t'+mutations+'\t'+read_csv2[index][5][0:3].upper()+'/'+newcodon+'\t'+gencode.get(read_csv2[index][5][0:3].upper())+aa+gencode.get(newcodon)+'\n'
                    out_file.write(riga)
                index+=1
        if read_csv2[len(read_csv2)-1][4].find('FRAME_SHIFT')==-1:
            out_file.write('\t'.join(read_csv2[len(read_csv2)-1])+'\n')
    else:
        out_file.write('error: no variants detected\n')
    out_file.close

if __name__ == "__main__":
    __main__()
