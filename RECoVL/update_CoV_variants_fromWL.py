import sys
import csv

csv_file = open('variantswatchlist.csv')
read_csv = list(csv.reader(csv_file, delimiter=","))

file=[]
lin=[]
for c in read_csv[1:]:
    riga=c[4].split(' ')[1]
    if int(c[1])>=10:
        lin.append(riga)
        file.append(c)

lineage=list(set(lin))

dictionary_mutations={
    '69del': 'IHV68I',
    '70del': 'IHV68I' ,
    '242del': 'LALH242H',
    '243del': 'LALH242H',
    '244del': 'LALH242H',
    '144del': 'YY144Y'
}

f=open('watchlist_unique.txt', 'w')
for l in lineage:
    mutations = []
    for c in file:
        if l in c[4]:
            for m in c[0].split('_'):
                if 'X' not in m:
                    if m in dictionary_mutations:
                        mutations.append(dictionary_mutations.get(m))
                    else:
                        mutations.append(m)

    mutations_unique=list(set(mutations))
    print(mutations_unique)

    if len(mutations_unique)>0:
        f.write(l+': ')
        for c in mutations_unique[:-1]:
            f.write(c+', ')
        f.write(mutations_unique[-1])
        f.write('\n')




