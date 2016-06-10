import pylab as plt
from matplotlib_venn import venn2
import csv, math
import seaborn

def make_list(file):
    sets = []
    setM = []
    foldChange = []
    with open(file) as csvfile:   
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        reader.next()
        for line in reader:
            sets.append(line[0].split(':')[0])
            if math.isinf(abs(float(line[5]))): pass
            else: foldChange.append(abs(float(line[5])))
            if 'chrM' in line[0].split(':')[0]: setM.append(line[0].split(':')[0])
    return set(sets), set(setM), foldChange
    
def write_intersec(write, lis, txt):
    write.write(txt)
    for i in lis:
        write.write(i+'\n')

file1 = 'NegH_NGEM_ensemble.txt'
file2 = 'NegH6D_NGEM_ensemble.txt'


txt1 = file1.split('.txt')[0]
txt2 = file2.split('.txt')[0]   

subsetA, subMA, foldChA = make_list(file1)
subsetB, subMB, foldChB = make_list(file2)

intersec = subsetA.intersection(subsetB)
intersecM = subMA.intersection(subMB)

venn2([subsetA, subsetB], set_labels = (txt1, txt2))
plt.savefig('Compare_%s-%s.png' %(txt1, txt2), dpi=300)

plt.show()

venn2([subMA, subMB], set_labels = (txt1+'_chrM', txt2+'_chrM'))
plt.savefig('Compare_%s-%s_chrM.png' %(txt1, txt2), dpi=300)

plt.show()

fig, ax = plt.subplots(2)
#print foldChA
ax[0].hist(foldChA)
ax[1].hist(foldChB)
ax[0].set_title(txt1, x = 1, y = 0.8)
ax[1].set_title(txt2, x = 1, y = 0.8)
plt.savefig('Fold_changes_distr_%s-%s.png' %(txt1, txt2), dpi=300)

plt.show()

out_file = open('Compare_%s-%s.out' %(txt1, txt2), 'w')
write_intersec(out_file, intersec, "Same genes for %s and %s\n" %(txt1, txt2))
out_file.close()

if len(intersecM) > 0:
    out_file = open('Compare_%s-%s_chrM.out' %(txt1, txt2), 'w')
    write_intersec(out_file, intersecM, "Same mitochondrial genes for %s and %s\n" %(txt1, txt2))
    out_file.close()    
else: 
    out_file = open('Unique_%s-%s_chrM.out' %(txt1, txt2), 'w')
    write_intersec(out_file, subMA, "Unique mitochondrial genes for %s\n" %txt1)
    write_intersec(out_file, subMB, "Unique mitochondrial genes for %s\n" %txt2)
    out_file.close()
