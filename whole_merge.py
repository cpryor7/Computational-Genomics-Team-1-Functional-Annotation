import sys

#unannotated gff loading
gff = sys.argv[7]
gff = open(gff)
gff = gff.readlines()
for i in range(0,len(gff)):
	gff[i] = gff[i].replace('\n','')
gff = [x for x in gff if (x.startswith('#')==False)]
gff_dic = {}
for i in gff:
	gff_dic[i.split('\t')[0]] = i

#name generating
name  = sys.argv[1]
a = open(name)
a = a.readlines()
id_list = []
id = {}
for i in a:
	id[i.split(' ')[0].split('>')[1]]=''
	id_list.append(i.split(' ')[0].split('>')[1])	

#egg format generating
egg = sys.argv[2]
b = open(egg)
data_egg = b.readlines()
axis = data_egg[4].replace('\n','')
axis = axis.split('\t')
b = data_egg[5:-3]
for i in range(0,len(b)):
	b[i] = b[i].replace('\n','')
for i in b:
	sample = i.split('\t')
	for j in range(1,len(sample)):
		if axis[j]=='Description':
			continue
		id[sample[0]] = id[sample[0]]+axis[j]+'='+sample[j]+';'

#rgi format generating
rgi = sys.argv[3]
c = open(rgi)
data_rgi = c.readlines()
axis = data_rgi[0].replace('\n','').split('\t')
c = data_rgi[1:]
for i in range(0,len(c)):
	c[i] = c[i].replace('\n','')
for i in c:
	sample = i.split('\t')
	query = sample[0].split('#')[0][:-1]
	for j in range(1,len(sample)):
		if axis[j]=='CARD_Protein_Sequence':
			continue
		if axis[j]=='Predicted_Protein':
			continue
		id[query] = id[query]+axis[j]+'='+sample[j]+';'

#vfdb format generating
vfdb = sys.argv[4]
d = open(vfdb)
data_vfdb = d.readlines()
d = data_vfdb
for i in range(0,len(d)):
	d[i] = d[i].replace('\n','')
for i in d:
	sample = i.split('\t')
	id[sample[0]]=id[sample[0]]+sample[-1]+';'

#signalP format generating
signalp = sys.argv[5]
e = open(signalp)
e = e.readlines()
e = e[1:]
for i in range(0,len(e)):
	e[i] = e[i].replace('\n','')
for i in e:
	sample = i.split(' ')
	id[sample[0]]=id[sample[0]]+'signalP_protein_type'+'='+sample[-1].split('\t')[2]+';'

#tmHmm format generating
tmhmm = sys.argv[6]
f = open(tmhmm)
f = f.readlines()
result = []
for i in range(0,len(f)):
	if f[i].startswith('#'):
		if 'Number of' in f[i]:
			result.append(f[i])
for i in range(0,len(result)):
	result[i] = result[i].replace('\n','')
for i in result:
	sample = i.split(' ')
	id[sample[1]]=id[sample[1]]+'#TMRs'+'='+sample[-1]

# output 
for i in range(0,len(gff)):
	gff[i] = gff[i] + id[id_list[i]]
	print(gff[i])