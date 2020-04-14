import csv
ortho_count = 0
ortho_output =[]
with open ("8.all.ort.group.tsv") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t")
    for line in tsvreader:
        #print (line[1])
        ortho_count += 1
        ortho_output.append(line)
#print (ortho_count)        
#print (ortho_output[1])

c831 = 0
o831 = []
cFac2 = 0
oFac2 = []
cFacM7 = 0
oFacM7 = []
cHom = 0
oHom = []
cInt5 = 0
oInt5 = []
cIntL = 0
oIntL = []
cIntX = 0
oIntX = []
cInu2 = 0
oInu2 = []
cInuD = 0
oInuD = []
cInuL = 0
oInuL = []

oList = [['R831', c831, o831], ['Rfac2', cFac2, oFac2], ['RfacM7', cFacM7, oFacM7], ['Rhom', cHom, oHom], ['Rint5', cInt5, oInt5], ['RintL', cIntL, oIntL], ['RintX', cIntX, oIntX], ['Rinu2', cInu2, oInu2], ['RinuD', cInuD, oInuD], ['RinuL', cInuL, oInuL]]
sList = ['R831','Rfac2','RfacM7','Rhom','Rint5','RintL','RintX','Rinu2','RinuD','RinuL']

R831pair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
F2pair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
FM7pair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Hompair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Int5pair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
IntLpair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
IntDpair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Inu2pair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
InuDpair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
InuLpair = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

pairwise_list = [R831pair, F2pair, FM7pair, Hompair, Int5pair, IntLpair, IntDpair, Inu2pair, InuDpair, InuLpair]

########-------Strain Homologs---------########
for strain in oList:
    org = strain[0]
    #print (org)
    for line in ortho_output:
        for item in line:
            if org in item:
                strain[1] += 1
                strain[2].append(line)

print ('Homologs for:')
print ('R. sp 831:', oList[0][1])
print ('R. faecis 2789:',oList[1][1])
print ('R. faecis M72:',oList[2][1])
print ('R. hominis A2-183:',oList[3][1])
print ('R. intestinalis 50/1:',oList[4][1])
print ('R. intestinalis LI-82:',oList[5][1])
print ('R. intestinalis XB6BA',oList[6][1])
print ('R. inulinivorans 2789:',oList[7][1])
print ('R. inulinivorans DSM16841:',oList[8][1])
print ('R. inulinivorans LI-83:',oList[9][1])

########-------Pairwise Homologs---------########
strain_count = -1                
for strain in oList:
    org = strain[0]
    strain_count += 1
    pair = pairwise_list[strain_count]
    #print (pair)
    paircount = 0
    for bact in sList:
        sub = bact
        for line in strain[2]:
            #print (strain[2])
            hit = next((s for s in line if sub in s),'-')
            if bact in hit:
                pair[paircount] += 1
                #print (pair)
        paircount += 1

print ('NR_ortholog Table:')
print ('R. sp 831:', pairwise_list[0])
print ('R. faecis 2789:',pairwise_list[1])
print ('R. faecis M72:',pairwise_list[2])
print ('R. hominis A2-183:',pairwise_list[3])
print ('R. intestinalis 50/1:',pairwise_list[4])
print ('R. intestinalis LI-82:',pairwise_list[5])
print ('R. intestinalis XB6BA',pairwise_list[6])
print ('R. inulinivorans 2789:',pairwise_list[7])
print ('R. inulinivorans DSM16841:',pairwise_list[8])
print ('R. inulinivorans LI-83:',pairwise_list[9])


########-------Faecis Homologs---------########
Rfac_orthologs = []
Rfac_orthologs_count = 0
for line in oList[1][2]:
    #print (line)
    sub ='RfacM7'
    hit = next((s for s in line if sub in s),'-')
    #print (hit)
    if 'RfacM7' in hit:
        Rfac_orthologs_count += 1
        Rfac_orthologs.append(line)

print ('# R. faecis homologs:',Rfac_orthologs_count)


########-------Intestinalis Homologs---------########
Rint5L_orthologs = []
Rint5L_orthologs_count = 0
for line in oList[4][2]:
    #print (line)
    sub ='RintL'
    hit = next((s for s in line if sub in s),'-')
    #print (hit)
    if 'RintL' in hit:
        Rint5L_orthologs_count += 1
        Rint5L_orthologs.append(line)

#print (Rint5L_orthologs_count)

Rint_orthologs = []
Rint_orthologs_count = 0
for line in Rint5L_orthologs:
    #print (line)
    sub ='RintX'
    hit = next((s for s in line if sub in s),'-')
    #print (hit)
    if 'RintX' in hit:
        Rint_orthologs_count += 1
        Rint_orthologs.append(line)

print ('# R. intestinalis homologs:',Rint_orthologs_count)


########-------Inulinivorans Homologs---------########
Rinu2D_orthologs = []
Rinu2D_orthologs_count = 0
for line in oList[7][2]:
    #print (line)
    sub ='RinuD'
    hit = next((s for s in line if sub in s),'-')
    #print (hit)
    if 'RinuD' in hit:
        Rinu2D_orthologs_count += 1
        Rinu2D_orthologs.append(line)

#print (Rinu2D_orthologs_count)

Rinu_orthologs = []
Rinu_orthologs_count = 0
for line in Rinu2D_orthologs:
    #print (line)
    sub ='RinuL'
    hit = next((s for s in line if sub in s),'-')
    #print (hit)
    if 'RinuL' in hit:
        Rinu_orthologs_count += 1
        Rinu_orthologs.append(line)

print ('# R. inulinivorans homologs:',Rinu_orthologs_count)



########-------Genus Homologs---------########
sp831 = o831
sp831_Fac2 = []
Fac2_FacM7 = []
FacM7_Hom = []
Hom_Int5 = []
Int5_IntL = []
IntL_IntX = []
IntX_Inu2 = []
Inu2_InuD = []
InuD_InuL = []
RoseOrthos = []

RoseList = [['R831', sp831], ['Rfac2', sp831_Fac2], ['RfacM7', Fac2_FacM7], ['Rhom', FacM7_Hom], ['Rint5', Hom_Int5], ['RintL', Int5_IntL], ['RintX', IntL_IntX], ['Rinu2', IntX_Inu2], ['RinuD', Inu2_InuD], ['RinuL', InuD_InuL]]
num_strain = len(RoseList)
#print (num_strain)

strain_count = 0                
for strain in RoseList:
    strain_count += 1
    sub = strain[0]
    #print (sub)
    #print (strain_count)
    for line in strain[1]:
        #print (strain[1])
        hit = next((s for s in line if sub in s),'-')
        #print (hit)
        if num_strain != strain_count:
            if sub in hit:
                RoseList[strain_count][1].append(line)
        else:
            if sub in hit:
                RoseOrthos.append(line)

print('# Rosebura Genus-wide Homologs:', len(RoseOrthos))

'''
########-------Interactions---------########
sp831 = o831
sp831_Fac2 = []#0Fac2
Fac2_FacM7 = []
FacM7_Hom = []
Hom_Int5 = [] #oHom
Int5_IntL = [] #oInt5
IntL_IntX = []
IntX_Inu2 = []
Inu2_InuD = [] #oInu2
InuD_InuL = []
RoseOrthos = []

RoseList = [['R831', sp831], ['Rfac2', sp831_Fac2], ['RfacM7', Fac2_FacM7], ['Rhom', FacM7_Hom], ['Rint5', Hom_Int5], ['RintL', Int5_IntL], ['RintX', IntL_IntX], ['Rinu2', IntX_Inu2], ['RinuD', Inu2_InuD], ['RinuL', InuD_InuL]]
num_strain = len(RoseList)
#print (num_strain)

strain_count = 0                
for strain in RoseList:
    strain_count += 1
    sub = strain[0]
    #print (sub)
    #print (strain_count)
    for line in strain[1]:
        #print (strain[1])
        hit = next((s for s in line if sub in s),'-')
        #print (hit)
        if num_strain != strain_count:
            if sub in hit:
                RoseList[strain_count][1].append(line)
        else:
            if sub in hit:
                RoseOrthos.append(line)

print('# Rosebura Genus-wide Homologs:', len(RoseOrthos))'''
