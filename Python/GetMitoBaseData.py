import requests
import pandas as pd
import sqlite3
from sqlite3 import OperationalError

# = requests.get('http://www.mitocheck.org/cgi-bin/mtc?action=get_data;ms_experiment=MSE00000205;data=summary_data;format=json')
#r = requests.get('http://www.mitocheck.org/cgi-bin/mtc?action=get_data;gene=ENSG00000178999;data=phenotypes;format=text')
def RedIntoDB(filename):
    fd = open(filename, 'r')
    sqlFile = fd.read()
    fd.close()
    conn = sqlite3.connect('mine.db')
    c = conn.cursor()
    sqlCommands = sqlFile.split(';')

    for i,command in enumerate(sqlCommands):
        try:
            c.execute(command)
            conn.commit()
            #print('I did '+command)
        except OperationalError as msg:
            if 'measurement' in command:
                print('I missed '+command)
                print(msg)
    return conn,c

#result = c.execute("select * from sys.tables;")
#print(result)


#conn,c=RedIntoDB ('mitocheck2.sql')
conn = sqlite3.connect('mine.db')
c = conn.cursor()

measurements=pd.read_sql("SELECT * FROM measurement;",conn)
measurementsU=measurements.name.unique()
#print(measurementsU)

# ----phenotype and images

siRNA=pd.read_sql("SELECT * FROM oligo_pair;",conn)
siRNA=siRNA[['external_ID','intended_target']]
siRNA=siRNA.rename(index=str, columns={'external_ID':'siRNA','intended_target':'symbol'})
siRNA.to_csv('siRNA_genes.csv', sep='\t')
phenotypes=pd.read_sql("SELECT * FROM phenotype;",conn)
image=pd.read_sql("SELECT * FROM image;",conn)
#print(image.filename[0])
#print(phenotypes[['phenotypeID','description']])
myphe='PHN00000009'
genes=pd.read_sql("SELECT * FROM gene_has_phenotype;",conn)
genes=genes.sort_values('phenotypeID')
genes_phe=genes[genes.phenotypeID==myphe]
list_genes=pd.read_sql("SELECT * FROM gene;",conn)
#print(list_genes.head())
list_genes=list_genes[['geneID','EnsemblID','symbol']]
#list_genes_phe=list_genes[list_genes.geneID in genes_phe.geneID]
#print(list_genes_phe)
genes_merged=pd.merge(genes_phe,list_genes, how='left')
genes_merged=genes_merged.drop_duplicates('geneID')
genes_merged['other_phenotypes']='NaN'
genes_merged['file']='NaN'
genesEns=genes_merged.EnsemblID
#print(genesEns)

c.close()
conn.close()
for i,gene in enumerate(genesEns):
    s='http://www.mitocheck.org/cgi-bin/mtc?action=get_data;gene={};data=images;format=json'.format(gene)
    r = requests.get(s)
    s=r.content.decode("utf-8").split('}')
    file=0
    k=0
    while file==0:
        lett=s[k]
        if 'small' in lett.lower():
            l=lett.replace('"','')
            l=l.partition('id:')
            l = l[2].partition(',')
            ff=l[0][1:]
            file=image.loc[image.image_setID==ff,'filename'].item()
            genes_merged.loc[i, 'file'] = '<video preload = "none" width="320" height="240" controls><source src="http://www.mitocheck.org/data'+file+'" type="video/mp4"></video>'
        k+=1
    try:
        geneid=genes_merged.loc[genes_merged.EnsemblID==gene,'geneID'].item()
    except:
        geneid = genes_merged.loc[genes_merged.EnsemblID == gene, 'geneID'].reset_index()
        geneid=geneid.geneID[0]
    phe=genes[genes.geneID==geneid]
    genesphe=pd.merge(phe,phenotypes,how='left').description
    s=''
    for phens in genesphe:
        s+='<p>'+phens+'</p>'
    genes_merged.loc[i, 'other_phenotypes'] =s
genes_merged.to_csv('genes_small.csv', sep='\t')



