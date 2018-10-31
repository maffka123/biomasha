import pandas as pd
from Bio import ExPASy, SwissProt
import mygene

def FixStr(s,flag):
    ss = ''
    if flag==1:
        s = s.partition('Full=')
        s1=s[2].split(';')
        ss+='<p>'+s1[0]+'</p>'
        s=s[2].partition('Short=')
        s1 = s[2].split(';')
        ss +='<p>'+s1[0]+'</p>'
    elif flag==2:
        s = s.replace("'", '')
        s = s.replace("[", '')
        s = s.replace("]", '')
        s = s.split(',')
        for k in s:
            if 'FUNC' in k or 'LOC' in k:
                ss += '<p>'+k+'</p>'
    return ss

def CombineData(mycsv):
    velos = pd.read_csv('Mitocheck_migration_scores_Ens89.txt', sep='\s+', header=0)
    # pp=pd.read_csv('genes_velo_description.csv',sep='\t', header=0)

    # ---Combine and get migration distance, sort
    final_data = pd.merge(mycsv, velos, left_on='symbol', right_on='Gene_symbol', how='left')
    final_data2 = final_data.groupby('EnsemblID')['Migration_(distance)'].mean().reset_index()
    final_data2 = final_data2.sort_values('Migration_(distance)', ascending=False).reset_index()
    final_data2 = pd.merge(final_data2, mycsv[['EnsemblID', 'symbol', 'file', 'other_phenotypes']], how='left')
    final_data2['UniProt_description'] = 'NaN'
    final_data2['UniProt_comments'] = 'NaN'
    final_data2.rename(columns={'index': 'uniprotID'}, inplace=True)
    return final_data2

def GetDescriptionFromWeb():
    # ---Get description of the genes
    try:
        final_data2=pd.read_csv('genes_binuc_description.csv',sep='\t')
        final_data2=final_data2[['EnsemblID','uniprotID','symbol','UniProt_description','UniProt_comments','file','other_phenotypes']]
        #final_data2.rename(columns={'Migration_(distance)':'Migration<p>(distance)</p>'}, inplace=True)
    except:
        # ---Read files
        mycsv = pd.read_csv('genes_binuc.csv', sep='\t', header=0)
        #final_data2=CombineData(mycsv)
        final_data2=mycsv[['geneID','EnsemblID','symbol','other_phenotypes','file']]
        Genenames = final_data2.EnsemblID
        mg = mygene.MyGeneInfo()
        for i, name in enumerate(Genenames):
            name2 = mg.query(name, scopes='ensembl.gene', fields='uniprot', species='human')
            try:
                name2 = name2['hits'][0]['uniprot']['Swiss-Prot']
                if type(name2) == list:
                    name2 = name2[0]
                print(str(i) + ' ' + name2)
                handle = ExPASy.get_sprot_raw(name2)
                record = SwissProt.read(handle)
                s = record.description
                s = FixStr(s, 1)
                final_data2.loc[i, 'UniProt_description'] = s
                s = str(record.comments)
                s=FixStr(s,2)
                final_data2.loc[i, 'UniProt_comments'] = s
                final_data2.loc[i, 'uniprotID'] = name2
            except:
                final_data2.loc[i, 'UniProt_description'] = 'no data found'
                final_data2.loc[i, 'UniProt_comments'] = 'no data found'
        final_data2.to_csv('genes_binuc_description.csv', sep='\t')
    return final_data2

#mycsv = pd.read_csv('genes_poly.csv', sep='\t', header=0)
final_data=GetDescriptionFromWeb()
#final_data2=pd.merge(final_data,mycsv[['EnsemblID','file','other_phenotypes']],left_on='EnsemblID', right_on='EnsemblID', how='left')
final_data=final_data[['EnsemblID','uniprotID','symbol','UniProt_description','UniProt_comments','file','other_phenotypes']]
#----add more data


pd.set_option('display.max_colwidth', -1)
f=open('genes_binuc.html','w')
out=final_data.to_html(na_rep = "", index = False,justify='center',escape=False)
out='<link rel="stylesheet" type="text/css" href="style.css">\n'+out
f.write(out)
f.close()
print(final_data.head(10))





#---------file readable for pandas
'''
for line in velos:
    v=line.split()
    writer = csv.writer(vel_new)
    writer.writerow(v)

vel_new.close()
'''