import pandas as pd

dataTomap=pd.read_csv('siRNA_genes.csv', sep='\t')
dataToadd=pd.read_csv('Mitocheck_primary_screen-spot_scores.txt', sep='\s+', header=0,low_memory=False)
nuclei=['Large','Grape','Binuclear','Polylobed']
s='SIRNA,'+','.join(nuclei)
dataToadd=dataToadd[s.split(',')]
dataToadd=dataToadd.rename(index=str, columns={'SIRNA':'siRNA'})
dataToadd.drop(dataToadd[dataToadd[nuclei[0]].str.match('positive')].index, inplace=True)
#dataToadd.drop(dataToadd[dataToadd['siRNA'].str.match('empty',na=False)].index, inplace=True)
#dataToadd.drop(dataToadd[dataToadd['siRNA'].str.match('scramble',na=False)].index, inplace=True)
dataToadd.drop(dataToadd[dataToadd['siRNA'].isin(['scramble','empty','bCOP'])].index, inplace=True)
dataToadd.dropna(subset=['siRNA'], inplace=True)
dataToadd[nuclei]=dataToadd[nuclei].apply(pd.to_numeric)
dataToadd=dataToadd.groupby('siRNA')[nuclei].mean().reset_index()
dataToadd=pd.merge(dataToadd,dataTomap)
for i,name in enumerate(['large','grape','binuc','poly']):
    basetable=pd.read_csv('genes_{}_description.csv'.format(name), sep='\t')
    subdataToadd=dataToadd[['symbol',nuclei[i]]]
    basetable=pd.merge(basetable,subdataToadd)
    basetable=basetable.sort_values(nuclei[i], ascending=False).reset_index()
    basetable=basetable[['EnsemblID','uniprotID','symbol',nuclei[i],'UniProt_description','UniProt_comments','file','other_phenotypes']]
    basetable.to_csv('genes_{}_description2.csv'.format(name), sep='\t')
    pd.set_option('display.max_colwidth', -1)
    f = open('genes_{}.html'.format(name), 'w')
    out = basetable.to_html(na_rep="", index=False, justify='center', escape=False)
    out = '<link rel="stylesheet" type="text/css" href="style.css">\n' + out
    f.write(out)
    f.close()
