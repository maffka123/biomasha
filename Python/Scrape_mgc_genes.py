from bs4 import BeautifulSoup
import requests
import pandas as pd
from time import sleep


def get_seq(trs):
    '''
    Read sequence in table cell
    :param trs: line in table whith DNA or Protein
    :return: sequence
    '''
    td=trs.find_all('td')[1].get_text()
    return(td)

def scrape_content(soup,df):
    '''
    Major function for getting content of the webpage
    :param soup: soup of the webpage
    :param df: dataframe with all data
    :return: dataframe
    '''

    #Some links are empty, in that case page displys message with: Invalid ID given
    check_invalid=soup.get_text()
    if 'Invalid ID given' in check_invalid:
        print('...and did not exist')
        return(df)
    if 'eCIS' not in check_invalid: # we will than anyway navigate by CIS
        print('...and did not exist')
        return (df)

    #When everything is fine we keep reading the page
    t1 = soup.find_all('table')[0]
    t2= t1.find_all('table')[2]
    t3 = t2.find_all('table')[0]
    t4 = t3.find_all('table')[0]
    trs = t4.find_all('tr') # There are bunch of nested tables

    source=trs[0].find_all('td')[1].get_text().strip()
    source='_'.join(source.split(' ')[:2])
    accession=trs[1].find_all('td')[1].get_text().strip()
    eCIS=trs[2].find_all('td')[1].get_text().strip().replace('(Show this eCIS locus)','')
    eCIS=eCIS.strip()
    if eCIS=='' or eCIS==' ': # Sometimes it can be empty
        print('...and did not exist')
        return (df)

    if 'Product' in trs[6].find_all('td')[0].get_text(): # Sometimes it is not there
        Product=trs[6].find_all('td')[1].get_text().strip()
        dnpos=7
    else:
        Product ='None'
        dnpos = 6
    dna=get_seq(trs[dnpos])
    prot=get_seq(trs[dnpos+1])

    #Feel in the data frame
    df = df.append({
        'Source':source,
        'accession': accession,
        'eCIS': eCIS,
        'Product': Product,
        'DNA': dna,
        'Protein': prot}, ignore_index=True)

    return(df)

def main():
    '''
    Loop over all possible links, there are 10500 of them, maybe
    :return: saved csv
    '''

    base_link = 'http://www.mgc.ac.cn/cgi-bin/dbeCIS/showgene.cgi?ID=G'  # max 10500
    df = pd.DataFrame(columns=['Source','accession', 'eCIS', 'Product', 'DNA', 'Protein']) # will be dataframe where to save data

    #Start looping
    for i in range(1,10500):
        num='{foo:05d}'.format(foo=i) #append 0000 in front, links are 5 digit
        link=base_link+num # form working link
        page = requests.get(link) # get the full webpage
        soup = BeautifulSoup(page.content, "html.parser") #read it
        print('scraping page N {f1:d} was {f2:d}'.format(f1=i,f2=page.status_code))
        df=scrape_content(soup,df) # get content out
        df.to_csv('all_genes_new.csv', index=False) #save
        sleep(0.6) # keep some pause, not to break their server

if __name__ == "__main__":
    main()
