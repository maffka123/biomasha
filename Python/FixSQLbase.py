fd = open('mitocheck.sql', 'r')
sqlFile = fd.readlines()
#sqlFile = fd.read()
#sqlCommands = sqlFile.split(';')
fd.close()
tables=[]
newbase=''
for i,line in enumerate(sqlFile):
    #if 'LOC000000' in line:
    #    newbase+=line

    if 'KEY' in line and 'PRIMARY' not in line and 'FOREIGN' not in line and 'ENABLE' not in line and 'DISABLE' not in line:
        sqlFile[i]=''
        sqlFile[i-1] = sqlFile[i-1][:-2]+'\n'
    elif 'ENGINE' in line:
        sqlFile[i] = sqlFile[i].replace(line, ");\n")
    elif 'COMMENT' in line:
        l2=line.split('COMMENT')
        sqlFile[i] = l2[0]+',\n'
    elif 'unsigned' in line:
        sqlFile[i] = sqlFile[i].replace('unsigned', '')
    elif 'enum' in line:
        sqlFile[i] = sqlFile[i].replace("enum('0','1') NOT NULL DEFAULT '0'", 'NOT NULL')
        #sqlFile[i] = sqlFile[i].replace("enum('0','1') DEFAULT NULL", 'NOT NULL')
    elif "'&mu" in line:
        sqlFile[i] = sqlFile[i].replace("&mu;m<sup>2</sup>/s", 'mu2/s')

    #if 'CREATE' in line:
    #    tables.append(line)


file = open('mitocheck2.sql','w')
file.write(''.join(sqlFile))
file.close()
'''
file = open('tables.txt','w')
file.write(''.join(tables))
file.close()
'''
file = open('newbase.sql','w')
file.write(''.join(newbase))
file.close()
