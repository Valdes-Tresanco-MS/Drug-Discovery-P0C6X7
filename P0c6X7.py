







import pandas as pd
import numpy as np
import urllib
import os 
import requests
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import seaborn as seabornInstance 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics
import re
import time





link1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/cgi2rcgi.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22protacxn%22:%22notnull%22},{%22cid%22:%22notnull%22},{%22repacxn%22:%22P0C6X7%22}]},%22order%22:[%22activity,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%" \
        "%22}"





for i in range(0,2):
    try:
        os.remove('downloaded1.csv')
        
        break
    except Exception as e:
        
        break
    else:
        break





for i in range(2):
    try:
        req = requests.get(link1)
        url_content = req.content
        csv_file = open('downloaded1.csv', 'wb')
        csv_file.write(url_content)
        csv_file.close()
        print("Completed the Request")
        break
    
     
      
    except Exception as a:
        print(str(a)+" is the error , Trying {} time".format(i))
        continue
    else:
        break
else:
    print("something Wrong , Try running Again [refer error code for more]")
            





data = pd.read_csv("downloaded1.csv",error_bad_lines=False)

data_df = pd.DataFrame(data)

new_data = data[['cid','acvalue']]


new_data = new_data.dropna()






except_val=0





cid_value = new_data['cid'].to_list()

PIC50_value = (-np.log10(new_data['acvalue']*10**-6)).to_list()





PIC50_value = pd.DataFrame(PIC50_value,columns = ["y"])





y_data = PIC50_value











link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
sub_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
sub_link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"

counter = 1





final_data_frame = pd.DataFrame()









x_data = pd.DataFrame()


cid_main_counter=0 








for x in cid_value:
    for iter_x in range(1000):
        try:    
            link = link_fixed
            link = link.replace("new_link_me",str(x))
            data1 = pd.read_csv(link)
            data1 = pd.DataFrame(data1)
            x_data = x_data.append(data1)
            data1 = 0
        except Exception as e:
            print("Exception Encountered as {} .! , Trying again !Iteration : {}".format(str(e),str(iter_x)))
            time.sleep(5)
            continue
        else:
            break
    
    else:
        print("Something Wrong with the Trials ! Restart The Algorithm!")





x_data_saved = x_data





x_data = x_data_saved











n_cid=[]






for i in cid_value:
    max_tries = 10
    
    for iter_ in range(max_tries):
        try:
            
            


            link = link_fixed
            link = link.replace("new_link_me",str(i))
            f_data = pd.read_csv(link) 
            
            sub_link = sub_link_fixed
            sub_link = sub_link.replace("replaceme",str(i))
            res = urllib.request.urlopen(sub_link)
            data_sub = res.read()
            data_sub = str(data_sub)
            data_sub = data_sub.replace("\\n",",") ; data_sub = data_sub.replace('b'," ")
            data_sub = data_sub.replace("'",""); data_sub = data_sub.replace(" ","")
            n_count_loop = 0 
            
            
            for j in data_sub.split(","):
                n_count_loop+=1
                if(len(j)>1):
                    if(n_count_loop<=15):
                        n_cid.append(j)
                    else:
                        break
                else:
                    pass
             
        except:
            if(except_val>=15):
                break
            else:
                print("Re-Trying")
                except_val+=1
                
                continue
        else:
            break
    






k_count = 0
t_count = 0
phase_count=0






f_data = pd.DataFrame()








for k in n_cid:
    k_count+=1
    t_count+=1
    phase_count+=1
    for sub_iter in range(10):
            try:
                if(t_count>=50):
                    
                    time.sleep(15)
                    t_count=0
                    
                
            
                link = link_fixed
                link = link.replace("new_link_me",str(k))
                
                f_data_df = pd.read_csv(link)
                f_data = f_data.append(f_data_df)
                
                break

            except Exception as e:
                print(str(e) + "Encountered , Please Wait :+ "+str(k_count))
                time.sleep(20)
                continue

    final_data_frame = final_data_frame.append(f_data)
    cid_main_counter+=1
    final_data_frame.to_csv("fdf1.csv")

    print("Data Fetching Continued : " + str(cid_main_counter))
     
    









final_data_frame.drop_duplicates(inplace=True)









x_data = pd.DataFrame(x_data)






y_data = pd.DataFrame(y_data)






train_file = pd.DataFrame()






x_data = x_data.astype("float64")
y_data = y_data.astype("float64")







new = y_data['y'].to_list()





train_file  = x_data
train_file['y'] = new





train_file
cid_reg_list=train_file['CID'].to_list()
train_file.drop('CID',axis=1,inplace =True)











for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')






x = list(train_file.columns)
x = x[:-1]

x_ = train_file[x].values
y_ = train_file['y'].values

















x_ = x_.astype("float64")
y_ = y_.astype("float64")




X_train, X_test, y_train, y_test = train_test_split(x_, y_, test_size=0.1, random_state=0)





print("Regression Starts")





regressor = LinearRegression()  
regressor.fit(X_train, y_train)






coef_dict = dict(zip(x, regressor.coef_))

y_pred = regressor.predict(X_test)





main_r2 = r2_score(y_test, y_pred, multioutput='uniform_average')
max_ = main_r2 



from itertools import combinations
comb_list = [[]]


def sub(arr,r):
    global comb_list
    for i in r:
        comb = list(combinations(arr,i))
        comb_list.append(comb)
    return comb_list











newone =0

newone = sub(x , [2,3,4,5,6])
del newone[0]


coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]






for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        regressor = LinearRegression()  
        regressor.fit(X_train, y_train)
        y_pred = regressor.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]






sec_index = max_r2.index(max(max_r2))



fir_index = index_r2[sec_index]






if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
maxi_r2 = max(max_r2)











reg_max_r2 = max(max_r2)









for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')











features  = list(features)





x_trained =train_file[features].values











y_trained = train_file['y'].values











X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)





regressor = LinearRegression()  
regressor.fit(X_train, y_train)





coef_dict = dict(zip(features, regressor.coef_))
coef_dict = pd.DataFrame(coef_dict , index = [0])





y_pred = regressor.predict(X_test)




saved_features = features

final_data_frame

final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]

final_data_frame = final_data_frameF2

pred_data_cid = final_data_frame["CID"].to_list()





features.append('CID')





final_data_frame1 = final_data_frame[features]





for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')





final_data_frame1.drop(columns='CID',inplace=True)


final_pred = regressor.predict(final_data_frame1)





final_pred = list(final_pred)





final_data_frame1['CID'] = pred_data_cid





final_data_frame1['y_'] = final_pred



saved_final_data_frame1 = final_data_frame1





sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)





final_larg= sorted_final_df.head(30)




final_larg = final_larg[['CID','y_']]





top_cid = final_larg["CID"].to_list()





for i in range(0,2):
    try:
        os.remove('downloaded1.csv')
        
        break
    except Exception as e:
        
        break
    else:
        break
for i in range(0,2):
    try:
        os.remove('fdf1.csv')
        
        break
    except Exception as e:
        
        break
    else:
        break





print("The Top 30 Drug Leads Which are identified with PubChem cid's are: ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))


import urllib.request
import os
l=[]
for itter in top_cid:
    l.append(itter)

def result2pdb(filename, ligname, recname):
    # check if exist docking result
    affinity = []
    all_result = open(f'{ligname}_{recname}result.pdb', 'w')
    models = 1
    result = open(filename, 'r')
    for line in result:
        line = line.strip('\n')
        if line[:5] == 'MODEL':
            all_result.write('MODEL' + ('%d' % models).rjust(9) + '\n')
            models += 1
        elif line[:5] == 'ATOM' or line[:6] == 'HETATM':
            all_result.write('HETATM' + line[7:67] + line[13:14].rjust(11) + '\n')
        elif line[:6] == 'ENDMDL':
            all_result.write(line + '\n')
        if re.search('REMARK VINA RESULT', line):
            affinity.append(float(line.split()[3].strip()))
            all_result.write('REMARK ENERGY   =   %s\n' % line.split()[3].strip())

    result.close()
    all_result.close()
    return affinity


try:
    if os.path.exists('ligands'):
        print('The ligand folder already exist. Renaming...')
        os.rename('ligands', f"ligands_old_{time.time()}")
    os.mkdir('ligands')
    os.chdir('ligands')
except IOError as e:
    raise e('It not possible create ligands folder...')

# info
lig_info = {}
# list of molecules
ms = []

for i in l:
    os.mkdir(f'{i}')
    os.mkdir(f'{i}/input')
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/SDF'.format(i)
    print(f'Downloading {i} compound file...')
    urllib.request.urlretrieve(url, f'{i}/input/PC_{i}.sdf')
    print(f'Downloading {i} compound file... Done.')
    print('Generating 3D mol files...')
    m = Chem.MolFromMolFile(f'{i}/input/PC_{i}.sdf')
    num_ha = m.GetNumAtoms()
    AllChem.Compute2DCoords(m)
    # add molecule to list for grid image
    ms.append(m)
    m2 = Chem.AddHs(m)
    # generate 3D coordinates
    AllChem.EmbedMolecule(m2, randomSeed=0xf00d)
    # write molecule in pdb format
    w = Chem.PDBWriter(f'{i}/input/PC_{i}.pdb')
    w.write(m2)
    print('Generating 3D mol files... Done')

    print('Preparing input...')

    p = subprocess.Popen([this_python, prepare_ligand4, '-l', f'{i}/input/PC_{i}.pdb', '-o', f'{i}/input/PC_{i}.pdbqt'])
    p.wait()
    print('Preparing input... Done')
    print('Molecular Docking...')
    p1 = subprocess.Popen([vina, '--config', '../configCLpro.txt', '--ligand', f'{i}/input/PC_{i}.pdbqt', '--out',
                           f'{i}/PC_{i}_CLout.pdbqt', '--log', f'{i}/PC_{i}_CLout.log', '--cpu', '8'])
    p1.wait()

    s = result2pdb(f'{i}/PC_{i}_CLout.pdbqt', i, 'CL')

    p2 = subprocess.Popen([vina, '--config', '../configPLpro.txt', '--ligand', f'{i}/input/PC_{i}.pdbqt', '--out',
                           f'{i}/PC_{i}_PLout.pdbqt', '--log', f'{i}/PC_{i}_PLout.log', '--cpu', '8'])
    p2.wait()

    s1 = result2pdb(f'{i}/PC_{i}_PLout.pdbqt', i, 'PL')
    print('Molecular Docking... Done.')

    # lig_info contain num_ha, list(affinity for CL), list(affinity for PL)
    lig_info[i] = [num_ha, s, s1]


Draw.DrawingOptions.atomLabelFontSize = 150
Draw.DrawingOptions.dotsPerAngstrom = 1500
Draw.DrawingOptions.bondLineWidth = 2
img = Draw.MolsToGridImage(ms, molsPerRow=5, subImgSize=(500, 450), legends=[x.GetProp("_Name") for x in ms],
                               legendFontSize=150)
img.save('2Dcompounds.png')
print('Drawing grind image... Done')

os.chdir(os.path.pardir)


print('=' * 139)
print(' Result for 3-Cysteine like Protease (3CLPro) [ PDB: 1p9u ]'.center(139))
print('-' * 139)
print('  Compound   |' + 'Poses/Ligand Efficiency'.center(125) + '|')
print('     ID      |' + '|'.join([('%s'.center(14)) % x for x in range(1, 10)]) + '|')
print('|'.join(['-' * 13 for x in range(0, 10)]) + '|')

for lig in lig_info:
    print(f'{lig}'.center(13) + '|' +
          '|'.join('{:.1f}/{:.2f}'.format(x, x / lig_info[lig][0]).center(13) for x in lig_info[lig][1]) + '|')
print('=' * 139)
print(' Result for Papain  like  protease  (PLpro)  [ PDB: 6w9c ]'.center(139))
print('-' * 139)
print('  Compound   |' + 'Poses/Ligand Efficiency'.center(125) + '|')
print('     ID      |' + '|'.join([('%s'.center(14)) % x for x in range(1, 10)]) + '|')
print('|'.join(['-' * 13 for x in range(0, 10)]) + '|')
for lig in lig_info:
    print(f'{lig}'.center(13) + '|' +
          '|'.join('{:.1f}/{:.2f}'.format(x, x / lig_info[lig][0]).center(13) for x in lig_info[lig][2]) + '|')
print('=' * 139)
