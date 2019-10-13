import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import copy
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap,BoundaryNorm #LinearSegmentedColormap
from itertools import chain

sns.set_context("paper")
sns.set_style("darkgrid")
colors = ['#e6f2ff', "#ff6666", '#88cc00','#ffa31a']
cmap = ListedColormap(colors)
boundaries = [0, 1, 2, 3, 4]
norm = BoundaryNorm(boundaries, cmap.N, clip=True)
fig, ax = plt.subplots(figsize=(20, 20))



def drawinit():
    x = np.arange(-50, 50)
    y = np.arange(-50, 50)
    xx, yy = np.meshgrid(x, y)
    s=xx.shape[0]
    ER = np.zeros((xx.shape[0],xx.shape[1]))
    ER[(np.square(xx)+np.square(yy))<40^2]=1
    u_sub=np.ones((xx.shape[0],xx.shape[1]))
    u_sub[(np.square(xx)+np.square(yy))<25^2]=0
    ER=ER*u_sub
    ER[45,45]=1
    MT = np.zeros((xx.shape[0], xx.shape[1]))
    MT[45, 45] = 1
    return ER, MT

def getdirchoices(dir_opt):
    if abs(dir_opt[0]) != abs(dir_opt[1]):
        if dir_opt[0] == 0:
            dirchoice = {1: [0, dir_opt[1]], 2: [1, dir_opt[1]], 3: [-1, dir_opt[1]]}
        else:
            dirchoice = {1: [dir_opt[0], 0], 2: [dir_opt[0], 1], 3: [dir_opt[0], -1]}
    else:
        dirchoice = {1: [dir_opt[0], dir_opt[1]], 2: [dir_opt[0], 0], 3: [0, dir_opt[1]]}
    return dirchoice

def GrowMTby1(MTs,mts,MT,bound=False):
    if bound:
        dir_opt = (np.array(MTs[mts][-1]) - MTs[mts][-2]).tolist()
        dirchoice={1:[dir_opt[0],dir_opt[1]],2:[dir_opt[0],dir_opt[1]],3:[dir_opt[0],dir_opt[1]]}
    else:
        if len(MTs[mts])==1: #after first step
            dir_opt=(np.array(MTs[mts][-1])-45).tolist()
            dirchoice=getdirchoices(dir_opt)
        else:
            dir_opt = (np.array(MTs[mts][-1]) - MTs[mts][-2]).tolist()
            dirchoice = getdirchoices(dir_opt)

    mtdir = dirchoice[np.random.randint(1, high=4, size=1).item()]
    toapp = np.array(MTs[mts][-1]) + mtdir
    if toapp[0]<1 or toapp[0]>99 or toapp[1]<1 or toapp[1]>99:
        toapp = np.array(MTs[mts][-1])
    MTs[mts].append(toapp.tolist())
    MT[toapp[0], toapp[1]] = 1
    return MT,MTs,toapp.tolist()

def appendMovie(MTs,p,ERtub,savelast=True):
    ER, MT = drawinit()

    b = np.array([item for sublist in np.array([*MTs.values()]) for item in sublist])
    c = np.array([item for sublist in np.array([*ERtub.values()]) for item in sublist])
    if b.size> 0:
        x = np.array(b[:, 0])
        y = np.array(b[:, 1])
        MT[x, y] = 1
    if c.size>0:
        x = np.array(c[:, 0])
        y = np.array(c[:, 1])
        ER[x, y] = 1
    if savelast:
        return MT, ER
    else:
        im = plt.imshow(MT+(ER*2), animated=True, cmap=cmap, norm=norm)
        plt.axis('off')
        ax = plt.gca()
        ax.grid(color='b', linestyle='-', linewidth=2)
        p.append([im])
        return p, MT,ER

def DrawFinal(MTs):
    plt.figure()
    for val in MTs:
        print('new MT')
        x = [45]
        x.extend(np.array(MTs[val])[:, 0])
        y = [45]
        y.extend(np.array(MTs[val])[:, 1])
        col = np.random.rand(1, 3)
        plt.scatter(x, y, alpha=0.3, c=col, s=3, marker="s")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.show()

#Create starting picture
ER, MT=drawinit()
p=[]

#Probabilities
ProbMTpickupER=0.5
ProbMTcollapse=0.03
ProbERfuse=0.5
ProbMTERcollapse=0.05

#Starting params
MTgif=np.empty((100, 100, 100))
MTgif[:,:,0]=copy.deepcopy(MT)
MTs = dict()
ERtub = dict()
MTER=dict()
ERcannot={}
mtpos=0
erpos=0
bound=dict()
MTstoloop=[]
gifid=1
forbiddStates=np.array([[-1,-1],[-1,0],[0,-1],[1,1],[0,1],[1,0],[1,-1],[-1,1]])

appendMovie(MTs,p,ERtub)

#Main time loop
for t in range(10000):
    #initial growth
    initdir=np.random.randint(-1,high=2,size=2)
    if initdir[0]!=0 and initdir[1]!=0:
        MTs[mtpos]=[[45+initdir[0],45+initdir[1]]]
        mtpos+=1

    MTstoloop=copy.deepcopy(MTs)
    #Loop over existing MTs and grow them
    for mts in MTstoloop:
        if mts in MTER: #if already bound
            if np.random.rand(1) > ProbMTERcollapse:  # if did not collapse
                MT, MTs, mtdir = GrowMTby1(MTs, mts, MT,bound=True)
                ERtub[MTER[mts]].append(mtdir)
                if ER[ERtub[MTER[mts]][-1][0],ERtub[MTER[mts]][-1][1]]==1 and np.random.rand(1)<ProbERfuse: #3-WJ form!
                    MTER.pop(mts)
            else:  # if collapsed
                ERcannot.pop(mts)
                MTs.pop(mts)
                ERtub.pop(MTER[mts])
                MTER.pop(mts)
        else:
            if np.random.rand(1) > ProbMTcollapse:  # if did not collapse
                MT, MTs, mtdir = GrowMTby1(MTs, mts, MT)
                if np.random.rand(1)<ProbMTpickupER and ER[MTs[mts][-1][0],MTs[mts][-1][1]]==1 and \
                        (np.array(mtdir)-1).tolist() not in list(chain.from_iterable(ERcannot.values())): #pull new ER tubule
                    ERtub[erpos]=[mtdir]
                    MTER[mts] = erpos
                    erpos+=1
                    ER[mtdir[0],mtdir[1]]=1
                    ERcannot[mts]=(np.array(mtdir)+forbiddStates).tolist()
            else:  # if collapsed
                MTs.pop(mts)


    p,MT,ER = appendMovie(MTs, p,ERtub,savelast=False)

#Plot everything in the end or save gif
if len(p)>0:
    ani = animation.ArtistAnimation(fig, p)
    ani.save('model.gif',fps=60, writer='PillowWriter')
else:
    im = plt.imshow(MT + (ER * 2), cmap=cmap, norm=norm)
    plt.axis('off')
    ax = plt.gca()
    ax.grid(color='b', linestyle='-', linewidth=2)
    plt.savefig('ProbFuse%0.2f_ProbCollapse%0.2f.png'%(ProbERfuse,ProbMTERcollapse))
