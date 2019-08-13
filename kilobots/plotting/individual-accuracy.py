import matplotlib.pyplot as plt
import matplotlib.cm
import csv
import numpy as np

data_file = "/Users/joefresna/DecisionsOnNetworks/kilobots/data/log_confsens_keep6.txt"

# Open the file
dataTab = []
with  open(data_file,"r") as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    for row in csvreader:
        dataTab.append(list(map(float, row)))

def plotAccOverL1(diff=None): 
    plt.grid()
    plt.plot( [d[0] for d in dataTab if diff == None or d[0]-d[1] in diff], [d[2]/d[3] for d in dataTab if diff == None or d[0]-d[1] in diff] )
    plt.show()

def plotAccOverDiff():
    robots=(len(dataTab[0])-4)/4
    allAcc=[]
    glbAcc=[]
    print([ (d[0],d[1]) for d in dataTab])
    for diff in np.arange(-3,4):
        #print("For Diff" + str(diff))
        acc=[]
        for r in np.arange(8, (int(robots)+2)*4, 4):
        #     print(r)
        #     print( np.mean( [d[r-1] for d in dataTab] ) )
            acc.append( np.mean( [d[r-1] for d in dataTab if d[0]-d[1] == diff ] ) )
            if (diff==-2): glbAcc.append( np.mean( [d[r-1] for d in dataTab ] ) )
        allAcc.append(acc)
    #print(len(np.arange(8, (int(robots)+2)*4, 4)))
    
    allAcc.append(glbAcc)
    plt.grid()
    plt.boxplot( allAcc )
    plt.show()

def plotAccOnLocation(diff=None):
    robots=(len(dataTab[0])-4)/4
    allAcc=[]
    glbAcc=[]
    robots=(len(dataTab[0])-4)/4
    for r in np.arange(8, (int(robots)+2)*4, 4):
        glbAcc.append( np.mean( [d[r-1] for d in dataTab if diff == None or d[0]-d[1] == diff ] ) )
    
    x = [dataTab[1][r] for r in np.arange(5, (int(robots)+1)*4, 4) ]
    y = [2000-dataTab[1][r] for r in np.arange(6, (int(robots)+1)*4, 4) ]
    colormap = matplotlib.cm.get_cmap('brg')
    for v,_ in enumerate(x):
        plt.text( x[v], y[v], round(glbAcc[v],2), fontsize=7, color=colormap(glbAcc[v]), weight="bold" )
    plt.xlim((0,2000))
    plt.ylim((0,2000))
    plt.axes().set_aspect('equal')
    plt.title("Diff:" + str(diff))
        
    #print(len(np.arange(8, (int(robots)+2)*4, 4)))
    
#     allAcc.append(glbAcc)
#     plt.boxplot( allAcc )
    plt.show()

if __name__ == '__main__':
    plotAccOverL1( [-1] )
    
    plotAccOverDiff()
#     #robots=(len(dataTab[0])-4)/4
#     print([ d[0] for d in dataTab ])
#     plotAccOnLocation(1)
#     plotAccOnLocation(-1)
#     plotAccOnLocation()
#print([d[1] for d in dataTab])