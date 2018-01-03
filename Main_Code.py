import numpy as np
from scipy.stats import spearmanr

def read_raw_data(csv,metacols=None):
    """
    Reads expression data into a ndarray of expression data and tuple of gene identifiers.
    
    Parameters
    ----------
    csv : str
        Comma-separated csv file to read from. The first two columns are gene identifiers
        and prob_ids respectively. Next columns are expression data of each gene (rows)
        from different samples.
    metacols : int or None
        If csv contains meta data columns as defined in csv parameter, use this integer to
        separate them from expression data.
    
    Returns
    -------
    expression_data : ndarray or False
        The expression data. False when IO error happens.
    meta_data : ndarray
        The meta data matrix according to metacols provided.
    """
    try:
        raw_exp = np.genfromtxt(csv,delimiter=',',names=True,dtype=None)
        data_cols = list(raw_exp.dtype.names[metacols:])
        meta_cols = list(raw_exp.dtype.names[0:metacols])
        return raw_exp[data_cols].view(np.float64).reshape(raw_exp.shape[0],len(data_cols)),raw_exp[meta_cols]
    except IOError:
        print("Error in Raw Expression read")
        return False,None

def read_or_calc_corr_data(dataset):
    """
    Reads correlation data from files or calculate them by first loading raw expression data. In
    this case, it also saves files to disk for later use.
    
    Parameters
    ----------
    dataset : str
        The string that is used to produce filenames. "%s-rs" is spearman correlation file,
        "%s-rs-p" is the pvalue of the correlation matrix and "%s.csv" is the raw expression
        data file.
        
    Returns
    -------
    rs : ndarray
        Correlation matrix
    pvalue : ndarray
        Correlation p-value matrix.
    """
    rs_filename = "%s-rs" % dataset
    pvalue_filename = "%s-rs-p" % dataset
    dataset_filename = "%s.csv" % dataset

    try:
        rs = np.fromfile(rs_filename)
        pvalue = np.fromfile(pvalue_filename)
        size = np.sqrt(rs.shape[0])
        rs = rs.reshape((size,size))
        pvalue = pvalue.reshape((size,size))
    except IOError:
        print("Need to calculate spearman correlation matrix. This takes a few minutes...")
        exp,meta = read_raw_data(dataset_filename)
        rs,pvalue = spearmanr(exp,axis=1)
        rs.tofile(rs_filename)
        pvalue.tofile(pvalue_filename)
    
    return rs,pvalue

normal_rs,normal_pvalue = read_or_calc_corr_data("normal")
tumor_rs,tumor_pvalue = read_or_calc_corr_data("tumor")


%matplotlib inline
import matplotlib.pyplot as plt
fig = plt.figure()
fig.set_size_inches(14, 9)

f, axarr = plt.subplots(2)

bins = np.arange(-1,1,0.01)
normal_hist,normal_bins = np.histogram(normal_rs,bins)
tumor_hist,tumor_bins = np.histogram(tumor_rs,bins)

width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

axarr[0].bar(center, normal_hist, align='center', width=width)
axarr[1].bar(center, tumor_hist, align='center', width=width)


import networkx as nx

def network(rs,pvalue,p=0.2,corr=0.7):
    """
    Builds a network from Spearman correlation coefficient matrix.
    
    Parameters
    ----------
    rs : ndarray
        Spearman correlation coefficient matrix
    pvalue : ndarray
        Spearman correlation coefficient respective p-values
    p : float
        Threshold for maximum possible p-value to choose from
    corr : float
        Threshold for minimum possible correlation to form a network
    
    Returns
    -------
    G : Graph
        A simple undirected graph.
    Adj: ndarray
        Adjacancy matrix of the graph G. Adj[i][j] shows correlation between nodes i and j.
    """
    zeros = np.zeros(rs.shape)
    rs_sig = np.where(pvalue < p,rs,zeros)
    rs_adj = np.where(np.absolute(rs_sig) > corr,rs_sig,zeros)
    np.fill_diagonal(rs_adj,0)
    G = nx.from_numpy_matrix(rs_adj)
    return G,rs_adj

import itertools

pvalue_cuts = [0.05,0.1,0.15,0.2]
corr_cuts = [0.6,0.7,0.8,0.9]
for (p,c) in itertools.product(pvalue_cuts,corr_cuts):
    tumor_graph,tumor_adj = network(tumor_rs,tumor_pvalue,p=p,corr=c)
    tumor_gc = max(nx.connected_component_subgraphs(tumor_graph), key=len)
    print (p,c)
    print(tumor_gc.size())

normal_graph,normal_adj = network(normal_rs,normal_pvalue,p=0.05,corr=0.5)
#normal_gc = max(nx.connected_component_subgraphs(normal_graph), key=len)

tumor_graph,tumor_adj = network(tumor_rs,tumor_pvalue,p=0.05,corr=0.25)
#tumor_gc = max(nx.connected_component_subgraphs(tumor_graph), key=len)
nx.write_weighted_edgelist(normal_graph,"normal-graph.edges")
nx.write_weighted_edgelist(tumor_graph,"tumor-graph.edges")
#nx.write_gml(normal_graph,"normal-graph.gml")
#nx.write_gml(tumor_graph,"tumor-graph.gml")

%matplotlib inline
import matplotlib.pyplot as plt
fig = plt.figure()
fig.set_size_inches(14, 9)

normal_degree_seq=sorted(nx.degree(normal_graph).values(),reverse=True) 
tumor_degree_seq=sorted(nx.degree(tumor_graph).values(),reverse=True)
plt.semilogy(normal_degree_seq,color='blue')
plt.title("Degree distribution")
plt.ylabel("degree")
plt.xlabel("rank")
plt.semilogy(tumor_degree_seq,color='red')

print nx.info(normal_graph)
print nx.info(tumor_graph)

pos=nx.spring_layout(normal_graph,weight=None)



fig = plt.figure()
fig.set_size_inches(20,20)
fig.set_dpi(300)
nx.draw_networkx(normal_gc,pos=pos,with_labels=False,node_size=50,label="Normal Cell",width=2)
fig.axes[0].get_xaxis().set_visible(False)
fig.axes[0].get_yaxis().set_visible(False)

def tri(G):
    n1,n2,n3,n4 = 0,0,0,0
    n5 = 0
    for i in G.edges_iter():
        intersections = list(set(G[i[0]]).intersection(set(G[i[1]])))
        n5+=len(intersections)
        for j in intersections:
            if G[i[0]][i[1]]['weight'] > 0 and G[i[0]][j]['weight'] > 0 and G[i[1]][j]['weight'] > 0:
                n1+=1
            elif G[i[0]][i[1]]['weight'] < 0 and G[i[0]][j]['weight'] < 0 and G[i[1]][j]['weight'] < 0:
                n4+=1
            elif G[i[0]][i[1]]['weight'] * G[i[0]][j]['weight'] * G[i[1]][j]['weight'] < 0:
                n2+=1
            elif G[i[0]][i[1]]['weight'] * G[i[0]][j]['weight'] * G[i[1]][j]['weight'] > 0:
                n3+=1
    return n1/3,n2/3,n3/3,n4/3,n5/3

def tri_iter(G):
    for i in G.edges_iter():
        intersections = list(set(G[i[0]]).intersection(set(G[i[1]])))
        for j in intersections:
            #print(type(G[i[0]][j]))
            yield (i[0],i[1],j)
            #yield (G[i[0]][i[1]],G[i[0]][j],G[i[1]][j])
            
for t in tri_iter(normal_graph):
    print t
    break

n1,n2,n3,n4,n5 = tri(normal_graph)
print n1,n2,n3,n4,n5
print normal_graph.number_of_edges()


n1,n2,n3,n4,n5 = tri(tumor_graph)
print n1,n2,n3,n4,n5
print tumor_graph.number_of_edges()



tumor_graph,tumor_adj = network(tumor_rs,tumor_pvalue,p=0.05,corr=0.408)
print tumor_graph.number_of_edges()
