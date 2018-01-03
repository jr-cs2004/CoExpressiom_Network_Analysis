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
