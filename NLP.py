import re
import os
import math

word_re = re.compile("\W?(\w+)\W?")

def get_all_files(myDirectory):
    allFiles = []
    for root, dirs, files in os.walk(myDirectory):
        for iFile in files:
            allFiles.append(os.path.join(root,iFile))
    return allFiles

def load_file_tokens(filepath):
    myFile = open(filepath,"r")
    allTokens = []
    for line in myFile:
        allTokens.extend([myToken.lower() for myToken in word_re.findall(line)])
    myFile.close()
    return allTokens

def load_collection_tokens(myDirectory):
    filePathList = get_all_files(myDirectory)
    allTokens = []
    for myFilePath in fileList:
        allTokens.extend(load_file_tokens(myFilePath))
    return allTokens

def get_flatlist(itemlist):
    return [item for source, sourceitems in itemlist for item in sourceitems]

def get_types(flatlist):
    return tuple(sorted(list[set(flatlist)]))

def get_tf(flatlist):
    types = get_types(flatlist)
    histogram = [0]*len(types)
    for iii, iType in enumerate(types):
        histogram[iii] = flatlist.count(iType)
    normFac = 1.0/max(histogram)
    tfDict = {}
    for iii, iType in enumerate(types):
        tfDict[iType] = histogram[iii]*normFac
    return tfDict

def get_idf(biglist):
    #Remember, biglist is list of lists
    #of items, each from a given source.
    flatlist = get_flatlist(biglist)
    types = get_types(flatlist)
    idf = {}
    NDoc = len(biglist)
    for iii, iType in enumerate(types):
        count = 0
        for iList in biglist:
            if iType in iList:
                count += 1
        idf[iType] = math.log((1.0*NDoc)/count)
    return idf

def get_tfidf_top(TF, IDF, topN):
    TFIDF = [0]*len(TF)
    for iii, iType in enumerate(TF):
        TFIDF[iii] = (TF[iType]*IDF[iType],iType)
    TFIDF = sorted(TFIDF)
    return TFIDF[:topN]
    

if (__name__=="__main__"):
    # Task 1: get term frequency (TF) from the set
    # of all documents relating to a single company.
    # These are organized into directories in
    # Corpus root = '/home1/c/cis530/hw1/data/corpus'.
    corpusRoot = "/home1/c/cis530/hw1/data/corpus"
    companyNames = [myDir for myDir in os.listdir(corpusRoot) if os.path.isdir(os.path.join(corpusRoot,myDir))]
    companyPaths = [os.path.join(corpusRoot,myDir) for myDir in companyNames]
    companyTFs = {}
    for iName, iPath in zip(companyNames,companyPaths):
        companyTFs[iName] = get_tf(load_collection_tokens(iPath))
    # Task 2: get inverse document frequency (IDF)
    # from the complete list of documents.
    # all data = '/home1/c/cis530/hw1/data/all_data/'
    allData = "/home1/c/cis530/hw1/data/all_data/"
    documentList = get_all_files(allData)
    bigList = [load_collection_tokens(iDocument) for iDocument in documentList]
    IDFdict = get_idf(bigList)
    # Task 3: get top words for identifying companies
    # using TF*IDF.
    for iCompany, iTFdict in companyTFs:
        print iCompany +": ",get_tfidft_top(iTFdict,IDF,7)

