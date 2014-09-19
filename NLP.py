import re
import os

word_re = re.compile("\W?(\w+)\W?")

def get_all_files(myDirectory):
    allFiles = []
    for root, dirs, files in os.walk(myDirectory):
        for iFile in files:
            allFiles.append(root+"/"+iFile)
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

def 