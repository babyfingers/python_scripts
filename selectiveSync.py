#!/usr/local/bin/python2.7

# selectiveSync.py:
# Use rsync -av to copy a directory tree, but allow
# for the exclusion of certain items/directories
# within that directory tree. The intended purpose is
# to make backing up user's home directory trees
# easier, esp. if they have a few enormous directories
# that do not need to be backed up.

itemsToSkip   = ["/data/home/alguire/inout_files/selectiveSync/testDir/donotcopy1/donotcopy2/donotcopy3","/data/home/alguire/inout_files/selectiveSync/testDir/donotcopy1/donotcopy2","/data/home/alguire/inout_files/selectiveSync/testDir/alsodonotcopy","/data/home/alguire/inout_files/selectiveSync/testDir/donotcopy1/okay_to_copy_yeah/literallypoison.txt"]
origCopyDir   = "/data/home/alguire/inout_files/selectiveSync/testDir"
backupCopyDir = "/data/home/alguire/inout_files/selectiveSync/backup"



import os

def dirStringSplit(dirString):
    """Takes a string describing the FULL PATH
    of a directory, splits it into an ordered
    list of directory names."""
    assert dirString[0]=="/"
    dirTree = dirString[1:].split("/")
    for listItem in dirTree:
        if (listItem==""): dirTree.remove("")
    return dirTree

def removeRedundantTrees(skipItemTreeList):
    """If we have one item tree that contains
    another, remove the redundant list."""
    treeRemovalIndexList = []
    NTree = len(skipItemTreeList)
    for iTreeIndex in range(NTree):
        iTree = skipItemTreeList[iTreeIndex]
        iTreeDepth = len(iTree)
        for jTreeIndex in range(NTree):
            redundancyFound = False
            jTree = skipItemTreeList[jTreeIndex]
            if (iTreeDepth>len(jTree) or iTreeIndex==jTreeIndex):
                continue
            print "Now comparing elements of the following trees:"
            print "iTree = "+str(iTree)
            print "jTree = "+str(jTree)
            for iDepth in range(iTreeDepth):
                if (iTree[iDepth]==jTree[iDepth]):
                    print "Element "+str(iDepth)+" agrees!"
                    redundancyFound = True
                else:
                    print "Element "+str(iDepth)+" does not agree!"
                    redundancyFound = False
                    break
            if redundancyFound:
                print "Awesome! We found a redundant tree, jTreeIndex = "+str(jTreeIndex)
                treeRemovalIndexList.append(jTreeIndex)
    treeRemovalIndexList = list(set(treeRemovalIndexList))
    treeRemovalIndexList = sorted(treeRemovalIndexList,key=lambda x: -x)
    for iTreeIndex in treeRemovalIndexList:
        print "Redundant tree removed: "+str(skipItemTreeList[iTreeIndex])
        skipItemTreeList.pop(iTreeIndex)
    print "Remaining skippable trees:"
    for s in skipItemTreeList:
        print s

def copyFilter(skipItemTreeList, treeIndex, treeDepth):
    myTree = skipItemTreeList[treeIndex]
    print "++++++++++++"
    print "Running copyFilter"
    print "++++++++++++"
    print "Current tree of interest = "+str(myTree)+" at depth = "+str(treeDepth)
    print "skipItemTreeList = "+str(skipItemTreeList)
    # Weed out any directory trees that have insufficient
    # depth, or have an item at this
    # depth, but which we don't expect to find in the
    # current directory. If we haven't reached the bottom
    # level of the tree, we only want to copy the
    # directory and not the contents.
    itemsToAvoid = []
    directoryCopyOnly = []
    for iTreeIndex in range(len(skipItemTreeList)):
        iTree = skipItemTreeList[iTreeIndex]
        print "Do we care about this tree? "+str(iTree)
        iTreeDepth = len(iTree)
        skipThisTree = False
        if (iTreeDepth-1<treeDepth):
            #myTreeIndexList.remove(iTreeIndex)
            print "Nope! Can't be found in this directory (we're too deep)."
            continue
        for depth in range(treeDepth):
            if myTree[depth] != iTree[depth]:
                #myTreeIndexList.remove(iTreeIndex)
                skipThisTree = True
                print "Nope! Can't be found in this directory (doesn't match)."
                print "iTree["+str(depth)+"] != myTree["+str(depth)+"]"
                print str(myTree[depth]) +"!="+str(iTree[depth])
                break
        if (skipThisTree): continue
#            print "Hey, I guess the tree's don't match, so we'll skip."
#            continue
        print "Yuup! It looks like we have to pay attention."
        if (treeDepth==iTreeDepth-1):
            itemsToAvoid.append(skipItemTreeList[iTreeIndex][treeDepth])
        else:
            directoryCopyOnly.append(skipItemTreeList[iTreeIndex][treeDepth])
    itemsToAvoid = tuple(set(itemsToAvoid))
    directoryCopyOnly = tuple(set(directoryCopyOnly))
    # itemsToAvoid should at least contain myTree[treeDepth]
    assert len(itemsToAvoid)+len(directoryCopyOnly)>0
    return itemsToAvoid, directoryCopyOnly

def refineItemList(itemList, itemsToAvoid, directoryCopyOnly):
    # Note that we expect that every item in a directory
    # should have a unique name, so the set operation
    # should not remove anything.
    itemList = list(set(itemList) -set(itemsToAvoid) -set(directoryCopyOnly))
    return itemList

def whereToStart(skipItemTreeList, treeIndex):
    """No need to duplicate our work; if a previous
    directory tree shares part of their path with the
    current one, skip the shared parts so we don't copy
    items in a directory more than once."""
    maxTreeDepth = 0
    for iTreeIndex in range(0,treeIndex):
        for iTreeDepth in range(len(skipItemTreeList[treeIndex])):
            if skipItemTreeList[treeIndex][iTreeDepth] != skipItemTreeList[iTreeIndex][iTreeDepth]:
                maxTreeDepth = max(maxTreeDepth,iTreeDepth)
                break
    return maxTreeDepth

def copyOver(copyItems,directoryCopyOnly,iPath,backupCopyDir):
    # Make sure the path exists in the backup
    # directory.
    copyDirectory = backupCopyDir+"/"+iPath
    print "copyDirectory = "+copyDirectory
    #if (not os.path.exists(copyDirectory)):
        #command = "rsync -av --exclude="+iPath+"/* "+iPath+" "+backupCopyDir+"/"+iPath 
        #command = "rsync -av -f\"+ "+iPath+"/\" -f\"- "+iPath+"/*\"" #TODO: copy directory iPath, but not contents
        #print command
        #os.system(command)
    assert os.path.isdir(copyDirectory)
    for item in copyItems:
        #command = "rsync -av "+iPath+"/"+item+" "+copyDirectory+"/"+iPath
        command = "rsync -av "+item+" "+copyDirectory
        print command
        os.system(command)
    for item in directoryCopyOnly:
        #command = "rsync -av --exclude="+iPath+"/"+item+"/* "+iPath+" "+backupCopyDir+"/"+iPath
        #command = "rsync -av --include='"+iPath+"/"+item+"' --exclude='"+iPath+"/"+item+"/*' "+iPath+"/"+item+" "+copyDirectory+"/"+iPath
        command = "rsync -avt --include='"+item+"' --exclude='"+item+"/*' "+item+" "+copyDirectory
        print command
        os.system(command)

def dirStringFromTree(itemTree,treeDepth):
    dirString = ""
    for iTreeDepth in range(treeDepth):
        dirString += itemTree[iTreeDepth]+"/"
    return dirString[:-1]

print "***************************************************"
print "**************** STARTING RUN *********************"
print "***************************************************"

assert os.path.exists(backupCopyDir)
assert os.path.isdir(backupCopyDir)
assert os.path.exists(origCopyDir)
assert os.path.isdir(origCopyDir)
origCopyDir = os.path.abspath(origCopyDir)
origCopyDirTree = dirStringSplit(origCopyDir)
origCopyTopDir = origCopyDirTree[-1]
backupCopyDir = os.path.abspath(backupCopyDir)

#Copy the top-level directory over first without any contents
#Note we must move into the containing directory to do this,
#because for some reason if you include the full path to the
#top-level directory it copies all the contents, too.
print "-------------------------------------------------"
print "Copying the top-level directory."
print "-------------------------------------------------"
os.chdir(origCopyDir+"/..")
command = "rsync -av --include='"+origCopyTopDir+"' --exclude='"+origCopyTopDir+"/*' "+origCopyTopDir+" "+backupCopyDir+"/"
print command
os.system(command)
backupCopyDir = os.path.abspath(backupCopyDir+"/"+origCopyTopDir)
os.chdir(origCopyTopDir)

# Prepare a standardized list of item trees
print "-------------------------------------------------"
print "Preparing a standardized list of skippable trees."
print "-------------------------------------------------"
skipItemTreeList = []
for skipItem in itemsToSkip:
    skipItem = os.path.abspath(skipItem)
    assert os.path.exists(skipItem)
    skipItemTree = dirStringSplit(skipItem)
    # Make sure the item to be skipped is in the directory
    # we are copying
    for iii in range(len(origCopyDirTree)):
        assert origCopyDirTree[iii] == skipItemTree[0]
        skipItemTree.pop(0)
    skipItemTreeList.append(skipItemTree)
removeRedundantTrees(skipItemTreeList)
#skipItemTreeList = tuple(map(None,*skipItemTreeList))


print "-------------------------------------------------"
print "Finally copying stuff."
print "-------------------------------------------------"
# Take the list of item trees and start rsyncing around
# it.
NTree = len(skipItemTreeList)
#for iTreeIndex in range(0):
for iTreeIndex in range(NTree):
    iTree = skipItemTreeList[iTreeIndex]
    print "==============================================="
    print "Working on iTree = "+str(iTree)
    print "==============================================="
    depthStart = whereToStart(skipItemTreeList, iTreeIndex)
    for iTreeDepth in range(depthStart,len(iTree)):
        iPath = dirStringFromTree(iTree, iTreeDepth)
        print "iPath = "+str(iPath)
        os.chdir(origCopyDir+"/"+iPath)
        allItems = os.listdir(".")
        print "allItems = "+str(allItems)
        avoidItems, directoryCopyOnly = copyFilter(skipItemTreeList, iTreeIndex, iTreeDepth)
        print "avoidItems = "+str(avoidItems)
        print "directoryCopyOnly = "+str(directoryCopyOnly)
        copyItems = refineItemList(allItems,avoidItems,directoryCopyOnly)
        print "copyItems = "+str(copyItems)
        copyOver(copyItems,directoryCopyOnly,iPath,backupCopyDir)

print "***************************************************"
print "**************** RUN COMPLETE *********************"
print "***************************************************"
