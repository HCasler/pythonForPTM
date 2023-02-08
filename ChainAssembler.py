#! usr/bin/env python
from ROOT import TChain
import os

class ChainAssembler:

    def __init__(self):
        self.jobDirPath = None # directory containing root files, or other directories with root files
        self.ntuplePath = None # the NTuple with the data you want
        self.outFilePaths = []
        self.chain = None

    def _isRootFile(self, filename):
        retval = False
        if filename[:4] == "nts." and filename[-5:] == ".root":
            retval = True
        return retval

    def _checkSubItemsRecursive(self, basepath, midpath):
        thisdir = basepath
        if midpath is not None:
            thisdir = basepath+"/"+midpath
        for f in os.listdir(thisdir):
            if os.path.isfile(thisdir+"/"+f) and self._isRootFile(f):
                outFilePath = thisdir+"/"+f 
                self.outFilePaths.append(outFilePath)
            elif os.path.isdir(thisdir+"/"+f):
                newMidpath = None
                if midpath is None:
                    newMidpath = f
                else:
                    newMidpath = midpath+"/"+f 
                self._checkSubItemsRecursive(basepath, newMidpath)


    def _collectOutFilePaths(self):
        # all objects inside job dir
        self._checkSubItemsRecursive(self.jobDirPath, None)


    def createChain(self):
        if self.jobDirPath is None:
            raise RuntimeError("ChainAssembler: must specify jobDirPath (file path to job output)")
        if self.ntuplePath is None:
            raise RuntimeError("ChainAssembler: must specify ntuplePath (ntuple name and directory in root files)")
        if len(self.outFilePaths) == 0: self._collectOutFilePaths()
        self.chain = TChain(self.ntuplePath)
        for filepath in self.outFilePaths:
            self.chain.Add(filepath)
        return self.chain
