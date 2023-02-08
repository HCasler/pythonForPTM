#! usr/bin/env python
from ChainAssembler import ChainAssembler
from PTMReader import VirtDetReader, PTMVirtDetReader, PTMDetectorReader
from ROOT import TH1F, TCanvas, TVector3
from bisect import bisect_left
from math import pi, sqrt

class PTMPlotMaker:

    @staticmethod
    def targetFrontTransform(x, y, z):
        # This is for a virtual detector flush with the upstream face of
        # the target. This is not where this virtual detector is located in
        # Offline/main. See https://github.com/Mu2e/Offline/commit/1be2dbf354f00f868e44c572fc655a707e56ba69
        # (this commit is in github.com/HCasler/Offline/tree/forTargetScans)
        x = x - 3930.6141
        y = y + 0.00#
        z = z + 6177.7583
        angle = -14*pi/180
        vector = TVector3(x, y, z)
        vector.RotateY(angle)
        x = vector.x()
        y = vector.y()
        z = vector.z()
        return x, y, z

    @staticmethod
    def targetBackTransform(x, y, z):
        # This is for a virtual detector flush with the downstream face of
        # the target. This is not where this virtual detector is located in
        # Offline/main. See https://github.com/Mu2e/Offline/commit/1be2dbf354f00f868e44c572fc655a707e56ba69
        # (this commit is in github.com/HCasler/Offline/tree/forTargetScans)
        x = x - 3877.3898
        y = y + 0.00#
        z = z + 6151.2429
        angle = -14*pi/180
        vector = TVector3(x, y, z)
        vector.RotateY(angle)
        x = vector.x()
        y = vector.y()
        z = vector.z()
        return x, y, z

    @staticmethod
    def ptmVirtDetTransform(x, y, z):
        # The "local" coordinates of the PTM virtual detectors undergo a
        # 166 degree rotation, rather than a 14 degree rotation, because of
        # how the PTM gets constructed. Means the Mu2e coordinates are correct
        # but the local coordinates flip x and -x.
        x = -1*x 
        return x, y, z

    @staticmethod
    def getClosestWire(position, wirePositions):
        if position in wirePositions:
            return wirePositions.index(position)
        # otherwise
        insertionPoint = bisect_left(wirePositions, position)
        before = insertionPoint - 1
        after = insertionPoint
        beforeGap = abs(position - wirePositions[before])
        afterGap = abs(wirePositions[after] - position)
        if beforeGap < afterGap:
            return before
        else: 
            return after

    @staticmethod
    def volIdToPosition(volId):
        numOnPlane = volId % 48
        centered = numOnPlane - 24
        scaledPosition = centered * 2
        return scaledPosition

    def __init__(self):
        # user input data
        self.dataPath = None
        self.makeScannerPlots = False
        self.makeTargetHists = False
        self.makePTMVirtualHists = False
        self.jobName = None
        self.verbose = False
        self.cleanupHists = True

        # By default, scaled so 1e6 protons in a narrow peak (missing the
        # target) gets a signal peak height of 9.5 V.
        self.signalConversionConst = 0.01759
        # By default, the error on the TOTAL integrated signal in one detector
        # plane is 5%. This is distributed evenly over all the wires in the 
        # plane.
        self.totalSignalErr = 0.05

        # internal data
        self.targetFrontChain = None
        self.targetBackChain = None
        self.PTMSensitiveChain = None
        self.nearPwcVdChain = None
        self.farPwcVdChain = None

        # histograms we hold onto; by default nothing is in here
        self.heldHists = {}

    def verbosePrint(self, printout):
        if self.verbose:
            print(printout)


    def gatherChains(self):
        self.verbosePrint("Gathering chains...")
        if self.makeTargetHists:
            targFAssb = ChainAssembler()
            targFAssb.jobDirPath = self.dataPath
            targFAssb.ntuplePath = "readvdPTFront/ntvd"
            self.targetFrontChain = targFAssb.createChain()
            self.verbosePrint("Created chain from {0} using {1} data files".format(targFAssb.ntuplePath, len(targFAssb.outFilePaths)))
            targBAssb = ChainAssembler()
            targBAssb.jobDirPath = self.dataPath
            targBAssb.ntuplePath = "readvdPTBack/ntvd"
            self.targetBackChain = targBAssb.createChain()
            self.verbosePrint("Created chain from {0} using {1} data files".format(targBAssb.ntuplePath, len(targBAssb.outFilePaths)))
        if self.makePTMVirtualHists:
            nearAssb = ChainAssembler()
            nearAssb.jobDirPath = self.dataPath
            nearAssb.ntuplePath = "readvdNr/ntvd"
            self.nearPwcVdChain = nearAssb.createChain()
            self.verbosePrint("Created chain from {0} using {1} data files".format(nearAssb.ntuplePath, len(nearAssb.outFilePaths)))
            farAssb = ChainAssembler()
            farAssb.jobDirPath = self.dataPath
            farAssb.ntuplePath = "readvdFr/ntvd"
            self.farPwcVdChain = farAssb.createChain()
            self.verbosePrint("Created chain from {0} using {1} data files".format(farAssb.ntuplePath, len(farAssb.outFilePaths)))
        if self.makeScannerPlots:
            pwcAssb = ChainAssembler()
            pwcAssb.jobDirPath = self.dataPath
            pwcAssb.ntuplePath = "readPTM/ntPTM"
            self.PTMSensitiveChain = pwcAssb.createChain()
            self.verbosePrint("Created chain from {0} using {1} data files".format(pwcAssb.ntuplePath, len(pwcAssb.outFilePaths)))
        self.verbosePrint("...Chains gathered")

    def saveTargetHists(self, canvas, cleanupHists=True):
        self.verbosePrint("Making and saving proton target histograms")
        reader = VirtDetReader()
        # protons incident on the front face of the target
        savename = self.jobName+"_POT.pdf"
        incProtHist = reader.getPositionHist(self.targetFrontChain, name=self.jobName+" POT", trackIDonly=[1], coordTransform=PTMPlotMaker.targetFrontTransform)
        incProtHist.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        self.verbosePrint("POT hist done")
        # primary beam protons that make it through/past the target
        savename = self.jobName+"_prots_out_targ_back.pdf"
        outProtHist = reader.getPositionHist(self.targetBackChain, name=self.jobName+" p+ out target back", trackIDonly=[1], coordTransform=PTMPlotMaker.targetBackTransform)
        outProtHist.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        self.verbosePrint("Primary protons out back done")
        # ALL particles coming off the back of the target
        savename = self.jobName+"_all_out_targ_back.pdf"
        backProtHist = reader.getPositionHist(self.targetBackChain, name=self.jobName+" all out target back", coordTransform=PTMPlotMaker.targetBackTransform)
        backProtHist.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        self.verbosePrint("All particles out back done")

        if cleanupHists:
            del incProtHist 
            del outProtHist
            del backProtHist
        else:
            self.heldHists[incProtHist.GetName()] = incProtHist
            self.heldHists[outProtHist.GetName()] = outProtHist
            self.heldHists[backProtHist.GetName()] = backProtHist
        self.verbosePrint("Proton target histograms done")

    def savePTMVirtualHists(self, canvas, cleanupHists=True):
        self.verbosePrint("Making and saving PTM histograms")
        reader = VirtDetReader()
        # for each virtual detector, save both the just-primary-proton info and the all-particle info
        savename = self.jobName + "_near_PWC_prots.pdf"
        nearProts = reader.getPositionHist(self.nearPwcVdChain, name=self.jobName+" beam p+ on near PWC", trackIDonly=[1], binsPerSide=100, coordTransform=PTMPlotMaker.ptmVirtDetTransform)
        self.verbosePrint("Made hist with name {0}".format(nearProts.GetName()))
        nearProts.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        savename = self.jobName + "_near_PWC_all.pdf"
        nearAll = reader.getPositionHist(self.nearPwcVdChain, name=self.jobName+" all particles on near PWC", binsPerSide=100, coordTransform=PTMPlotMaker.ptmVirtDetTransform)
        nearAll.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        self.verbosePrint("Near PWC 2D histograms done")

        savename = self.jobName + "_far_PWC_prots.pdf"
        farProts = reader.getPositionHist(self.farPwcVdChain, name=self.jobName+" beam p+ on far PWC", trackIDonly=[1], binsPerSide=100)
        self.verbosePrint("Made hist with name {0}".format(farProts.GetName()))
        farProts.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        savename = self.jobName + "_far_PWC_all.pdf"
        farAll = reader.getPositionHist(self.farPwcVdChain, name=self.jobName+" all particles on far PWC", binsPerSide=100)
        farAll.Draw('colz')
        canvas.Print(savename, "pdf")
        canvas.Clear()
        self.verbosePrint("Far PWC 2D histograms done")

        if cleanupHists:
            del nearProts 
            del nearAll
            del farProts 
            del farAll
        else:
            self.heldHists[nearProts.GetName()] = nearProts
            self.heldHists[nearAll.GetName()] = nearAll
            self.heldHists[farProts.GetName()] = farProts
            self.heldHists[farAll.GetName()] = farAll
            self.verbosePrint("Keys in heldHists: {0}".format(self.heldHists.keys()))
        self.verbosePrint("PTM histograms done")

    


    def addBinErrs(self, hist):
        binSum = 0.0
        for i in range(1, hist.GetNbinsX()+1):
            binVal = hist.GetBinContent(i)
            binSum += binVal
        binErr = self.totalSignalErr * binSum / sqrt(48.0)
        self.verbosePrint("Bin err: {0:.4f} V".format(binErr))
        binErrSqrSum = 0.0
        for i in range(1, hist.GetNbinsX()+1):
            hist.SetBinError(i, binErr)
            binErrSqrSum += binErr**2
        binErrTotal = sqrt(binErrSqrSum)
        totalErrFrac = binErrTotal / binSum
        self.verbosePrint("Hist sum error frac: {0:.4f}".format(totalErrFrac))

    def saveScannerPlots(self, canvas, cleanupHists=True):
        self.verbosePrint("Making and saving scanner plots")
        reader = PTMDetectorReader()
        ionizingProfiles = reader.getIonizingProfiles(self.PTMSensitiveChain, self.jobName+"PTM_ionizing_")
        # Make the titles look nicer, and save the ionizing e dep data
        horizIon1 = ionizingProfiles["horiz1"]
        horizIon1.SetTitle("PTM PWC #1 horizontal: ionizing E dep")
        horizIon1.Draw('hist')
        canvas.Print(horizIon1.GetName()+".pdf", "pdf")
        canvas.Clear()
        horizIon2 = ionizingProfiles["horiz2"]
        horizIon2.SetTitle("PTM PWC #2 horizontal: ionizing E dep")
        horizIon2.Draw('hist')
        canvas.Print(horizIon2.GetName()+".pdf", "pdf")
        canvas.Clear()
        vertIon1 = ionizingProfiles["vert1"]
        vertIon1.SetTitle("PTM PWC #1 vertical: ionizing E dep")
        vertIon1.Draw('hist')
        canvas.Print(vertIon1.GetName()+".pdf", "pdf")
        canvas.Clear()
        vertIon2 = ionizingProfiles["vert2"]
        vertIon2.SetTitle("PTM PWC #2 vertical: ionizing E dep")
        vertIon2.Draw('hist')
        canvas.Print(vertIon2.GetName()+".pdf", "pdf")
        canvas.Clear()
        self.verbosePrint("Ionizing energy deposit profiles done")

        # Now to make the voltage signal plots
        horizSig1 = TH1F(horizIon1)
        horizSig1.Scale(self.signalConversionConst)
        self.addBinErrs(horizSig1)
        horizSig1.SetName(self.jobName+"_horizSignal_1")
        horizSig1.GetYaxis().SetTitle("scanner signal (V)")
        horizSig1.SetTitle("PTM PWC #1 horizontal: scanner signal")
        horizSig1.Draw("hist e1")
        canvas.Print(horizSig1.GetName()+".pdf", "pdf")
        canvas.Clear()

        horizSig2 = TH1F(horizIon2)
        horizSig2.Scale(self.signalConversionConst)
        self.addBinErrs(horizSig2)
        horizSig2.SetName(self.jobName+"_horizSignal_2")
        horizSig2.GetYaxis().SetTitle("scanner signal (V)")
        horizSig2.SetTitle("PTM PWC #2 horizontal: scanner signal")
        horizSig2.Draw("hist e1")
        canvas.Print(horizSig2.GetName()+".pdf", "pdf")
        canvas.Clear()

        vertSig1 = TH1F(vertIon1)
        vertSig1.Scale(self.signalConversionConst)
        self.addBinErrs(vertSig1)
        vertSig1.SetName(self.jobName+"_vertSignal_1")
        vertSig1.GetYaxis().SetTitle("scanner signal (V)")
        vertSig1.SetTitle("PTM PWC #1 vertical: scanner signal")
        vertSig1.Draw("hist e1")
        canvas.Print(vertSig1.GetName()+".pdf", "pdf")
        canvas.Clear()

        vertSig2 = TH1F(vertIon2)
        vertSig2.Scale(self.signalConversionConst)
        self.addBinErrs(vertSig2)
        vertSig2.SetName(self.jobName+"_vertSignal_2")
        vertSig2.GetYaxis().SetTitle("scanner signal (V)")
        vertSig2.SetTitle("PTM PWC #2 vertical: scanner signal")
        vertSig2.Draw("hist e1")
        canvas.Print(vertSig2.GetName()+".pdf", "pdf")
        canvas.Clear()

        # cleanup
        if cleanupHists:
            del horizIon1
            del horizIon2
            del vertIon1
            del vertIon2
            del horizSig1
            del horizSig2
            del vertSig1
            del vertSig2
        else:
            self.heldHists[horizIon1.GetName()] = horizIon1
            self.heldHists[horizIon2.GetName()] = horizIon2
            self.heldHists[vertIon1.GetName()] = vertIon1
            self.heldHists[vertIon2.GetName()] = vertIon2
            self.heldHists[horizSig1.GetName()] = horizSig1
            self.heldHists[horizSig2.GetName()] = horizSig2
            self.heldHists[vertSig1.GetName()] = vertSig1
            self.heldHists[vertSig2.GetName()] = vertSig2
        self.verbosePrint("All PWC scanner plots done")

    def makeAllPlots(self):
        self.verbosePrint("About to make all plots for job {0}".format(self.jobName))
        self.gatherChains()
        canvas = TCanvas()
        if self.makeTargetHists:
            self.saveTargetHists(canvas, cleanupHists=self.cleanupHists)
        if self.makePTMVirtualHists:
            self.savePTMVirtualHists(canvas, cleanupHists=self.cleanupHists)
        if self.makeScannerPlots:
            self.saveScannerPlots(canvas, cleanupHists=self.cleanupHists)
        self.verbosePrint("Finished all plots for job {0}".format(self.jobName))

    def redrawPlots(self, canvas, gpopt=None):
        if self.heldHists is None or len(self.heldHists) == 0:
            print("No held hists to re-save")
        else:
            for k in self.heldHists.keys():
                theHist = self.heldHists[k]
                if gpopt is None:
                    theHist.Draw()
                else:
                    theHist.Draw(gpopt)
                canvas.Print(theHist.GetName()+".pdf", "pdf")
                canvas.Clear()
            self.verbosePrint("Re-saved all held hists")

    def clearData(self):
        self.verbosePrint("Clearing internal data so this instance can be used again")
        self.dataPath = None
        self.jobName = None

        # internal data
        self.targetFrontChain = None
        self.targetBackChain = None
        self.PTMSensitiveChain = None
        self.nearPwcVdChain = None
        self.farPwcVdChain = None

        # histograms we hold onto; by default nothing is in here
        self.heldHists = {}



