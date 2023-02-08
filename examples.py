#! usr/bin/env python

"""PTMPlotMaker is intended to make it quick and easy to write scripts for
 making lots of PTM plots and extracting simple stats from them.
 The classes in PTMReader (including PTMDetectorReader, VirtDetReader, etc)
 are used by PTMPlotMaker to make histograms from the data in root files.
 You can also call PTMReader classes on their own, and they include some
 functionality not used by PTMPlotMaker.

 Some of this functionality, like making the production target histograms,
 assumes you're using data generated by a branch of Offline with the 
 production target virtual detectors moved to be flush of the ends of the
 target, such as: github.com/HCasler/Offline/tree/forTargetScans """

from PTMPlotMaker import PTMPlotMaker
from ChainAssembler import ChainAssembler
from PTMReader import PTMDetectorReader
from ROOT import TCanvas

# PTMPlotMaker examples

print("Examples of using PTMPlotMaker:\n\n")

# making all plots at once
plotMaker = PTMPlotMaker()
plotMaker.dataPath = "sampleData"
plotMaker.jobName = "Example1_makeAllHistsAtOnce"
plotMaker.makeTargetHists = True
plotMaker.makePTMVirtualHists = True
plotMaker.makeScannerPlots = True
plotMaker.verbose = False
plotMaker.cleanupHists = False
plotMaker.makeAllPlots()

# since cleanupHists is False, we can access the plots already made
print("Plots still available:")
if len(plotMaker.heldHists.keys()) == 0:
    print("None")
for k in plotMaker.heldHists.keys():
    print(k)
# let's pick one and mess with it
potPlot = plotMaker.heldHists["Example1_makeAllHistsAtOnce POT"]
potPlot.Scale(18)
# This is a TH2I, so you can do all the stats-related stuff you normally could
xmean = potPlot.GetMean(1)
ymean = potPlot.GetMean(2)
potPlot.SetTitle("POT, scaled up 18x, mean = ({0:.1f}, {1:.1f})".format(xmean,  ymean))
c2 = TCanvas()
potPlot.Draw("lego")
c2.Print("Example2_changeSavedPlot.pdf", "pdf")
c2.Clear()

# this clears the ntuple tchains, stored plots, data path, and job name.
plotMaker.clearData()
# At this point, we could populate this with completely different data
plotMaker.dataPath = "sampleData2"
plotMaker.jobName = "Example3_makeOnlySomeHists"
# The following values are kept from before the data was cleared, so if we 
# want them to be different, we have to change them
plotMaker.makeTargetHists = False
plotMaker.makePTMVirtualHists = False
plotMaker.cleanupHists = True
# increase the e-dep-to-volts conversion factor, similar to increasing the bias
# voltage on a real detector
plotMaker.signalConversionConst = 0.135
plotMaker.makeAllPlots()

# since cleanupHists is True, the plots we made are gone after saving as pdfs
print("Plots still available:")
if len(plotMaker.heldHists.keys()) == 0:
    print("None")
for k in plotMaker.heldHists.keys():
    print(k)



# PTMReader is what actually pulls the data out of the ROOT files and makes
# the histograms for PTMPlotMaker.


print("\n\nExample of using PTMDetectorReader:\n\n")

# you have to create a TChain to read the data from.
# IMPORTANT: ChainAssembler assumes root file names are in the form of
# nts.otherStuffHere.root
assemb = ChainAssembler()
assemb.jobDirPath = "sampleData"
assemb.ntuplePath = "readPTM/ntPTM"
chain = assemb.createChain()

reader = PTMDetectorReader()

# Energy deposit histogram built into PTMDetectorReader, but note it requires 
# more fiddly input.
# vert1_volIds contains the lowest and highest volume IDs associated with the 
# first PWC's vertical signal plane.
eDepHist = reader.getIonizingEDepHist(chain, reader.vert1_volIds, name="Example4_usingPTMDetectorReader", maxVal=0.015)
eDepHist.SetTitle("Ionizing energy deposit in vertical plane of first PWC")
eDepHist.Draw('hist')
c2.Print(eDepHist.GetName()+".pdf", "pdf")
c2.Clear()