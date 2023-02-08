#! usr/bin/env python
from ROOT import TH1F, TH1I, TH2I, TH2F
from array import array
from math import sqrt, pi

class PTMDetectorReader:
    """ This is for making histograms from the PTM gas volume sensitive 
    detectors. """

    def __init__(self):
        self.vert1_volIds = [0, 47]
        self.horiz1_volIds = [48, 95]
        self.vert2_volIds = [96, 143]
        self.horiz2_volIds = [144, 191]

    def _volIdToPosition(self, volId):
        numOnPlane = volId % 48
        centered = numOnPlane - 24
        scaledPosition = centered * 2
        return scaledPosition

    def _makeUniqueParticleId(self, eventId, trackId, pdgId, fnum):
        return "{0}_{1}_{2}_{3}".format(eventId, trackId, pdgId, fnum)

    def getIonizingProfiles(self, chain, namebase, pdgIDonly=[]):
        horizHits_1 = array('d', [])
        horizWeights_1 = array('d', [])
        vertHits_1 = array('d', [])
        vertWeights_1 = array('d', [])
        horizHits_2 = array('d', [])
        horizWeights_2 = array('d', [])
        vertHits_2 = array('d', [])
        vertWeights_2 = array('d', [])
        for entry in chain:
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                spatialPosition = self._volIdToPosition(entry.volId)
                iEdep = entry.iedep
                if entry.volId >= self.vert1_volIds[0] and entry.volId <= self.vert1_volIds[1]:
                    vertHits_1.append(spatialPosition)
                    vertWeights_1.append(iEdep)
                elif entry.volId >= self.vert2_volIds[0] and entry.volId <= self.vert2_volIds[1]:
                    vertHits_2.append(spatialPosition)
                    vertWeights_2.append(iEdep)
                elif entry.volId >= self.horiz1_volIds[0] and entry.volId <= self.horiz1_volIds[1]:
                    horizHits_1.append(spatialPosition)
                    horizWeights_1.append(iEdep)
                elif entry.volId >= self.horiz2_volIds[0] and entry.volId <= self.horiz2_volIds[1]:
                    horizHits_2.append(spatialPosition)
                    horizWeights_2.append(iEdep)
        horiz1 = TH1F(namebase+"horiz1", namebase+"horiz1", 48, -48, 48)
        horiz1.FillN(len(horizHits_1), horizHits_1, horizWeights_1, 1)
        horiz1.GetXaxis().SetTitle("horiz position (mm)")
        horiz1.GetYaxis().SetTitle("ionizing E dep (MeV)")
        horiz2 = TH1F(namebase+"horiz2", namebase+"horiz2", 48, -48, 48)
        horiz2.FillN(len(horizHits_2), horizHits_2, horizWeights_2, 1)
        horiz2.GetXaxis().SetTitle("horiz position (mm)")
        horiz2.GetYaxis().SetTitle("ionizing E dep (MeV)")
        vert1 = TH1F(namebase+"vert1", namebase+"vert1", 48, -48, 48)
        vert1.FillN(len(vertHits_1), vertHits_1, vertWeights_1, 1)
        vert1.GetXaxis().SetTitle("vert position (mm)")
        vert1.GetYaxis().SetTitle("ionizing E dep (MeV)")
        vert2 = TH1F(namebase+"vert2", namebase+"vert2", 48, -48, 48)
        vert2.FillN(len(vertHits_2), vertHits_2, vertWeights_2, 1)
        vert2.GetXaxis().SetTitle("vert position (mm)")
        vert2.GetYaxis().SetTitle("ionizing E dep (MeV)")

        outDict = {"horiz1": horiz1, "horiz2": horiz2, "vert1": vert1, "vert2": vert2}
        return outDict


    def getHitCountProfiles(self, chain, namebase, pdgIDonly=[]):
        horizHits_1 = array('d', [])
        vertHits_1 = array('d', [])
        horizHits_2 = array('d', [])
        vertHits_2 = array('d', [])
        for entry in chain:
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                spatialPosition = self._volIdToPosition(entry.volId)
                iEdep = entry.iedep
                if entry.volId >= self.vert1_volIds[0] and entry.volId <= self.vert1_volIds[1]:
                    vertHits_1.append(spatialPosition)
                elif entry.volId >= self.vert2_volIds[0] and entry.volId <= self.vert2_volIds[1]:
                    vertHits_2.append(spatialPosition)
                elif entry.volId >= self.horiz1_volIds[0] and entry.volId <= self.horiz1_volIds[1]:
                    horizHits_1.append(spatialPosition)
                elif entry.volId >= self.horiz2_volIds[0] and entry.volId <= self.horiz2_volIds[1]:
                    horizHits_2.append(spatialPosition)
        horiz1 = TH1F(namebase+"horiz1", namebase+"horiz1", 48, -48, 48)
        horiz1.FillN(len(horizHits_1), horizHits_1, array('d', [1 for i in horizHits_1]), 1)
        horiz1.GetXaxis().SetTitle("horiz position (mm)")
        horiz1.GetYaxis().SetTitle("ionizing E dep (MeV)")
        horiz2 = TH1F(namebase+"horiz2", namebase+"horiz2", 48, -48, 48)
        horiz2.FillN(len(horizHits_2), horizHits_2, array('d', [1 for i in horizHits_2]), 1)
        horiz2.GetXaxis().SetTitle("horiz position (mm)")
        horiz2.GetYaxis().SetTitle("ionizing E dep (MeV)")
        vert1 = TH1F(namebase+"vert1", namebase+"vert1", 48, -48, 48)
        vert1.FillN(len(vertHits_1), vertHits_1, array('d', [1 for i in vertHits_1]), 1)
        vert1.GetXaxis().SetTitle("vert position (mm)")
        vert1.GetYaxis().SetTitle("ionizing E dep (MeV)")
        vert2 = TH1F(namebase+"vert2", namebase+"vert2", 48, -48, 48)
        vert2.FillN(len(vertHits_2), vertHits_2, array('d', [1 for i in vertHits_2]), 1)
        vert2.GetXaxis().SetTitle("vert position (mm)")
        vert2.GetYaxis().SetTitle("ionizing E dep (MeV)")

        outDict = {"horiz1": horiz1, "horiz2": horiz2, "vert1": vert1, "vert2": vert2}
        return outDict

    def getIonizingEDepHist(self, chain, volIds, name=None, pdgIDonly=[], numBins=100, maxVal=None):
        if name is None:
            name = "Ionizing E Dep"
            if len(pdgIDonly) > 0:
                name += " for pdgIds: {0}".format(pdgIDonly)
        checkValues = False
        totalIEdeps = {}
        fnum = 1
        fpt = chain.GetFile()
        for entry in chain:
            if entry.volId < min(volIds) or entry.volId > max(volIds):
                continue
            thisFpt = chain.GetFile()
            if thisFpt != fpt:
                fnum += 1
                fpt = thisFpt
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                uniqueId = self._makeUniqueParticleId(entry.evt, entry.trk, entry.pdg, fnum)
                thisEdep = entry.iedep 
                if uniqueId in totalIEdeps:
                    totalIEdeps[uniqueId] += thisEdep
                else:
                    totalIEdeps[uniqueId] = thisEdep

        # if minVal is not None or maxVal is not None:
        #     checkValues = True
        # allEdeps = []
        # if checkValues:
        #     if minVal is None: minVal = 0
        #     if maxVal is None: maxVal = max(allEdeps)+0.5
        #     for u in totalIEdeps.values():
        #         if u >= minVal and u <= maxVal:
        #             allEdeps.append(u)
        # else: allEdeps = totalIEdeps.values()
        allEdeps = array('d', totalIEdeps.values())
        theMax = max(allEdeps) if maxVal is None else maxVal
        edepHist = TH1I(name, name, numBins, 0.0, maxVal)
        edepHist.FillN(len(allEdeps), allEdeps, array('d', [1 for i in allEdeps]), 1)
        edepHist.GetXaxis().SetTitle("ionizing E dep (MeV)")
        edepHist.GetYaxis().SetTitle("count")
        return edepHist

class VirtDetReader:
    """ This is for making histograms from virtual detectors. """

    def _getXEnds(self, xvals):
        xRange = max(xvals) - min(xvals)
        xmin = min(xvals) - 0.025*xRange
        xmax = max(xvals) + 0.025*xRange
        return xmin, xmax

    def _getYEnds(self, yvals):
        yRange = max(yvals) - min(yvals)
        ymin = min(yvals) - 0.025*yRange
        ymax = max(yvals) + 0.025*yRange
        return ymin, ymax

    def getPositionHist(self, chain, name=None, pdgIDonly=[], trackIDonly=[], binsPerSide=100, coordTransform=None):
        if name is None:
            name = "Virtual Detector Hit Positions"
            if len(pdgIDonly) > 0:
                name += " for pdgIds: {0}".format(pdgIDonly)
        xvals = array('d', [])
        yvals = array('d', [])
        for entry in chain:
            if (len(pdgIDonly) == 0 or entry.pdg in pdgIDonly) and (len(trackIDonly) == 0 or entry.trk in trackIDonly):
                x = entry.xl
                y = entry.yl 
                z = entry.zl
                if coordTransform is not None:
                    x, y, z = coordTransform(x, y, z)
                xvals.append(x)
                yvals.append(y)
        xmin, xmax = self._getXEnds(xvals)
        ymin, ymax = self._getYEnds(yvals)
        hitHist = TH2I(name, name, binsPerSide, xmin, xmax, binsPerSide, ymin, ymax)
        hitHist.FillN(len(xvals), xvals, yvals, array('d', [1 for i in xvals]), 1)
        hitHist.GetXaxis().SetTitle("x position (mm)")
        hitHist.GetYaxis().SetTitle("y position (mm)")
        return hitHist

    def getKEWieghtedPositionHist(self, chain, name=None, pdgIDonly=[], binsPerSide=100, coordTransform=None):
        if name is None:
            name = "Virtual Detector KE-Weighted Hit Positions"
            if len(pdgIDonly) > 0:
                name += " for pdgIds: {0}".format(pdgIDonly)
        xvals = array('d', [])
        yvals = array('d', [])
        kes = array('d', [])
        for entry in chain:
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                x = entry.xl
                y = entry.yl 
                z = entry.zl
                if coordTransform is not None:
                    x, y, z = coordTransform(x, y, z)
                xvals.append(x)
                yvals.append(y)
                kes.append(entry.ke)
        xmin, xmax = self._getXEnds(xvals)
        ymin, ymax = self._getYEnds(yvals)
        hitHist = TH2F(name, name, binsPerSide, xmin, xmax, binsPerSide, ymin, ymax)
        hitHist.FillN(len(xvals), xvals, yvals, kes, 1)
        hitHist.GetXaxis().SetTitle("x position (mm)")
        hitHist.GetYaxis().SetTitle("y position (mm)")
        return hitHist

    def getIncidentKEHist(self, chain, name=None, pdgIDonly=[], numBins=100):
        if name is None:
            name = "Virtual Detector Incident Kinetic Energy"
            if len(pdgIDonly) > 0:
                name += " for pdgIds: {0}".format(pdgIDonly)
        kes = array('d', [])
        for entry in chain:
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                kes.append(entry.ke)
        keHist = TH1I(name, name, numBins, 0.0, max(kes))
        keHist.FillN(len(kes), kes, array('d', [1 for i in kes]), 1)
        keHist.GetXaxis().SetTitle("incident KE (MeV)")
        keHist.GetYaxis().SetTitle("count")
        return keHist

    def getTotalParticleCount(self, chain, pdgIDonly=[]):
        totalCount = 0
        for entry in chain:
            if len(pdgIDonly) == 0 or entry.pdg in pdgIDonly:
                totalCount += 1
        return totalCount

    def getParticlesAccounting(self, chain):
        # returns a dict where the keys are PDGids and the values contain 
        # the total count of that particle and the mean and stdev of the 
        # kinetic energy
        # like:
        # {'2212': {'count':9000, 'keMean':7999.8, 'keStdev':103.7},
        #  '104' : {'count':33,   'keMean':910.76, 'keStdev':29.5}}
        accounting = {}
        for entry in chain:
            if str(entry.pdg) in accounting:
                thisEntry = accounting[str(entry.pdg)]
                thisEntry["keSum"] += entry.ke 
                thisEntry["keSqrSum"] += entry.ke**2
                thisEdep["count"] += 1
                accounting[str(entry.pdg)] = thisEntry
            else:
                thisEntry = {}
                thisEntry["keSum"] = entry.ke 
                thisEntry["keSqrSum"] = entry.ke**2
                thisEdep["count"] = 1
                accounting[str(entry.pdg)] = thisEntry
        for thisKey in accounting.keys:
            thisEntry = accounting[thisKey]
            keMean = thisEntry["keSum"] / thisEntry["count"]
            keVariance = (thisEntry["keSqrSum"] / thisEntry["count"]) - keMean**2
            keStdev = math.sqrt(keVariance)
            thisEntry["keMean"] = keMean
            thisEntry["keStdev"] = keStdev
            accounting[thisKey] = thisEntry
        return accounting



class PTMVirtDetReader(VirtDetReader):
    def _getXEnds(self, xvals):
        return -48, 48

    def _getYEnds(self, yvals):
        return -48, 48
