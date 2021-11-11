import logging
from bamboo.analysisutils import loadPlotIt
import os.path
from bamboo.analysismodules import AnalysisModule, HistogramsModule


class CMSPhase2SimRTBModule(AnalysisModule):
    """ Base module for processing Phase2 flat trees """

    def __init__(self, args):
        super(CMSPhase2SimRTBModule, self).__init__(args)
        self._h_genwcount = {}

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import decorateCMSPhase2SimTree
        from bamboo.dataframebackend import DataframeBackend
        t = decorateCMSPhase2SimTree(tree, isMC=True)
        be, noSel = DataframeBackend.create(t)
        from bamboo.root import gbl
        self._h_genwcount[sample] = be.rootDF.Histo1D(
            gbl.ROOT.RDF.TH1DModel("h_count_genweight",
                                   "genweight sum", 1, 0., 1.),
            "_zero_for_stats",
            "genweight"
        )
        return t, noSel, be, tuple()

    def mergeCounters(self, outF, infileNames, sample=None):
        outF.cd()
        self._h_genwcount[sample].Write("h_count_genweight")

    def readCounters(self, resultsFile):
        return {"sumgenweight": resultsFile.Get("h_count_genweight").GetBinContent(1)}

class CMSPhase2SimHistoModule(CMSPhase2SimRTBModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """

    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Customised cutflow reports and plots """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList if (isinstance(
            ap, Plot) or isinstance(ap, DerivedPlot)) and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            printCutFlowReports(config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                        readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt,
                      eras=(eraMode, eras), verbose=self.args.verbose)

################################
  ## Actual analysis module ##
################################


class CMSPhase2Sim(CMSPhase2SimHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        plots = []

        # select photons
        photons = op.select(t.gamma, lambda ph: op.AND(op.abs(ph.eta) < 3, op.NOT(op.in_range(1.442, op.abs(ph.eta), 1.566)), ph.pt > 10))
        
        # select loose ISO photons
        ISOphotons = op.select(photons, lambda ph: ph.isopass & (1 << 0))
        
        # select loose ID photons
        IDphotons = op.select(ISOphotons, lambda ph: ph.idpass & (1 << 0))
                
        # Photon selections
        hasTwoISOPh = noSel.refine("hasTwoISOPh", cut=op.rng_len(ISOphotons) >= 2)
        
        hasTwoIDPh = noSel.refine("hasTwoIDPh", cut=op.rng_len(IDphotons) >= 2)

        # electrons
        electrons = op.select(t.elec, lambda el: op.AND(op.abs(el.eta) < 3, op.NOT(
            op.in_range(1.442, op.abs(el.eta), 1.566)), el.pt > 10.))
        
        # select loose ISO electrons
        ISOelectrons = op.select(electrons, lambda el: el.isopass & (1 << 0))

        # select loose ID electrons
        IDelectrons = op.select(ISOelectrons, lambda el: el.idpass & (1 << 0))
        
        # Electron selections
        hasOneEl = noSel.refine("hasOneElec", cut=[op.rng_len(ISOelectrons) >= 1])
        
        hasOneIDel = hasOneEl.refine("hasOneIDelec", cut=[op.rng_len(IDelectrons) >= 1])

        # muons
        muons = op.select(t.muon, lambda mu: op.AND(
            mu.pt > 10., op.abs(mu.eta) < 3))

        # select loose ISO muons
        ISOmuons = op.select(muons, lambda mu: mu.isopass & (1 << 0))

        # select loose ID muons
        IDmuons = op.select(
            ISOmuons, lambda mu: mu.idpass & (1 << 0))
        
        # muon selections
        hasOneMu = noSel.refine("hasOneMuon", cut=[op.rng_len(ISOmuons) >= 1])
        hasOneIDmu = hasOneMu.refine(
            "hasOneIDmuon", cut=[op.rng_len(IDmuons) >= 1])


       # plots

        plots.append(Plot.make1D("LeadingPhotonPTnonID", ISOphotons[0].pt, hasTwoISOPh, EqB(
            50, 0., 250.), title="Leading Photon pT"))
        plots.append(Plot.make1D("LeadingPhotonPtID_Eff", IDphotons[0].pt, hasTwoIDPh, EqB(
            50, 0., 250.), title="Leading Photon pT"))
        
        plots.append(Plot.make1D("LeadingElectronPTnonID", ISOelectrons[0].pt, hasOneEl, EqB(
            50, 0., 250.), title="Leading Electron pT"))
        plots.append(Plot.make1D("LeadingElectronPtID_Eff", IDelectrons[0].pt, hasOneIDel, EqB(
            50, 0., 250.), title="Leading Electron pT"))
        
        plots.append(Plot.make1D("LeadingMuonPTnonID", ISOmuons[0].pt, hasOneMu, EqB(
            50, 0., 250.), title="Leading Muon pT"))
        plots.append(Plot.make1D("LeadingMuonPtID_Eff", IDmuons[0].pt, hasOneIDmu, EqB(
            50, 0., 250.), title="Leading Muon pT"))
        
        return plots
