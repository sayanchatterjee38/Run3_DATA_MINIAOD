import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023

process = cms.Process("SayanCMW")

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('MergingProducer.generalAndHiPixelTracks.MergingPixAndGenProducer_cfi')

# __________________ General _________________

# Configure the logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
)

# Configure the number of maximum event the analyser run on in interactive mode
# -1 == ALL
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(1000) 
)


runOnDATA = True ###IMPORTANT: please, to run in MC set to "False"
if runOnDATA :
   rootFiles = open("/eos/cms/store/group/phys_heavyions/sayan/TrackingTools_run3_datataking/pbpb_miniaod_datasets/run374810_ls0001_ls1752_streamPhysicsHIPhysicsRawPrime0_miniaod.txt", "r").readlines()

else :
   rootFiles = '/store/himc/HINPbPbSpring21MiniAOD/MinBias_Hydjet_Drum5F_2018_5p02TeV/MINIAODSIM/NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/280000/091945fa-aa74-4dde-9891-b620f03a6749.root'
process.source = cms.Source("PoolSource",
                fileNames = cms.untracked.vstring(rootFiles),
                skipBadFiles=cms.untracked.bool(True),
	        duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

##json
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/eos/cms/store/group/phys_heavyions/sayan/HIN_run3_pseudo_JSON/HIPhysicsRawPrime/Golden_Online_live.json').getVLuminosityBlockRange()

### centrality ###
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

# Set the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnDATA :
   process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Express_v4', '') #centrality table information from this GT for DATA
   process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
   process.GlobalTag.toGet.extend([
      cms.PSet(record = cms.string("HeavyIonRcd"),
               tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
               connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
               label = cms.untracked.string("HFtowers")
      ),
   ])
   
else :
   process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_HI_v2', '') #centrality table information from this GT for MC


# __________________ Event selection _________________
# Add PbPb collision event selection

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')


# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.phfCoincFilter2Th4 *
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)


# Define the output
process.TFileService = cms.Service("TFileService",fileName = cms.string('pbpb_miniAOD_nonsym_C2_cent080_run3_genAndHiPixTrk.root'))

###trigger selection for data
#process.hltMB.andOr = cms.bool(True)  # True = OR, False = AND between the HLT paths
#process.hltMB.throw = cms.bool(False) # throw exception on unknown path names

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltfilter = hltHighLevel.clone(
    HLTPaths = [
        "HLT_HIMinimumBiasHF1ANDZDC1nOR_v1",
    ]
)


# Load you analyzer with initial configuration
process.load("Analyzers.SayanCMW.SayanCMW_cfi")
process.defaultCPDC.vertex = 'offlineSlimmedPrimaryVertices'
process.defaultCPDC.nonsym = False  ##Set True if you want non-symmetric filling C2 correlation
process.defaultCPDC.genTrk = False  ##It will be always False when running genAndHiPixel tracks
process.defaultCPDC.genAndHiPixTrk = True  ##It will be always True when running genAndHiPixel tracks

process.defaultCPDC.ptmin_trigg = 0.3
process.defaultCPDC.ptmax_trigg = 4.0
process.defaultCPDC.ptmin_ass = 0.3
process.defaultCPDC.ptmax_ass = 4.0
## Important : trigger_ptbin should be equal to associate_ptbin
process.defaultCPDC.trigger_ptbin = cms.untracked.vdouble(0.3, 1.0, 2.0, 3.0, 4.0)   ##trigger_ptbin[0] = ptmin_trigg   ##trigger_ptbin[trigger_ptbin.size() -1] = ptmax_trigg   
process.defaultCPDC.associate_ptbin = cms.untracked.vdouble(0.3, 1.0, 2.0, 3.0, 4.0)  ##associate_ptbin[0] = ptmin_ass   ##associate_ptbin[associate_ptbin.size() -1] = ptmax_ass



process.load('MergingProducer.generalAndHiPixelTracks.MergingPixAndGenProducer_cfi')
process.generalAndHiPixelTracks.vertexSrc = 'offlineSlimmedPrimaryVertices'


process.load("Analyzers.SayanCMW.SayanCMW_cff")

process.defaultAnalysis_05     = process.CPDC05.clone()
process.defaultAnalysis_510    = process.CPDC510.clone()

process.defaultAnalysis_1020   = process.CPDC1020.clone()
process.defaultAnalysis_2030   = process.CPDC2030.clone()
process.defaultAnalysis_3040   = process.CPDC3040.clone()
process.defaultAnalysis_4050   = process.CPDC4050.clone()
process.defaultAnalysis_5060   = process.CPDC5060.clone()
process.defaultAnalysis_6070   = process.CPDC6070.clone()
process.defaultAnalysis_7080   = process.CPDC7080.clone()

process.defaultAnalysis_1015   = process.CPDC1015.clone()
process.defaultAnalysis_1520   = process.CPDC1520.clone()
process.defaultAnalysis_2025   = process.CPDC2025.clone()
process.defaultAnalysis_2530   = process.CPDC2530.clone()

process.defaultAnalysis_3035   = process.CPDC3035.clone()
process.defaultAnalysis_3540   = process.CPDC3540.clone()
process.defaultAnalysis_4045   = process.CPDC4045.clone()
process.defaultAnalysis_4550   = process.CPDC4550.clone()

process.defaultAnalysis_5055   = process.CPDC5055.clone()
process.defaultAnalysis_5560   = process.CPDC5560.clone()
process.defaultAnalysis_6065   = process.CPDC6065.clone()
process.defaultAnalysis_6570   = process.CPDC6570.clone()

process.defaultAnalysis_7075   = process.CPDC7075.clone()
process.defaultAnalysis_7580   = process.CPDC7580.clone()
process.defaultAnalysis_8085   = process.CPDC8085.clone()
process.defaultAnalysis_8590   = process.CPDC8590.clone()



process.p = cms.Path(process.eventFilter *
                     process.hltfilter *
                     process.centralityBin *                 # Compute centrality
                     #process.hiEvtAnalyzer *
                     process.generalAndHiPixelTracks *
#                     process.defaultAnalysis_05 *
#                     process.defaultAnalysis_510 *
#                     process.defaultAnalysis_1020 *
#                     process.defaultAnalysis_2030 *
                     process.defaultAnalysis_3040 *
#                     process.defaultAnalysis_4050 *
#                     process.defaultAnalysis_5060 *
                     process.defaultAnalysis_6070 
#                     process.defaultAnalysis_7080
)

process.schedule = cms.Schedule(process.p)

