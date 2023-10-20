from CRABClient.UserUtilities import config # , getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'PbPb2023_miniAOD_C2_32bin_0p5pttrg4p0_0p5ptass4p0_etamod2p4_cent080_Mmix10_374810_Oct19_2023_new3'
config.General.workArea = 'PbPb2023_miniAOD_nonsym_C2_32bin_0p5pttrg4p0_0p5ptass4p0_etamod2p4_cent080_Maximmix10_run374810_Oct19_2023'
config.General.transferOutputs = True
config.General.transferLogs = False


config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 4000
config.JobType.psetName = '../cfg/SayanCMW_central_newcfg_cent080_temp.py'
config.JobType.allowUndistributedCMSSW =True
config.JobType.inputFiles = ['../cfg/CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']


config.section_('Data')
config.Data.unitsPerJob = 1        #LumiBased  #40 is good
config.Data.totalUnits = -1
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#heavy ion quota
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/sayan/C2_Jet'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/user/sayan/'
config.Data.publication = False
config.Data.inputDataset = '/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD'
config.Data.runRange = '374810'
#config.Data.userInputFiles = open('/eos/cms/store/group/phys_heavyions/sayan/TrackingTools_run3_datataking/pbpb_miniaod_datasets/run374810_ls0001_ls1752_streamPhysicsHIPhysicsRawPrime0_miniaod_crab.txt').readlines()
config.Data.lumiMask = '/eos/cms/store/group/phys_heavyions/sayan/HIN_run3_pseudo_JSON/HIPhysicsRawPrime/Golden_Online_live.json'
#config.Data.lumiMask = 'run374810.json'
config.Data.outputDatasetTag = 'PbPb2023_miniAOD_nonsym_C2_32bin_0p5pttrg4p0_0p5ptass4p0_etamod2p4_cent080_Maximmix10_run374810_Oct19_2023_new3'


config.section_('Site')
config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.whitelist = ['T2_CH_CERN']




