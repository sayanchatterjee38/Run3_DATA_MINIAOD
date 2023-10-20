import FWCore.ParameterSet.Config as cms

from Analyzers.SayanCMW.SayanCMW_cfi import *

#### Standard analysis for pixel Rereco PbPb 2015 ####

#### centrality 0-5% ####
CPDC05 = defaultCPDC.clone()
CPDC05.centmin = cms.untracked.int32(0)
CPDC05.centmax = cms.untracked.int32(5) ##5
#### centrality 5-10% ####
CPDC510 = defaultCPDC.clone()
CPDC510.centmin = cms.untracked.int32(5)
CPDC510.centmax = cms.untracked.int32(10)
#### centrality 10-20% ####
CPDC1020 = defaultCPDC.clone()
CPDC1020.centmin = cms.untracked.int32(10)
CPDC1020.centmax = cms.untracked.int32(20)
#### centrality 15-20% ####
CPDC1520 = defaultCPDC.clone()
CPDC1520.centmin = cms.untracked.int32(15)
CPDC1520.centmax = cms.untracked.int32(20)
#### centrality 20-30% ####
CPDC2030 = defaultCPDC.clone()
CPDC2030.centmin = cms.untracked.int32(20)
CPDC2030.centmax = cms.untracked.int32(30)
#### centrality 30-40% ####
CPDC3040 = defaultCPDC.clone()
CPDC3040.centmin = cms.untracked.int32(30)
CPDC3040.centmax = cms.untracked.int32(40)
#### centrality 40-50% ####
CPDC4050 = defaultCPDC.clone()
CPDC4050.centmin = cms.untracked.int32(40)
CPDC4050.centmax = cms.untracked.int32(50)
#### centrality 50-70% ####
CPDC5070 = defaultCPDC.clone()
CPDC5070.centmin = cms.untracked.int32(50)
CPDC5070.centmax = cms.untracked.int32(70)
#### centrality 50-60% ####
CPDC5060 = defaultCPDC.clone()
CPDC5060.centmin = cms.untracked.int32(50)
CPDC5060.centmax = cms.untracked.int32(60)
#### centrality 60-70% ####
CPDC6070 = defaultCPDC.clone()
CPDC6070.centmin = cms.untracked.int32(60)
CPDC6070.centmax = cms.untracked.int32(70)

#### centrality 70-80% ####
CPDC7080 = defaultCPDC.clone()
CPDC7080.centmin = cms.untracked.int32(70)
CPDC7080.centmax = cms.untracked.int32(80)
#### centrality 70-90% ####
CPDC70100 = defaultCPDC.clone()
CPDC70100.centmin = cms.untracked.int32(70)
CPDC70100.centmax = cms.untracked.int32(100)

#### centrality 0-5% ####
CPDC010 = defaultCPDC.clone()
CPDC010.centmin = cms.untracked.int32(0)
CPDC010.centmax = cms.untracked.int32(10) 

CPDC6080 = defaultCPDC.clone()
CPDC6080.centmin = cms.untracked.int32(60)
CPDC6080.centmax = cms.untracked.int32(80)


CPDC1015 = defaultCPDC.clone()
CPDC1015.centmin = cms.untracked.int32(10)
CPDC1015.centmax = cms.untracked.int32(15) ##5                                                                                                                                                             

CPDC1520 = defaultCPDC.clone()
CPDC1520.centmin = cms.untracked.int32(15)
CPDC1520.centmax = cms.untracked.int32(20)

CPDC2025 = defaultCPDC.clone()
CPDC2025.centmin = cms.untracked.int32(20)
CPDC2025.centmax = cms.untracked.int32(25) ##5                                                                                                                                                             

CPDC2530 = defaultCPDC.clone()
CPDC2530.centmin = cms.untracked.int32(25)
CPDC2530.centmax = cms.untracked.int32(30) ##5                                                                                                                                                             

CPDC3035 = defaultCPDC.clone()
CPDC3035.centmin = cms.untracked.int32(30)
CPDC3035.centmax = cms.untracked.int32(35) ##5                                                                                                                                                             

CPDC3540 = defaultCPDC.clone()
CPDC3540.centmin = cms.untracked.int32(35)
CPDC3540.centmax = cms.untracked.int32(40) ##5                                                                                                                                                             

CPDC4045 = defaultCPDC.clone()
CPDC4045.centmin = cms.untracked.int32(40)
CPDC4045.centmax = cms.untracked.int32(45) ##5                                                                                                                                                             

CPDC4550 = defaultCPDC.clone()
CPDC4550.centmin = cms.untracked.int32(45)
CPDC4550.centmax = cms.untracked.int32(50) ##5                                                                                                                                                             

CPDC5055 = defaultCPDC.clone()
CPDC5055.centmin = cms.untracked.int32(50)
CPDC5055.centmax = cms.untracked.int32(55) ##5                                                                                                                                                             

CPDC5560 = defaultCPDC.clone()
CPDC5560.centmin = cms.untracked.int32(55)
CPDC5560.centmax = cms.untracked.int32(60) ##5                                                                                                                                                             

CPDC6065 = defaultCPDC.clone()
CPDC6065.centmin = cms.untracked.int32(60)
CPDC6065.centmax = cms.untracked.int32(65) ##5                                                                                                                                                             

CPDC6570 = defaultCPDC.clone()
CPDC6570.centmin = cms.untracked.int32(65)
CPDC6570.centmax = cms.untracked.int32(70) ##5                                                                                                                                                             

CPDC7075 = defaultCPDC.clone()
CPDC7075.centmin = cms.untracked.int32(70)
CPDC7075.centmax = cms.untracked.int32(75) ##5                                                                                                                                                             

CPDC7580 = defaultCPDC.clone()
CPDC7580.centmin = cms.untracked.int32(75)
CPDC7580.centmax = cms.untracked.int32(80) ##5                                                                                                                                                             

CPDC8085 = defaultCPDC.clone()
CPDC8085.centmin = cms.untracked.int32(80)
CPDC8085.centmax = cms.untracked.int32(85) ##5                                                                                                                                                             

CPDC8590 = defaultCPDC.clone()
CPDC8590.centmin = cms.untracked.int32(85)
CPDC8590.centmax = cms.untracked.int32(90) ##5                                                                                                                                                             

