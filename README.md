# Run3_DATA_MINIAOD

To setup the code, please, follow:

ssh -XY username@lxplus8.cern.ch

export SCRAM_ARCH=el8_amd64_gcc11

voms-proxy-init -voms cms

cmsrel CMSSW_13_2_5_patch1

cd CMSSW_13_2_5_patch1/src

cmsenv

git clone -b CMSSW_13_2_5_patch1_Data_Run3 https://github.com/sayanchatterjee38/Run3_DATA_MINIAOD.git

cd Run3_DATA_MINIAOD

mv * ../

cd ../

rm -rf Run3_DATA_MINIAOD 

scram build clean

scram b -j 12

To run the code locally, go to inside cfg:  cd Analyzers/SayanCMW/cfg

for generalAndHiPixelTracks (genAndHiPix): use SayanCMW_central_newcfg_cent080_genAndHiPixTrk.py

emacs SayanCMW_central_newcfg_cent080_genAndHiPixTrk.py

Verify : 

process.defaultCPDC.genTrk = False  ## It will always be False when running genAndHiPixel tracks

process.defaultCPDC.genAndHiPixTrk = True  ## It will always be True when running genAndHiPixel tracks 

process.defaultCPDC.ptmin_trigg = 0.3    ## minimum trigger pT cut, for genAndHiPix tracks, and it should not be less than 0.3 GeV/c

process.defaultCPDC.ptmin_ass = 0.3    ## minimum associate pT cut, for genAndHiPix tracks, and it should not be less than 0.3 GeV/c

process.defaultCPDC.ptmax_trigg = 4.0  ## maximum trigger pT cut, for both genAndHiPix tracks and genTrk, and it should not be greater than 16.0 GeV/c

process.defaultCPDC.ptmax_ass = 4.0  ## maximum associate pT cut, for both genAndHiPix tracks and genTrk, and it should not be greater than 16.0 GeV/c


do pt binning: Not more than 20 binning each (trigger_pt, associate_pt)

## Important : trigger_ptbin should be equal to associate_ptbin

process.defaultCPDC.trigger_ptbin = cms.untracked.vdouble(0.3, 1.0, 2.0, 3.0, 4.0)   ##trigger_ptbin[0] = ptmin_trigg   ##trigger_ptbin[trigger_ptbin.size() -1] = ptmax_trigg   

process.defaultCPDC.associate_ptbin = cms.untracked.vdouble(0.3, 1.0, 2.0, 3.0, 4.0)  ##associate_ptbin[0] = ptmin_ass   ##associate_ptbin[associate_ptbin.size() -1] = ptmax_ass



for generalTracks (genTrk): use SayanCMW_central_newcfg_cent080_genTrk.py

emacs SayanCMW_central_newcfg_cent080_genTrk.py

Verify : 

process.defaultCPDC.genTrk = True  ## It will always be True when running general tracks

process.defaultCPDC.genAndHiPixTrk = False  ## It will always be False when running general tracks 

process.defaultCPDC.ptmin_trigg = 0.5    ## minimum trigger pT cut, for general tracks, and it should not be less than 0.5 GeV/c

process.defaultCPDC.ptmin_ass = 0.5    ## minimum associate pT cut, for general tracks, and it should not be less than 0.5 GeV/c

process.defaultCPDC.ptmax_trigg = 4.0  ## maximum trigger pT cut, for both genAndHiPix tracks and genTrk, and it should not be greater than 16.0 GeV/c

process.defaultCPDC.ptmax_ass = 4.0  ## maximum associate pT cut, for both genAndHiPix tracks and genTrk, and it should not be greater than 16.0 GeV/c

do pt binning: Not more than 20 binning each (trigger_pt, associate_pt)

## Important : trigger_ptbin should be equal to associate_ptbin

process.defaultCPDC.trigger_ptbin = cms.untracked.vdouble(0.5, 1.0, 2.0, 3.0, 4.0)   ##trigger_ptbin[0] = ptmin_trigg   ##trigger_ptbin[trigger_ptbin.size() -1] = ptmax_trigg   

process.defaultCPDC.associate_ptbin = cms.untracked.vdouble(0.5, 1.0, 2.0, 3.0, 4.0)  ##associate_ptbin[0] = ptmin_ass   ##associate_ptbin[associate_ptbin.size() -1] = ptmax_ass



To run on Crab: go to inside crab,  cd Analyzers/SayanCMW/crab

for generalAndHiPixelTracks (genAndHiPix): use SayanCMW_crabJobs_ptbin_cent080_genAndHiPixTrk_run3_DAS.py

emacs SayanCMW_crabJobs_ptbin_cent080_genAndHiPixTrk_run3_DAS.py

Verify :   config.JobType.psetName = '../cfg/SayanCMW_central_newcfg_cent080_genAndHiPixTrk.py'


for generalTracks (genTrk): use SayanCMW_crabJobs_ptbin_cent080_genTrk_run3_DAS.py

emacs SayanCMW_crabJobs_ptbin_cent080_genTrk_run3_DAS.py

Verify :   config.JobType.psetName = '../cfg/SayanCMW_central_newcfg_cent080_genTrk.py'


Submit crab Job for run 374810: Dataset: /HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD

crab submit -c SayanCMW_crabJobs_ptbin_cent080_genAndHiPixTrk_run3_DAS.py (for generalAndHiPixelTracks)

crab submit -c SayanCMW_crabJobs_ptbin_cent080_genTrk_run3_DAS.py (for generalTracks)
