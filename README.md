# Run3_DATA_MINIAOD

To setup the code, please, follow:

ssh -XY username@lxplus8.cern.ch

export SCRAM_ARCH=el8_amd64_gcc11

voms-proxy-init -voms cms

cmsrel CMSSW_13_2_5_patch1

cd CMSSW_13_2_5_patch1/src

cmsenv

git clone -b CMSSW_13_2_X_trkAnalysis git@github.com:CesarBernardes/TrackingCode.git .

scram build clean

scram b -j 12
