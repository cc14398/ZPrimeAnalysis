# ZPrimeAnalysis

Based off of the exercises for the E/gamma short exercise for CMS DAS 2020 (https://github.com/cms-egamma/EgammaDAS2020)


Setup for EGamma exercises (needs a touch of modification to include this code):


#in your working area

source /cvmfs/cms.cern.ch/cmsset_default.sh 

export SCRAM_ARCH=slc7_amd64_gcc700

scram proj -n CMSSW_1068_CMSDAS CMSSW CMSSW_10_6_8

cd CMSSW_1068_CMSDAS/src

cmsenv

git cms-init

git clone git@github.com:cms-egamma/EgammaDAS2020.git EgammaUser/EgammaDAS2020 

git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools

#normally we check out the master but there was a fix put in for this school and it is not yet merged

#hence we created a separate temporary branch which will disappear soon

git clone git@github.com:cms-egamma/EgammaPostRecoTools.git EgammaUser/EgammaPostRecoTools -b CMSDAS2020

#small fix for the scales and smearing

git cms-merge-topic jainshilpi:ULV1_backport106X_forUsers

git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b UL2017SSV2 EgammaAnalysis/ElectronTools/data/

scram b -j 4
