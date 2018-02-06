from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DY2017v2' ##
config.General.workArea = 'DY2017v2' ##
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RunGammaGammaMuMu_cfg.py'
config.JobType.inputFiles = ["MCPileupHighStats.root","MyDataPileupHistogram2017.root"]


config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM' ##
#config.Data.inputDataset = '/DoubleEG/Run2016C-PromptReco-v2/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.unitsPerJob = 50
#config.Data.lumiMask = '/afs/cern.ch/work/k/kshcheli/private/WW17/LL/CMSSW_9_4_0/src/JSON/combined_RPIN_CMS.json'
config.Data.outLFNDirBase = '/store/user/kshcheli/WW2017DileptonsWithZ/MC2/'
config.Data.publication = False
config.Data.outputDatasetTag = 'DY2017v2' ##

config.Site.storageSite = 'T2_CH_CERN'
