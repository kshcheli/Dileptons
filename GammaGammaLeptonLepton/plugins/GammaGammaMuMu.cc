// -*- C++ -*-
//
// Package:    GammaGammaMuMu
// Class:      GammaGammaMuMu
// 
/**\class GammaGammaMuMu GammaGammaMuMu.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/src/GammaGammaMuMu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme,40 4-B20,+41227671567,
//         Created:  Thu Sep 13 15:17:14 CET 2012
// $Id: GammaGammaMuMu.cc,v 1.3 2013/04/28 08:40:45 lforthom Exp $
//
//

// WP80 cuts
const float MAX_MissingHits      = 0.0;
const float MIN_Dist             = 0.02;
const float MIN_Dcot             = 0.02;
const float cut_EB_trackRel03    = 0.09;
const float cut_EB_ecalRel03     = 0.07;
const float cut_EB_hcalRel03     = 0.10;
const float cut_EB_sigmaIetaIeta = 0.01;
const float cut_EB_deltaPhi      = 0.06;
const float cut_EB_deltaEta      = 0.004;
const float cut_EB_HoverE        = 0.04;
const float cut_EE_trackRel03    = 0.04;
const float cut_EE_ecalRel03     = 0.05;
const float cut_EE_hcalRel03     = 0.025;
const float cut_EE_sigmaIetaIeta = 0.03;
const float cut_EE_deltaPhi      = 0.03;
const float cut_EE_deltaEta      = 0.007; 
const float cut_EE_HoverE        = 0.025;

#include "/afs/cern.ch/work/k/kshcheli/private/WW17/LL/CMSSW_9_4_0/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMu.h"
//
// constructors and destructor
//
GammaGammaMuMu::GammaGammaMuMu(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  _fetchMuons = false;
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outfilename", "output.root");
  
  hltMenuLabel_ = iConfig.getParameter<std::string>("HLTMenuLabel");
  triggersList_ = iConfig.getParameter<std::vector<std::string> >("TriggersList");
  _hlts = new HLTmatches(triggersList_);
  nHLT = triggersList_.size();
	
  recoVertexLabel_ = iConfig.getParameter<edm::InputTag>("RecoVertexLabel");
  
  // Generator level
  sqrts_ = iConfig.getParameter<Double_t>("SqrtS");
  runOnMC_ = iConfig.getUntrackedParameter<bool>("RunOnMC", true);
  runOnProtons_ = iConfig.getUntrackedParameter<bool>("RunOnProtons", true);
  genLabel_ = iConfig.getParameter<edm::InputTag>("GenParticlesCollectionLabel");
  minPtMC_ = iConfig.getUntrackedParameter<Double_t>("MCAcceptPtCut", 20.);
  minEtaMC_ = iConfig.getUntrackedParameter<Double_t>("MCAcceptEtaCut", 2.5);
  
  // Pileup input tags

  pileupLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo", std::string("addPileupInfo"));
  mcPileupFile_ = iConfig.getUntrackedParameter<std::string>("mcpufile", "MCPileupHighStats.root");
  mcPileupPath_ = iConfig.getUntrackedParameter<std::string>("mcpupath", "h_trueNumInter");
  dataPileupFile_ = iConfig.getUntrackedParameter<std::string>("datapufile", "MyDataPileupHistogram2017.root");
  dataPileupPath_ = iConfig.getUntrackedParameter<std::string>("datapupath", "pileup");
//LumiWeights = new edm::LumiReWeighting("MCPileupHighStats.root","MyDataPileupHistogram0to75_MuonPhys.root","h_trueNumInter","pileup");
  
  // Leptons input tags
  leptonsType_ = iConfig.getParameter< std::vector<std::string> >("LeptonsType");
  for (i=0; i<leptonsType_.size(); i++) {
    if (leptonsType_[i]=="Muon") _fetchMuons = true;
  }
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalMuonCollectionLabel", std::string("muons"));
  //  rhoIsoLabel_ = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
  printCandidates_ = iConfig.getUntrackedParameter<bool>("PrintCandidates", false);
  trackLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("TrackCollectionLabel", std::string("generalTracks"));
  
  file = new TFile(outputFile_.c_str(), "recreate");
  file->cd();
  // tree definition
  tree = new TTree("ntp1", "ntp1");

  // PU reweighting
  if (runOnMC_) {
    LumiWeights = new edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_);
  }

//std::cout << "Still no idea"<<std::endl;


  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));  
  consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));  
  consumes<reco::MuonCollection>(edm::InputTag("muons"));   
  consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));    
  consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
  consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  //  tokenRPLocalTrack = consumes< edm::DetSetVector<TotemRPLocalTrack> >(edm::InputTag("tagLocalTrack")); 
 // tokenRPLocalTrack = consumes< edm::DetSetVector<TotemRPLocalTrack> >(edm::InputTag("totemRPLocalTrackFitter"));


//PPS 2017 stuff

  //tokenStatus_      =  consumes< edm::DetSetVector<TotemVFATStatus> >       ( iConfig.getParameter<edm::InputTag>( "tagStatus" ) ); 
  //tokenDigi_        = consumes< edm::DetSetVector<CTPPSDiamondDigi> >      ( iConfig.getParameter<edm::InputTag>( "tagDigi" ) ) ;
  tokenDiamondHit_  = consumes< edm::DetSetVector<CTPPSDiamondRecHit> >    ( iConfig.getParameter<edm::InputTag>( "tagDiamondRecHits" ) ) ;
  tokenDiamondTrack_= consumes< edm::DetSetVector<CTPPSDiamondLocalTrack> >( iConfig.getParameter<edm::InputTag>( "tagDiamondLocalTracks" ) ); 
  tokenLocalTrack_  = consumes< edm::DetSetVector<TotemRPLocalTrack> >     ( iConfig.getParameter<edm::InputTag>( "tagLocalTrack" ) ); 
  //tokenFEDInfo_     = consumes< std::vector<TotemFEDInfo> >                ( iConfig.getParameter<edm::InputTag>( "tagFEDInfo" ) ); 
  tokenPixelDigi_   = consumes<edm::DetSetVector<CTPPSPixelDigi> >(iConfig.getParameter<edm::InputTag>("tagRPixDigi") ) ;
  tokenPixelCluster_= consumes<edm::DetSetVector<CTPPSPixelCluster> >(iConfig.getParameter<edm::InputTag>("tagRPixCluster") ) ;
  tokenPixelRecHit_ =  consumes<edm::DetSetVector<CTPPSPixelRecHit> >(iConfig.getParameter<edm::InputTag>("tagRPixRecHit") ) ;
  tokenPixelLocalTrack_ =  consumes<edm::DetSetVector<CTPPSPixelLocalTrack> >(iConfig.getParameter<edm::InputTag>("tagRPixLocalTrack") ) ;




  //  consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices","","")); 

}


GammaGammaMuMu::~GammaGammaMuMu()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  file->Write();
  file->Close();

  //  delete _hlts;
  //  delete tree;

}


//
// member functions
//

void
GammaGammaMuMu::LookAtTriggersMuMu(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Int_t trigNum;
	
  // Get the trigger information from the event
  iEvent.getByLabel(edm::InputTag("TriggerResults","",hltMenuLabel_),hltResults_) ; 
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltResults_);

  for (unsigned int i=0; i<trigNames.size(); i++) {
    //std::cout << "--> " << trigNames.triggerNames().at(i) << std::endl;
    trigNum = _hlts->TriggerNum(trigNames.triggerNames().at(i));
    if (trigNum==-1) continue; // Trigger didn't match the interesting ones
    HLT_Accept[trigNum] = hltResults_->accept(i) ? 1 : 0;
    //    HLT_Prescl[trigNum] = hltConfig_.prescaleValue(iEvent, iSetup, trigNames.triggerNames().at(i));
    
    //LF FIXME need to think about that implementation...
    /*if (trigNames.triggerNames().at(i).find("CaloIdL")) {} // Leading lepton
      else if (trigNames.triggerNames().at(i).find("CaloIdT")) {} // Trailing lepton*/
    //std::cout << "*-------> " << trigNames.triggerNames().at(i).substr(0, trigNames.triggerNames().at(i).find_last_of("_"));
    //HLT_LeadingLepton_Prescl[] = hltConfig_.prescaleValue(event, iSetup, "HLT_Mu8Ele17L");  
  }
}

// ------------ method called for each event  ------------
void
GammaGammaMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //std::cout << "Beginning First init" << std::endl;

  // First initialization of the variables
  nCandidatesInEvent = 0;
  nPrimVertexCand = nFilteredPrimVertexCand = nMuVertexCand = 0;
  nMuonCand = nLeptonCand = 0;
  nExtraTracks = nQualityExtraTrack = 0;
  nJetCand = 0;
  nGenMuonCand = nGenMuonCandOutOfAccept = 0;
  nGenProtCand = 0;
  

  GenPair_p = GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = GenPair_Y = -999.;
  GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;

  closesttrkdxyz = closesthighpuritytrkdxyz = 999.;
  closestkalmantrkdxyz = 999.;
  for (i=0; i<MAX_LL; i++) {
    MuonCand_p[i] = MuonCand_px[i] = MuonCand_py[i] = MuonCand_pz[i] = -999.;
    MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = -999.;
    MuonCand_innerTrackPt[i] = MuonCand_innerTrackEta[i] = MuonCand_innerTrackPhi[i] = -999.;

    MuonCand_charge[i] = -999;
    MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
    MuonCand_npxlhits[i] = MuonCand_nstatseg[i] = MuonCand_ntrklayers[i] = -999;
    MuonCand_dxy[i] = MuonCand_dz[i] = -999.;
    MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = MuonCand_ispfmuon[i] = -999;
    MuonCand_istight[i] = -999;
    MuonCandTrack_nmuchits[i] = -999;
    MuonCandTrack_chisq[i] = -999.;
  }
  for (i=0; i<MAX_ET; i++) {
    ExtraTrack_p[i] = ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
    ExtraTrack_pt[i] = ExtraTrack_eta[i] = ExtraTrack_phi[i] = -999.;
    ExtraTrack_charge[i] = ExtraTrack_ndof[i] = -999;
    ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
    ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
    ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
  }
  for (i=0; i<MAX_PAIRS; i++) {
    Pair_candidates[i][0] = Pair_candidates[i][1] = -1;
    Pair_mindist[i] = Pair_p[i] = Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = Pair_Y[i] = -999.;
    Pair_dphi[i] = Pair_dpt[i] = Pair_3Dangle[i] = -999.;
    Pair_extratracks0pt5mm[i] = 0;
    Pair_extratracks1mm[i] = Pair_extratracks2mm[i] = Pair_extratracks3mm[i] = 0;
    Pair_extratracks4mm[i] = Pair_extratracks5mm[i] = Pair_extratracks1cm[i] = 0;
    Pair_extratracks2cm[i] = Pair_extratracks3cm[i] = Pair_extratracks4cm[i] = 0;
    Pair_extratracks5cm[i] = Pair_extratracks10cm[i] = 0;    
  }
  for (i=0; i<MAX_VTX; i++) {
    PrimVertexCand_id[i] = -1;
    PrimVertexCand_tracks[i] = PrimVertexCand_matchedtracks[i] = PrimVertexCand_unmatchedtracks[i] = 0;
    PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
    PrimVertexCand_chi2[i] = PrimVertexCand_ndof[i] = -999.;
  }

// PPS 2017

  nArmsTimingRecHits = 0;
  nArmsTimingTracks = 0;


  nTracksTiming = 0;
  nRecHitsTiming = 0;

  nArmsStrips = 0;
  nTracksStrips = 0;
  nPixelRecHits = 0;

  nArmsPixRecHits = 0;
  nLayersArm1PixRecHits = 0;
  nLayersArm2PixRecHits = 0;
  nPixelTracks = 0;
  nArmsPixelTracks = 0;
  nPixelTracksArm1 = 0;
  nPixelTracksArm2 = 0;


  int is45 = 0;
  int is56 = 0;
  int is45strips = 0;
  int is56strips = 0;

  int is45RecHitsTiming = 0;
  int is56RecHitsTiming = 0;


  int is45TimingTracks = 0;
  int is56TimingTracks = 0;

  int is45RecHitsPixel = 0;
  int is56RecHitsPixel = 0;

//  int layer1arm1 = 0; int layer2arm1 = 0; int layer3arm1 = 0; int layer4arm1 = 0;
//  int layer1arm2 = 0; int layer2arm2 = 0; int layer3arm2 = 0; int layer4arm2 = 0;
//  int nlayersarm1 = 0;
//  int nlayersarm2 = 0;

  unsigned int Channel = -1;
  unsigned int Arm = -1;
  double LE = 0.0; double TE = 0.0;
  unsigned int MH = -1;

  for(int i = 0; i < 100; i++)
    {
      ArmTiming[i] = -1;
      LeadingEdge[i] = 0;
      TrailingEdge[i] = 0;
      ToT[i] = 0;
      MultiHit[i] = 0;
      OOTIndex[i] = 0;
      ChannelTiming[i] = -999;
      PlaneTiming[i] = -999;
      XTiming[i] = -999;
      YTiming[i] = -999;
      TimingTrackT[i] = 0;
      TimingTrackTErr[i] = 0;
      TimingTrackX[i] = -999; 
      TimingTrackY[i] = -999;
      TimingTrackZ[i] = -999;
      TimingTrackArm[i] = -1;
      TimingTrackOOTIndex[i] = -999;
      TimingTrackMultiHit[i] = -999;
      TimingTrackChi2[i] = -999;
      TimingRecHitT[i] = -999;
      TimingRecHitX[i] = -999;
      TimingRecHitY[i] = -999;
      TimingRecHitOOTIndex[i] = -999; 
      TimingRecHitMultiHit[i] = -999;
      TimingRecHitToT[i] = 0;
      TimingRecHitChannel[i] = -999;
      TimingRecHitArm[i] = -999;
      TimingRecHitPlane[i] = -999;
      ArmStrips[i] = -1;
      StripTrackX[i] = 0; 
      StripTrackY[i] = 0;
      StripTrackZ[i] = 0;
      StripTrackTx[i] = 0;
      StripTrackTy[i] = 0;

  StationStrips[i] = 0;
  RPStrips[i] = 0;
  StripsMysteryFlag[i] = 0;

    }

  for(int i = 0; i < 1000; i++)
    {
      PixRecHitX[i] = 0;
      PixRecHitY[i] = 0;
      PixRecHitZ[i] = 0;
      PixRecHitArm[i] = -1;
      PixRecHitPlane[i] = -1;
      PixTrackX[i] = 0;
      PixTrackY[i] = 0;
      PixTrackTx[i] = 0;
      PixTrackTy[i] = 0;
      PixTrackChi2[i] = 0;
      PixTrackZ[i] = 0;
      PixTrackArm[i] = -1;
    }




  // JH
  KalmanVertexCand_x = KalmanVertexCand_y = KalmanVertexCand_z = -999.;
  KalmanVertexCand_chi2 = -999.;
  KalmanVertexCand_ndof = 0;
  ClosestExtraTrackKalman_vtxdxyz = 999.;
  
  Weight = 1.;
  
  foundPairInEvent = false;
  foundGenCandPairInEvent = false;
  
  muonsMomenta.clear();
  electronsMomenta.clear();
  
  _leptonptmp = new TLorentzVector();
  
  // Run and BX information
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();
  
  // High level trigger information retrieval  
  LookAtTriggersMuMu(iEvent, iSetup);
  
  // beam spot information
  //  iEvent.getByLabel(beamSpotLabel_, beamspot_h);
  //  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // Get the vertex collection from the event
  iEvent.getByLabel(recoVertexLabel_, recoVertexColl);
  const reco::VertexCollection* vertices = recoVertexColl.product();

  // rho for isolation 
  //  iEvent.getByLabel(rhoIsoLabel_, rhoIso_h); 
  //  rhoIso = *(rhoIso_h.product()); 
  
  //std::cout << "Passed Isolation" << std::endl;


  // Generator level information
  if (runOnMC_) {
    iEvent.getByLabel(genLabel_, genPartColl);
    
    for (genPart=genPartColl->begin(); genPart!=genPartColl->end(); genPart++) {
      if (genPart->pt()<minPtMC_ || (minEtaMC_!=-1. && fabs(genPart->eta())>minEtaMC_)) {
        if (fabs(genPart->pdgId())==13) nGenMuonCandOutOfAccept++;
        if (fabs(genPart->pdgId())==22) nGenPhotCandOutOfAccept++;
        continue;
      }
      if (genPart->pdgId()!=2212 && genPart->status()!=1) continue;
      if (fabs(genPart->pdgId())==13 && nGenMuonCand<MAX_GENMU) {
        GenMuonCand_p[nGenMuonCand] = genPart->p();
        GenMuonCand_px[nGenMuonCand] = genPart->px();
        GenMuonCand_py[nGenMuonCand] = genPart->py();
        GenMuonCand_pz[nGenMuonCand] = genPart->pz();
        GenMuonCand_pt[nGenMuonCand] = genPart->pt();
        GenMuonCand_eta[nGenMuonCand] = genPart->eta();
        GenMuonCand_phi[nGenMuonCand] = genPart->phi();
        
        nGenMuonCand++;
      }

      if (genPart->pdgId()==2212 && nGenProtCand<MAX_GENPRO) { 
	GenProtCand_p[nGenProtCand] = genPart->p(); 
	GenProtCand_e[nGenProtCand] = genPart->energy();  
	GenProtCand_px[nGenProtCand] = genPart->px(); 
	GenProtCand_py[nGenProtCand] = genPart->py(); 
	GenProtCand_pz[nGenProtCand] = genPart->pz(); 
	GenProtCand_pt[nGenProtCand] = genPart->pt(); 
	GenProtCand_eta[nGenProtCand] = genPart->eta(); 
	GenProtCand_phi[nGenProtCand] = genPart->phi(); 
	GenProtCand_status[nGenProtCand] = genPart->status(); 
 
	double xi = 1 - ((genPart->energy())/(sqrts_/2.0)); 
	double t = -(std::pow(genPart->pt(), 2)); 
 
	GenProtCand_xi[nGenProtCand] = xi; 
	GenProtCand_t[nGenProtCand] = t; 
       
	nGenProtCand++; 
      } 


      foundGenCandPairInEvent = false;
      if (_fetchMuons) { // Looks at dimuons
        if(nGenMuonCand!=2) continue; // FIXME maybe a bit tight according to the newer PU conditions?
        l1.SetXYZM(GenMuonCand_px[0], GenMuonCand_py[0], GenMuonCand_pz[0], MASS_MU);      	
        l2.SetXYZM(GenMuonCand_px[1], GenMuonCand_py[1], GenMuonCand_pz[1], MASS_MU);      	
        foundGenCandPairInEvent = true;      	
      }
      if (foundGenCandPairInEvent) {
        pair = l1+l2;
        GenPair_p = pair.P();
        GenPair_pt = pair.Pt();
        GenPair_mass = pair.M();
        GenPair_phi = pair.Phi();
        GenPair_eta = pair.Eta();
	GenPair_Y = pair.Rapidity();
        dphi = fabs(l1.Phi()-l2.Phi());
        GenPair_dphi = (dphi<pi) ? dphi : 2.*pi-dphi; // dphi lies in [-pi, pi]
        GenPair_dpt = fabs(l1.Pt()-l2.Pt());
        GenPair_3Dangle = (l1.Angle(l2.Vect()))/pi;
      }
      //      if(genPart->pdgId()==2212 && fabs(genPart->pz())>3000.) {
        // Kinematic quantities computation
        // xi = fractional momentum loss
	//        if (genPart->pz()>0.) xi = 1.-genPart->pz()/sqrts_;
	//        else xi = 1.+genPart->pz()/sqrts_;
	//        t = -(std::pow(genPart->pt(), 2)+std::pow(MASS_P*xi, 2))/(1.-xi);
      //      }
    }
  }

  // Pileup information
  if (runOnMC_) {
    iEvent.getByLabel(pileupLabel_, pileupInfo);

    // This part is optional if the distributions are already generated
    sum_nvtx = 0;
    npv = npvtrue = npvm1true = npvp1true = npv0true = npv0 = 0;

    for(PVI=pileupInfo->begin(); PVI!=pileupInfo->end(); PVI++) {
      beamXing = PVI->getBunchCrossing();
      npv = PVI->getPU_NumInteractions();
      npvtrue = PVI->getTrueNumInteractions();
      sum_nvtx += npvtrue;

      if(beamXing == -1) npvm1true+=npvtrue;
      if(beamXing == 0) {
        npv0 += npv;
        npv0true += npvtrue;
      }
      if(beamXing == 1) npvp1true+=npvtrue;
    }
    nTruePUforPUWeightBX0 = npv0true;
    //

    //    std::cout << "\t\tJH: Creating iEventB" << std::endl;
    //    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    //   Weight = LumiWeights->weight(*iEventB);
    Weight = LumiWeights->weight( nTruePUforPUWeightBX0 );
  }
  
  std::vector<reco::TransientTrack> translepttrks;  
  reco::TrackCollection * lepttrks = new reco::TrackCollection;  
  edm::ESHandle<TransientTrackBuilder> theKalVtx;  
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theKalVtx);  
  iEvent.getByLabel(trackLabel_, trackColl);  

  // RP
  if(runOnProtons_)
    {

// PPS 2017
/////////////////////////////////////////////////////////
/*
edm::Handle< edm::DetSetVector<TotemVFATStatus> > diamondVFATStatus;
  iEvent.getByToken( tokenStatus_, diamondVFATStatus );

  edm::Handle< edm::DetSetVector<CTPPSDiamondDigi> > diamondDigis;
  iEvent.getByToken( tokenDigi_, diamondDigis );

  edm::Handle< std::vector<TotemFEDInfo> > fedInfo;
  iEvent.getByToken( tokenFEDInfo_, fedInfo );
*/
  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > diamondRecHits;
  iEvent.getByToken( tokenDiamondHit_, diamondRecHits );

  edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > diamondLocalTracks;
  iEvent.getByToken( tokenDiamondTrack_, diamondLocalTracks );

  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > stripTracks;
  iEvent.getByToken( tokenLocalTrack_, stripTracks );

/*
  edm::Handle< edm::DetSetVector<CTPPSPixelDigi> > pixDigi;
  iEvent.getByToken(tokenPixelDigi_, pixDigi);

  edm::Handle< edm::DetSetVector<CTPPSPixelCluster> > pixClusters;
  iEvent.getByToken(tokenPixelCluster_, pixClusters);
*/
  edm::Handle< edm::DetSetVector<CTPPSPixelRecHit> > pixRecHits;
  iEvent.getByToken(tokenPixelRecHit_, pixRecHits);

  edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixLocalTracks;
  iEvent.getByToken(tokenPixelLocalTrack_, pixLocalTracks);


  /* RecHits - timing */
  for ( const auto& rechits_ds : *diamondRecHits ) {
    const CTPPSDiamondDetId detidforrh( rechits_ds.detId() );
    for ( const auto& rechit : rechits_ds ) {

      TimingRecHitArm[nRecHitsTiming] = detidforrh.arm();
      if(TimingRecHitArm[nRecHitsTiming] == 0)
	is45RecHitsTiming = 1;
      if(TimingRecHitArm[nRecHitsTiming] == 1)
	is56RecHitsTiming = 1;

      TimingRecHitChannel[nRecHitsTiming] = detidforrh.channel();      
      TimingRecHitPlane[nRecHitsTiming] = detidforrh.plane();
      TimingRecHitT[nRecHitsTiming] = rechit.getT();
      TimingRecHitX[nRecHitsTiming] = rechit.getX();                                                                                                            
      TimingRecHitY[nRecHitsTiming] = rechit.getY();                                                                                                            
      TimingRecHitOOTIndex[nRecHitsTiming] = rechit.getOOTIndex();                                                                                              
      TimingRecHitMultiHit[nRecHitsTiming] = rechit.getMultipleHits();
      TimingRecHitToT[nRecHitsTiming] = rechit.getToT();
      nRecHitsTiming++;
    }
  }


  nArmsTimingRecHits = is45RecHitsTiming+is56RecHitsTiming;


  /* Diamond tracks */
  for ( const auto& ds2 : *diamondLocalTracks ) 
    {
      for ( const auto& tr2 : ds2 ) 
	{
	  if ( ! tr2.isValid() ) continue;
	  
	  CTPPSDetId diamId2( ds2.detId() );
	  unsigned int arm1 = diamId2.arm();
     
	  //	  unsigned int arm2 = diamId2.arm();
	  TimingTrackT[nTracksTiming] = tr2.getT();
	  TimingTrackTErr[nTracksTiming] = tr2.getTSigma();
	  TimingTrackX[nTracksTiming] = tr2.getX0();
	  TimingTrackY[nTracksTiming] = tr2.getY0();
	  TimingTrackZ[nTracksTiming] = tr2.getZ0();
	  TimingTrackOOTIndex[nTracksTiming] = tr2.getOOTIndex();
	  TimingTrackMultiHit[nTracksTiming] = tr2.getMultipleHits();
	  TimingTrackChi2[nTracksTiming] = tr2.getChiSquared();
	  TimingTrackArm[nTracksTiming] = arm1;

	  nTracksTiming++;


      if(arm1 == 0)
	is45TimingTracks = 1;
      if(arm1 == 1)
	is56TimingTracks = 1;
	}
    }
nArmsTimingTracks = is45TimingTracks+is56TimingTracks;

  /* Digis - timing */
/*
  Channel = -1;                                                                                                                                    
  Arm = -1;                                                                                                                                        
  LE = 0;                                                                                                                                          
  TE = 0;                                                                                                                                          
  MH = 0;                                                                                                                                          

  for ( const auto& digis : *diamondDigis ) {                                                                                                      
    const CTPPSDiamondDetId detId( digis.detId() );                                                                                                
    CTPPSDiamondDetId detId_pot( digis.detId() );                                                                                                  
    for ( const auto& digi : digis ) {                                                                                                             
      detId_pot.setPlane( 0 );                                                                                                                     
      detId_pot.setChannel( 0 );                                                                                                                   

      LE = digi.getLeadingEdge() * 0.025;                                                                                                          
      TE = digi.getTrailingEdge() * 0.025;                                                                                                         

      Channel = detId.channel();                                                                                                                   
      MH = digi.getMultipleHit();                                                                                                                  
      Arm = detId.arm();                                                                                                                           

      // Clock
      if(Channel != 30)
	{
	  ArmTiming[nHitsTiming] = Arm;                                                                                                            
	  if(ArmTiming[nHitsTiming] == 0 && LE > 0)                                                                                                
	    {                                                                                                                                       
	      is45=1;                                                                                                                              
	      if(detId.plane() == 0) layer1arm1 = 1;                                                                                               
	      if(detId.plane() == 1) layer2arm1 = 1;                                                                                               
	      if(detId.plane() == 2) layer3arm1 = 1;                                                                                               
	      if(detId.plane() == 3) layer4arm1 = 1;                                                                                               
	    }                                                                                                                                      
 	  
	  if(ArmTiming[nHitsTiming] == 1 && LE > 0)                                                                                                
	    {                                                                                                                                      
	      is56=1;                                                                                                                              
	      if(detId.plane() == 0) layer1arm2 = 1;                                                                                               
	      if(detId.plane() == 1) layer2arm2 = 1;                                                                                               
	      if(detId.plane() == 2) layer3arm2 = 1;                                                                                               
	      if(detId.plane() == 3) layer4arm2 = 1;                                                                                               
	    }                                                                                                                                      
	  
	  ChannelTiming[nHitsTiming] = Channel;                                                                                                    
	  PlaneTiming[nHitsTiming] = detId.plane();                                                                                                
	  LeadingEdge[nHitsTiming] = LE;                                                                                                           
             
	  TrailingEdge[nHitsTiming] = TE;
	  ToT[nHitsTiming] = TE - LE;
	  MultiHit[nHitsTiming] = MH;
	  //      OOTIndex[nHitsTiming] = rechit.getOOTIndex();                                                                                 

	  nHitsTiming++;
	}         
    }
  }                                                                                                                                                
  */
 
  /* Strips */
  for ( const auto& ds1 : *stripTracks ) {
    for ( const auto& tr1 : ds1 ) {
      if ( ! tr1.isValid() )  continue;

      CTPPSDetId rpId1( ds1.detId() );
      unsigned int arm1 = rpId1.arm();
      unsigned int stNum1 = rpId1.station();
      unsigned int rpNum1 = rpId1.rp();
      if (stNum1 != 0 || ( rpNum1 != 2 && rpNum1 != 3 ) ){ StripsMysteryFlag[nTracksStrips]=1;}

      StripTrackX[nTracksStrips] = (tr1.getX0())/1000.;
      StripTrackY[nTracksStrips] = (tr1.getY0())/1000.;
      StripTrackZ[nTracksStrips] = (tr1.getZ0())/1000.;

      StripTrackTx[nTracksStrips] = tr1.getTx();
      StripTrackTy[nTracksStrips] = tr1.getTy();
      ArmStrips[nTracksStrips] = arm1;
      StationStrips[nTracksStrips] = stNum1;
      RPStrips[nTracksStrips] = rpNum1;

      if(arm1 == 0)
	is45strips = 1;
      if(arm1 == 1)
	is56strips = 1;

      nTracksStrips++;
    }
  }

  nArmsStrips = is45strips+is56strips;

  /* Pixel digis */
/*
  int pixelsarm1 = 0;
  int pixelsarm2 = 0;
  int planespixelsarm1 = 0;
  int planespixelsarm2 = 0;

  int is45pixels = 0;
  int is56pixels = 0;
  int prevplane1 = -1;
  int prevplane2 = -1;

  if(pixDigi.isValid())
    {
      for(const auto &ds_digi : *pixDigi)
        {
          int idet = (ds_digi.id>>DetId::kDetOffset)&0xF;
          if(idet != DetId::VeryForward) {
            continue;
          }

          int plane = ((ds_digi.id>>16)&0x7);

          CTPPSDetId theId(ds_digi.id);
          int pixelsarm = theId.arm()&0x1;

          if(pixelsarm == 0)
            {
              if(plane != prevplane1)
                {
                  planespixelsarm1++;
                  prevplane1 = plane;
                }
              pixelsarm1++;
	      is45pixels = 1;
            }
          if(pixelsarm == 1)
            {
              if(plane != prevplane2)
                {
                  planespixelsarm2++;
                  prevplane2 = plane;
                }
              pixelsarm2++;
	      is56pixels = 1;
            }
	  //          int station = theId.station()&0x3;
	  //          int rpot = theId.rp()&0x7;
        }
    }

  nArmsPixelDigis = is45pixels + is56pixels; 
  nLayersArm1PixelDigis = planespixelsarm1;
  nLayersArm2PixelDigis = planespixelsarm2;
*/
  /* Pixel RecHits */
  int pixelsrharm1 = 0;
  int pixelsrharm2 = 0;
  int planespixelsrharm1 = 0;
  int planespixelsrharm2 = 0;

  int is45pixelsrh = 0;
  int is56pixelsrh = 0;
  unsigned int pixelsrhprevplane1 = -1;
  unsigned int pixelsrhprevplane2 = -1;

  if(pixRecHits.isValid())
    {
      for ( const auto& rechits_px : *pixRecHits ) 
	{
	  const CTPPSPixelDetId pxDetid( rechits_px.detId() );
	  for ( const auto& rechitpx : rechits_px ) 
	    {
	      PixRecHitX[nPixelRecHits] = rechitpx.getPoint().x();
	      PixRecHitY[nPixelRecHits] = rechitpx.getPoint().y();
	      PixRecHitZ[nPixelRecHits] = rechitpx.getPoint().z();
	      PixRecHitArm[nPixelRecHits] = pxDetid.arm();
              PixRecHitPlane[nPixelRecHits] = pxDetid.plane();
	      nPixelRecHits++;

	      if(pxDetid.arm() == 0)
		{
		  if(pxDetid.plane() != pixelsrhprevplane1)
		    {
		      planespixelsrharm1++;
		      pixelsrhprevplane1 = pxDetid.plane();
		    }
		  pixelsrharm1++;
		  is45pixelsrh = 1;
		}
	      if(pxDetid.arm() == 1)
		{
		  if(pxDetid.plane() != pixelsrhprevplane2)
		    {
		      planespixelsrharm2++;
		      pixelsrhprevplane2 = pxDetid.plane();
		    }
		  pixelsrharm2++;
		  is56pixelsrh = 1;
		}
	    }
	}
    }

  nArmsPixRecHits = is45pixelsrh + is56pixelsrh;
  nLayersArm1PixRecHits = planespixelsrharm1;
  nLayersArm2PixRecHits = planespixelsrharm2;

  /* Pixel tracks */
  for ( const auto& dspxtr1 : *pixLocalTracks ) {
    for ( const auto& pxtr1 : dspxtr1 ) {
      if ( ! pxtr1.isValid() )  continue;

      CTPPSDetId pxrpId1( dspxtr1.detId() );
      unsigned int pxarm1 = pxrpId1.arm();

      PixTrackX[nPixelTracks] = pxtr1.getX0();
      PixTrackY[nPixelTracks] = pxtr1.getY0();
      PixTrackTx[nPixelTracks] = pxtr1.getTx();
      PixTrackTy[nPixelTracks] = pxtr1.getTy();
      PixTrackChi2[nPixelTracks] = pxtr1.getChiSquared();
      PixTrackZ[nPixelTracks] = pxtr1.getZ0();
      PixTrackArm[nPixelTracks] = pxarm1;

      if(pxarm1 == 0)
	nPixelTracksArm1++;
      if(pxarm1 == 1)
	nPixelTracksArm2++;

      nPixelTracks++;
    }
  }

nArmsPixelTracks = (nPixelTracksArm1 > 0) + (nPixelTracksArm2 > 0);

////////////////////////////////////////////////////////


      // Placeholder for more awesome stuff
/*
      edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rplocaltracks; 
      iEvent.getByToken(tokenRPLocalTrack, rplocaltracks); 
      
      for (auto &ds1 : *rplocaltracks)
	{
	  for (auto &tr1 : ds1)
	    {
	      if (! tr1.isValid())
		continue;
	      
	      LocalProtCand_x[nLocalProtCand] = (tr1.getX0())/1000.0;
	      LocalProtCand_y[nLocalProtCand] = (tr1.getY0())/1000.0; 
	      LocalProtCand_z[nLocalProtCand] = (tr1.getZ0())/1000.0; 
	      LocalProtCand_xSigma[nLocalProtCand] = (tr1.getX0Sigma())/1000.0;
              LocalProtCand_ySigma[nLocalProtCand] = (tr1.getY0Sigma())/1000.0; 
              LocalProtCand_Tx[nLocalProtCand] = tr1.getTx();  
              LocalProtCand_Ty[nLocalProtCand] = tr1.getTy();   
	      LocalProtCand_TxSigma[nLocalProtCand] = tr1.getTxSigma(); 
              LocalProtCand_TySigma[nLocalProtCand] = tr1.getTySigma();  

	      nLocalProtCand++; 
	    }
	}
*/
    }

  // Get the muons collection from the event
  if (_fetchMuons) {
    iEvent.getByLabel(muonLabel_, muonColl);
    //    std::cout << "Looping on muon coll of size = " << muonColl->size() << std::endl;
    for (muon=muonColl->begin(); muon!=muonColl->end() && nMuonCand<MAX_MUONS; muon++) {
      MuonCand_p[nMuonCand] = muon->p();
      MuonCand_px[nMuonCand] = muon->px();
      MuonCand_py[nMuonCand] = muon->py();
      MuonCand_pz[nMuonCand] = muon->pz();
      MuonCand_pt[nMuonCand] = muon->pt();
      MuonCand_eta[nMuonCand] = muon->eta();
      MuonCand_phi[nMuonCand] = muon->phi();
      MuonCand_charge[nMuonCand] = muon->charge();
      //      MuonCand_dxy[nMuonCand] = muon->dB();
      MuonCand_dxy[nMuonCand] = 0;
      MuonCand_nstatseg[nMuonCand] = muon->numberOfMatchedStations();
      
      MuonCand_isglobal[nMuonCand] = muon->isGlobalMuon();
      MuonCand_istracker[nMuonCand] = muon->isTrackerMuon();
      MuonCand_isstandalone[nMuonCand] = muon->isStandAloneMuon();
      MuonCand_ispfmuon[nMuonCand] = muon->isPFMuon();

      if((MuonCand_isglobal[nMuonCand] == 1) || (MuonCand_istracker[nMuonCand] == 1)) 
	{
	  MuonCand_innerTrackPt[nMuonCand] = muon->innerTrack()->pt(); 
	  MuonCand_innerTrackEta[nMuonCand] = muon->innerTrack()->eta();  
	  MuonCand_innerTrackPhi[nMuonCand] = muon->innerTrack()->phi();  

	  _leptonptmp->SetXYZM(muon->innerTrack()->px(), muon->innerTrack()->py(), muon->innerTrack()->pz(), muon->mass());
          // JH - use 2 highest pT muons for vertexing
	  if(nMuonCand < 2)
	    {
	      lepttrks->push_back (*(muon->innerTrack()));
	      reco::TransientTrack tmptrk = (*theKalVtx).build( *(muon->innerTrack() ) ); 
	      translepttrks.push_back( tmptrk ); 
	    }
	}
      else
	{
	  _leptonptmp->SetXYZM(muon->px(), muon->py(), muon->pz(), muon->mass());
	  MuonCand_innerTrackPt[nMuonCand] = -999.;
	  MuonCand_innerTrackEta[nMuonCand] = -999.;
	  MuonCand_innerTrackPhi[nMuonCand] = -999.;
	}

      muonsMomenta.insert(std::pair<Int_t,TLorentzVector>(nMuonCand, *_leptonptmp));
      
      MuonCand_vtxx[nMuonCand] = muon->vertex().x();
      MuonCand_vtxy[nMuonCand] = muon->vertex().y();
      MuonCand_vtxz[nMuonCand] = muon->vertex().z();
      
      if (MuonCand_istracker[nMuonCand]) {
	MuonCand_npxlhits[nMuonCand] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
	MuonCand_ntrklayers[nMuonCand] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      }

      if (MuonCand_isglobal[nMuonCand] && MuonCand_istracker[nMuonCand]) {
	MuonCandTrack_nmuchits[nMuonCand] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
	MuonCandTrack_chisq[nMuonCand] = muon->globalTrack()->normalizedChi2();
	istight = true;
	istight&= MuonCand_ispfmuon[nMuonCand];
	istight&= (MuonCandTrack_chisq[nMuonCand]<10.);
	istight&= (MuonCandTrack_nmuchits[nMuonCand]>=1);
	istight&= (MuonCand_nstatseg[nMuonCand]>=2);
	istight&= (MuonCand_dxy[nMuonCand]<.2);
	istight&= (MuonCand_npxlhits[nMuonCand]>0);
	istight&= (MuonCand_ntrklayers[nMuonCand]>5);
	MuonCand_istight[nMuonCand] = istight;

      } 

      nMuonCand++;
    }
  }
  
  // JH 
  // If 2 muons, make a vertex and compute pair quantities from the leading 2
  // For now we do this "outside-in", rather than starting from the list of offlinePrimaryVertices and looking for muons
  if(translepttrks.size() >= 2)
    {
      l1.SetXYZM(MuonCand_px[0], 
		 MuonCand_py[0],
		 MuonCand_pz[0],
		 MASS_MU); 
      l2.SetXYZM(MuonCand_px[1],
		 MuonCand_py[1],
		 MuonCand_pz[1],
		 MASS_MU); 

      pair = l1+l2; 

      Pair_p[0] = pair.P(); 
      Pair_pt[0] = pair.Pt(); 
      Pair_mass[0] = pair.M(); 
      Pair_phi[0] = pair.Phi(); 
      Pair_eta[0] = pair.Eta(); 
      Pair_Y[0] = pair.Rapidity();
      dphi = fabs(l1.Phi()-l2.Phi()); 
      Pair_dphi[0] = (dphi<pi) ? dphi : 2.*pi-dphi; // dphi lies in [-pi, pi] 
      Pair_dpt[0] = fabs(l1.Pt()-l2.Pt()); 
      Pair_3Dangle[0] = (l1.Angle(l2.Vect()))/pi; 
      nMuVertexCand++;

      KalmanVertexFitter fitter(true); 
      TransientVertex dileptonVertex = fitter.vertex(translepttrks); 
      if(dileptonVertex.isValid()) { 
	foundPairInEvent = true;
	KalmanVertexCand_x = dileptonVertex.position().x(); 
	KalmanVertexCand_y = dileptonVertex.position().y(); 
	KalmanVertexCand_z = dileptonVertex.position().z(); 
	KalmanVertexCand_chi2 = dileptonVertex.totalChiSquared();
	KalmanVertexCand_ndof = dileptonVertex.degreesOfFreedom();

	etind = 0; 

	//  Count nearby tracks
	for ( track = trackColl->begin(); track != trackColl->end() && etind<MAX_ET; ++track ) 
	  {
	    if((track->pt() == MuonCand_innerTrackPt[0]) || (track->pt() == MuonCand_innerTrackPt[1]))
	      continue;

	    vtxdst = sqrt(std::pow((track->vertex().x()-KalmanVertexCand_x),2)+ 
			  std::pow((track->vertex().y()-KalmanVertexCand_y),2)+ 
			  std::pow((track->vertex().z()-KalmanVertexCand_z),2)); 
           
	    nExtraTracks++;
	    
	    if (vtxdst<0.05) Pair_extratracks0pt5mm[0]++;
	    if (vtxdst<0.1) Pair_extratracks1mm[0]++; 
	    if (vtxdst<0.2) Pair_extratracks2mm[0]++; 
	    if (vtxdst<0.3) Pair_extratracks3mm[0]++; 
	    if (vtxdst<0.4) Pair_extratracks4mm[0]++; 
	    if (vtxdst<0.5) Pair_extratracks5mm[0]++; 
	    if (vtxdst<1.0) Pair_extratracks1cm[0]++; 
	    if (vtxdst<2.0) Pair_extratracks2cm[0]++; 
	    if (vtxdst<3.0) Pair_extratracks3cm[0]++; 
	    if (vtxdst<4.0) Pair_extratracks4cm[0]++; 
	    if (vtxdst<5.0) Pair_extratracks5cm[0]++; 
	    if (vtxdst<10.) Pair_extratracks10cm[0]++; 
	    if (vtxdst<closesttrkdxyz) { 
	      closesttrkdxyz = vtxdst; 
	    } 
	    if (track->quality(reco::TrackBase::highPurity)==1) { 
	      if (vtxdst<closesthighpuritytrkdxyz) { 
		closesthighpuritytrkdxyz = vtxdst; 
	      } 
	    } 

	    // Save track properties if within 5mm
	    if(vtxdst<0.5) 
	      {
		ExtraTrack_purity[etind] = track->quality(reco::TrackBase::highPurity);  
		ExtraTrack_nhits[etind] = track->numberOfValidHits();  
            
		ExtraTrack_p[etind] = track->p();  
		ExtraTrack_px[etind] = track->px();  
		ExtraTrack_py[etind] = track->py();  
		ExtraTrack_pz[etind] = track->pz();  
		ExtraTrack_pt[etind] = track->pt();  
		ExtraTrack_eta[etind] = track->eta();  
		ExtraTrack_phi[etind] = track->phi();  
		ExtraTrack_charge[etind] = track->charge();  
		ExtraTrack_chi2[etind] = track->chi2();  
		ExtraTrack_ndof[etind] = track->ndof();  
		ExtraTrack_vtxdxyz[etind] = vtxdst;  
		ExtraTrack_vtxT[etind] = sqrt(std::pow(track->vertex().x()-KalmanVertexCand_x,2)+  
					      std::pow(track->vertex().y()-KalmanVertexCand_y,2));  
		ExtraTrack_vtxZ[etind] = fabs(track->vertex().z()-KalmanVertexCand_z);  
		ExtraTrack_x[etind] = track->vertex().x();  
		ExtraTrack_y[etind] = track->vertex().y();  
		ExtraTrack_z[etind] = track->vertex().z();  

		etind++; 
	      }

	  }
	ClosestExtraTrack_vtxdxyz = closesttrkdxyz; 
	ClosestHighPurityExtraTrack_vtxdxyz = closesthighpuritytrkdxyz;
      } 
    }

  for (vertex=vertices->begin(); vertex!=vertices->end() && nPrimVertexCand<MAX_VTX; ++vertex) 
    { 
      PrimVertexCand_id[nPrimVertexCand] = nPrimVertexCand; 
      PrimVertexCand_x[nPrimVertexCand] = vertex->x(); 
      PrimVertexCand_y[nPrimVertexCand] = vertex->y(); 
      PrimVertexCand_z[nPrimVertexCand] = vertex->z(); 
      PrimVertexCand_tracks[nPrimVertexCand] = vertex->nTracks(); 
      PrimVertexCand_chi2[nPrimVertexCand] = vertex->chi2(); 
      PrimVertexCand_ndof[nPrimVertexCand] = vertex->ndof(); 
      nPrimVertexCand++;
    }

  // end JH 


  delete _leptonptmp;

  if(foundPairInEvent)
    {
      if((MuonCand_pt[0] >= 40) && (MuonCand_pt[1] >= 40))
	tree->Fill();
    }

	
  //	tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMu::beginJob()
{
  // Booking the ntuple

  tree->Branch("Run", &Run, "Run/I");
  tree->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree->Branch("BX", &BX, "BX/I");
  tree->Branch("EventNum", &EventNum, "EventNum/I");
  /*tree->Branch("AvgInstDelLumi", &AvgInstDelLumi, "AvgInstDelLumi/D");
  tree->Branch("BunchInstLumi", &BunchInstLumi, "BunchInstLumi[3]/D");*/

  tree->Branch("nHLT", &nHLT, "nHLT/I");
  tree->Branch("HLT_Accept", HLT_Accept, "HLT_Prescl[nHLT]/I");
  tree->Branch("HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I");
  
  if (_fetchMuons) {
    tree->Branch("nMuonCand", &nMuonCand, "nMuonCand/I");
    tree->Branch("MuonCand_px", MuonCand_px, "MuonCand_px[nMuonCand]/D");
    tree->Branch("MuonCand_py", MuonCand_py, "MuonCand_py[nMuonCand]/D");
    tree->Branch("MuonCand_pz", MuonCand_pz, "MuonCand_pz[nMuonCand]/D");
    tree->Branch("MuonCand_p", MuonCand_p, "MuonCand_p[nMuonCand]/D");
    tree->Branch("MuonCand_pt", MuonCand_pt, "MuonCand_pt[nMuonCand]/D");
    tree->Branch("MuonCand_eta", MuonCand_eta, "MuonCand_eta[nMuonCand]/D");
    tree->Branch("MuonCand_phi", MuonCand_phi, "MuonCand_phi[nMuonCand]/D");
    tree->Branch("MuonCand_charge", MuonCand_charge, "MuonCand_charge[nMuonCand]/I");
    tree->Branch("MuonCand_vtxx", MuonCand_vtxx, "MuonCand_vtxx[nMuonCand]/D");
    tree->Branch("MuonCand_vtxy", MuonCand_vtxy, "MuonCand_vtxy[nMuonCand]/D");
    tree->Branch("MuonCand_vtxz", MuonCand_vtxz, "MuonCand_vtxz[nMuonCand]/D");
    tree->Branch("MuonCand_dxy", MuonCand_dxy, "MuonCand_dxy[nMuonCand]/D");
    tree->Branch("MuonCand_dz", MuonCand_dz, "MuonCand_dz[nMuonCand]/D");
    tree->Branch("MuonCand_nstatseg", MuonCand_nstatseg, "MuonCand_nstatseg[nMuonCand]/I");
    tree->Branch("MuonCand_ntrklayers", MuonCand_ntrklayers, "MuonCand_ntrklayers[nMuonCand]/I");
    tree->Branch("MuonCand_npxlhits", MuonCand_npxlhits, "MuonCand_npxlhits[nMuonCand]/I");
    tree->Branch("MuonCand_isglobal", MuonCand_isglobal, "MuonCand_isglobal[nMuonCand]/I");
    tree->Branch("MuonCand_istracker", MuonCand_istracker, "MuonCand_istracker[nMuonCand]/I");
    tree->Branch("MuonCand_isstandalone", MuonCand_isstandalone, "MuonCand_isstandalone[nMuonCand]/I");
    tree->Branch("MuonCand_ispfmuon", MuonCand_ispfmuon, "MuonCand_ispfmuon[nMuonCand]/I");
    tree->Branch("MuonCand_istight", MuonCand_istight, "MuonCand_istight[nMuonCand]/I");
    tree->Branch("MuonCandTrack_nmuchits", MuonCandTrack_nmuchits, "MuonCandTrack_nmuchits[nMuonCand]/I");
    tree->Branch("MuonCandTrack_chisq", MuonCandTrack_chisq, "MuonCandTrack_chisq[nMuonCand]/D");
    tree->Branch("MuonCand_innerTrackPt", MuonCand_innerTrackPt, "MuonCand_innerTrackPt[nMuonCand]/D");
    tree->Branch("MuonCand_innerTrackEta", MuonCand_innerTrackEta, "MuonCand_innerTrackEta[nMuonCand]/D"); 
    tree->Branch("MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, "MuonCand_innerTrackPhi[nMuonCand]/D"); 

    if (runOnMC_) {
      tree->Branch("nGenMuonCand", &nGenMuonCand, "nGenMuonCand/I");
      tree->Branch("nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/I");    
      tree->Branch("GenMuonCand_p", GenMuonCand_p, "GenMuonCand_p[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_px", GenMuonCand_px, "GenMuonCand_px[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_py", GenMuonCand_py, "GenMuonCand_py[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_pz", GenMuonCand_pz, "GenMuonCand_pz[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D");
      tree->Branch("GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D");

      tree->Branch("nGenProtCand", &nGenProtCand, "nGenProtCand/I");     
      tree->Branch("GenProtCand_p", GenProtCand_p, "GenProtCand_p[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_px", GenProtCand_px, "GenProtCand_px[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_py", GenProtCand_py, "GenProtCand_py[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_pz", GenProtCand_pz, "GenProtCand_pz[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_e", GenProtCand_e, "GenProtCand_e[nGenProtCand]/D");  
      tree->Branch("GenProtCand_pt", GenProtCand_pt, "GenProtCand_pt[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_eta", GenProtCand_eta, "GenProtCand_eta[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_phi", GenProtCand_phi, "GenProtCand_phi[nGenProtCand]/D"); 
      tree->Branch("GenProtCand_status", GenProtCand_status, "GenProtCand_status[nGenProtCand]/I"); 
      tree->Branch("GenProtCand_xi", GenProtCand_xi, "GenProtCand_xi[nGenProtCand]/D");   
      tree->Branch("GenProtCand_t", GenProtCand_t, "GenProtCand_t[nGenProtCand]/D");   

    }
  }
  
  if(runOnProtons_){
/* 2016
    tree->Branch("nLocalProtCand", &nLocalProtCand, "nLocalProtCand/I");
    tree->Branch("LocalProtCand_x", LocalProtCand_x, "LocalProtCand_x[nLocalProtCand]/D");
    tree->Branch("LocalProtCand_y", LocalProtCand_y, "LocalProtCand_y[nLocalProtCand]/D"); 
    tree->Branch("LocalProtCand_z", LocalProtCand_z, "LocalProtCand_z[nLocalProtCand]/D"); 
    tree->Branch("LocalProtCand_xSigma", LocalProtCand_xSigma, "LocalProtCand_xSigma[nLocalProtCand]/D");
    tree->Branch("LocalProtCand_ySigma", LocalProtCand_ySigma, "LocalProtCand_ySigma[nLocalProtCand]/D"); 
    tree->Branch("LocalProtCand_Tx", LocalProtCand_Tx, "LocalProtCand_Tx[nLocalProtCand]/D");
    tree->Branch("LocalProtCand_Ty", LocalProtCand_Ty, "LocalProtCand_Ty[nLocalProtCand]/D"); 
    tree->Branch("LocalProtCand_TxSigma", LocalProtCand_TxSigma, "LocalProtCand_TxSigma[nLocalProtCand]/D"); 
    tree->Branch("LocalProtCand_TySigma", LocalProtCand_TySigma, "LocalProtCand_TySigma[nLocalProtCand]/D");  
*/

// PPS 2017
/*
  tree->Branch("nHitsTiming", &nHitsTiming, "nHitsTiming/I");
  tree->Branch("LeadingEdge", &LeadingEdge, "LeadingEdge[nHitsTiming]/D");
  tree->Branch("TrailingEdge", &TrailingEdge, "TrailingEdge[nHitsTiming]/D");
  tree->Branch("ToT", &ToT, "ToT[nHitsTiming]/D");
  tree->Branch("ChannelTiming", &ChannelTiming, "ChannelTiming[nHitsTiming]/I");
  tree->Branch("PlaneTiming", &PlaneTiming, "PlaneTiming[nHitsTiming]/I");
  tree->Branch("ArmTiming", &ArmTiming, "ArmTiming[nHitsTiming]/I");
  tree->Branch("XTiming", &XTiming, "XTiming[nHitsTiming]/D");
  tree->Branch("YTiming", &YTiming, "YTiming[nHitsTiming]/D");
*/
  tree->Branch("nRecHitsTiming", &nRecHitsTiming, "nRecHitsTiming/I");
  tree->Branch("TimingRecHitArm",&TimingRecHitArm,"TimingRecHitArm[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitChannel",&TimingRecHitChannel,"TimingRecHitChannel[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitPlane",&TimingRecHitPlane,"TimingRecHitPlane[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitT",&TimingRecHitT,"TimingRecHitT[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitX",&TimingRecHitX,"TimingRecHitX[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitY",&TimingRecHitY,"TimingRecHitY[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitToT",&TimingRecHitToT,"TimingRecHitToT[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitOOTIndex",&TimingRecHitOOTIndex,"TimingRecHitOOTIndex[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitMultiHit",&TimingRecHitMultiHit,"TimingRecHitMultiHit[nRecHitsTiming]/I");

  tree->Branch("nTracksTiming", &nTracksTiming, "nTracksTiming/I");
  tree->Branch("TimingTrackT", &TimingTrackT, "TimingTrackT[nTracksTiming]/D");
  tree->Branch("TimingTrackTErr", &TimingTrackTErr, "TimingTrackTErr[nTracksTiming]/D");
  tree->Branch("TimingTrackX", &TimingTrackX, "TimingTrackX[nTracksTiming]/D");
  tree->Branch("TimingTrackY", &TimingTrackY, "TimingTrackY[nTracksTiming]/D");
  tree->Branch("TimingTrackZ", &TimingTrackZ, "TimingTrackZ[nTracksTiming]/D");
  tree->Branch("TimingTrackOOTIndex", &TimingTrackOOTIndex, "TimingTrackOOTIndex[nTracksTiming]/I");
  tree->Branch("TimingTrackMultiHit", &TimingTrackMultiHit, "TimingTrackMultiHit[nTracksTiming]/I");
  tree->Branch("TimingTrackChi2", &TimingTrackChi2, "TimingTrackChi2[nTracksTiming]/D");
  tree->Branch("TimingTrackArm", &TimingTrackArm, "TimingTrackArm[nTracksTiming]/I");

  //  tree->Branch("IsClock", &IsClock, "IsClock[nHitsTiming]/I");
 
  tree->Branch("nArmsTimingRecHits", &nArmsTimingRecHits, "nArmsTimingRecHits/I");
  tree->Branch("nArmsTimingTracks", &nArmsTimingTracks, "nArmsTimingTracks/I");

 
  tree->Branch("nTracksStrips", &nTracksStrips, "nTracksStrips/I");
  tree->Branch("StripTrackX", &StripTrackX, "StripTrackX[nTracksStrips]/D");
  tree->Branch("StripTrackY", &StripTrackY, "StripTrackY[nTracksStrips]/D");
  tree->Branch("StripTrackZ", &StripTrackZ, "StripTrackZ[nTracksStrips]/D");

  tree->Branch("StripTrackTx", &StripTrackTx, "StripTrackTx[nTracksStrips]/D");
  tree->Branch("StripTrackTy", &StripTrackTy, "StripTrackTy[nTracksStrips]/D");
  tree->Branch("ArmStrips", &ArmStrips, "ArmStrips[nTracksStrips]/I");
  tree->Branch("StationStrips", &StationStrips, "StationStrips[nTracksStrips]/I");
  tree->Branch("RPStrips", &RPStrips, "RPStrips[nTracksStrips]/I");
  tree->Branch("StripsMysteryFlag", &StripsMysteryFlag, "StripsMysteryFlag[nTracksStrips]/I");

  tree->Branch("nArmsStrips", &nArmsStrips, "nArmsStrips/I");

  //tree->Branch("nArmsPixelDigis", &nArmsPixelDigis, "nArmsPixelDigis/I");
  //tree->Branch("nLayersArm1PixelDigis", &nLayersArm1PixelDigis, "nLayersArm1PixelDigis/I");
  //tree->Branch("nLayersArm2PixelDigis", &nLayersArm2PixelDigis, "nLayersArm2PixelDigis/I");
  tree->Branch("nPixelRecHits", &nPixelRecHits, "nPixelRecHits/I");
  tree->Branch("PixRecHitX", &PixRecHitX, "PixRecHitX[nPixelRecHits]/D");
  tree->Branch("PixRecHitY", &PixRecHitY, "PixRecHitY[nPixelRecHits]/D");
  tree->Branch("PixRecHitZ", &PixRecHitZ, "PixRecHitZ[nPixelRecHits]/D");
  tree->Branch("PixRecHitArm", &PixRecHitArm, "PixRecHitArm[nPixelRecHits]/I");
  tree->Branch("PixRecHitPlane", &PixRecHitPlane, "PixRecHitPlane[nPixelRecHits]/I");
  tree->Branch("nArmsPixRecHits", &nArmsPixRecHits, "nArmsPixRecHits/I");
  tree->Branch("nLayersArm1PixRecHits", &nLayersArm1PixRecHits, "nLayersArm1PixRecHits/I");
  tree->Branch("nLayersArm2PixRecHits", &nLayersArm2PixRecHits, "nLayersArm2PixRecHits/I");

  tree->Branch("nPixelTracks", &nPixelTracks, "nPixelTracks/I");
  tree->Branch("PixTrackX", &PixTrackX, "PixTrackX[nPixelTracks]/D");
  tree->Branch("PixTrackY", &PixTrackY, "PixTrackY[nPixelTracks]/D");
  tree->Branch("PixTrackTx", &PixTrackTx, "PixTrackTx[nPixelTracks]/D");
  tree->Branch("PixTrackTy", &PixTrackTy, "PixTrackTy[nPixelTracks]/D");
  tree->Branch("PixTrackChi2", &PixTrackChi2, "PixTrackChi2[nPixelTracks]/D");
  tree->Branch("PixTrackZ", &PixTrackZ, "PixTrackZ[nPixelTracks]/D");
  tree->Branch("PixTrackArm", &PixTrackArm, "PixTrackArm[nPixelTracks]/I");
  tree->Branch("nArmsPixelTracks", &nArmsPixelTracks, "nArmsPixelTracks/I");


  }

	// Primary vertices' information
  tree->Branch("nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I");
  tree->Branch("PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_matchedtracks", PrimVertexCand_matchedtracks, "PrimVertexCand_matchedtracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_unmatchedtracks", PrimVertexCand_unmatchedtracks, "PrimVertexCand_unmatchedtracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/I");
  tree->Branch("nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/I");

  // Kalman dilepton vertex information
  tree->Branch("KalmanVertexCand_x", &KalmanVertexCand_x, "KalmanVertexCand_x/D");
  tree->Branch("KalmanVertexCand_y", &KalmanVertexCand_y, "KalmanVertexCand_y/D");
  tree->Branch("KalmanVertexCand_z", &KalmanVertexCand_z, "KalmanVertexCand_z/D");
  tree->Branch("KalmanVertexCand_chi2", &KalmanVertexCand_chi2, "KalmanVertexCand_chi2/D");
  tree->Branch("KalmanVertexCand_ndof", &KalmanVertexCand_ndof, "KalmanVertexCand_ndof/I"); 
  tree->Branch("ClosestExtraTrackKalman_vtxdxyz", &ClosestExtraTrackKalman_vtxdxyz, "ClosestExtraTrackKalman_vtxdxyz/D");

  // Lepton pairs' information
  tree->Branch("nMuVertexCand", &nMuVertexCand, "nMuVertexCand/I");
  tree->Branch("Pair_candidates", Pair_candidates, "Pair_candidates[nMuVertexCand][2]/I");
  tree->Branch("Pair_mindist", Pair_mindist, "Pair_mindist[nMuVertexCand]/D");
  tree->Branch("Pair_p", Pair_p, "Pair_p[nMuVertexCand]/D");
  tree->Branch("Pair_pt", Pair_pt, "Pair_pt[nMuVertexCand]/D");
  tree->Branch("Pair_dpt", Pair_dpt, "Pair_dpt[nMuVertexCand]/D");
  tree->Branch("Pair_mass", Pair_mass, "Pair_mass[nMuVertexCand]/D");
  tree->Branch("Pair_eta", Pair_eta, "Pair_eta[nMuVertexCand]/D");
  tree->Branch("Pair_phi", Pair_phi, "Pair_phi[nMuVertexCand]/D");
  tree->Branch("Pair_dphi", Pair_dphi, "Pair_dphi[nMuVertexCand]/D");
  tree->Branch("Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nMuVertexCand]/D");
  tree->Branch("Pair_Y", Pair_Y, "Pair_Y[nMuVertexCand]/D"); 
  tree->Branch("Pair_extratracks0pt5mm",Pair_extratracks0pt5mm, "Pair_extratracks0pt5mm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks1mm", Pair_extratracks1mm, "Pair_extratracks1mm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks2mm", Pair_extratracks2mm, "Pair_extratracks2mm[nMuVertexCand]/I");

  tree->Branch("Pair_extratracks3mm", Pair_extratracks3mm, "Pair_extratracks3mm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks4mm", Pair_extratracks4mm, "Pair_extratracks4mm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks5mm", Pair_extratracks5mm, "Pair_extratracks5mm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks1cm", Pair_extratracks1cm, "Pair_extratracks1cm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks2cm", Pair_extratracks2cm, "Pair_extratracks2cm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks3cm", Pair_extratracks3cm, "Pair_extratracks3cm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks4cm", Pair_extratracks4cm, "Pair_extratracks4cm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks5cm", Pair_extratracks5cm, "Pair_extratracks5cm[nMuVertexCand]/I");
  tree->Branch("Pair_extratracks10cm", Pair_extratracks10cm, "Pair_extratracks10cm[nMuVertexCand]/I");
  if (runOnMC_) {
    tree->Branch("GenPair_p", &GenPair_p, "GenPair_p/D");
    tree->Branch("GenPair_pt", &GenPair_pt, "GenPair_pt/D");
    tree->Branch("GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D");
    tree->Branch("GenPair_mass", &GenPair_mass, "GenPair_mass/D");
    tree->Branch("GenPair_eta", &GenPair_eta, "GenPair_eta/D");
    tree->Branch("GenPair_Y", &GenPair_Y, "GenPair_Y/D");
    tree->Branch("GenPair_phi", &GenPair_phi, "GenPair_phi/D");
    tree->Branch("GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D");
    tree->Branch("GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle/D");
  }
  // Extra tracks on vertex's information
  tree->Branch("nExtraTracks", &nExtraTracks, "nExtraTracks/I");
  tree->Branch("ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I");
  tree->Branch("ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I");
  tree->Branch("ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I");
  tree->Branch("ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I");
  tree->Branch("ExtraTrack_p", ExtraTrack_p, "ExtraTrack_p[nExtraTracks]/D");
  tree->Branch("ExtraTrack_pt", ExtraTrack_pt, "ExtraTrack_pt[nExtraTracks]/D");
  tree->Branch("ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D");
  tree->Branch("ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D");
  tree->Branch("ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D");
  tree->Branch("ExtraTrack_eta", ExtraTrack_eta, "ExtraTrack_eta[nExtraTracks]/D");
  tree->Branch("ExtraTrack_phi", ExtraTrack_phi, "ExtraTrack_phi[nExtraTracks]/D");
  tree->Branch("ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D");
  tree->Branch("ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D");
  tree->Branch("ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D");
  tree->Branch("ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D");
  tree->Branch("nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/I");
  tree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D");
  tree->Branch("ClosestHighPurityExtraTrack_vtxdxyz",&ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz/D");
  tree->Branch("ClosestTrack_vtxdxyz",&ClosestTrack_vtxdxyz,"ClosestTrack_vtxdxyz/D");

  // Pileup reweighting
  tree->Branch("nTruePUforPUWeight",&nTruePUforPUWeight,"nTruePUforPUWeight/I");
  tree->Branch("nTruePUafterPUWeight",&nTruePUafterPUWeight,"nTruePUafterPUWeight/D");
  tree->Branch("nTruePUforPUWeightBXM1", &nTruePUforPUWeightBXM1, "nTruePUforPUWeightBXM1/I");
  tree->Branch("nTruePUafterPUWeightBXM1", &nTruePUafterPUWeightBXM1, "nTruePUafterPUWeightBXM1/D");
  tree->Branch("nTruePUforPUWeightBXP1", &nTruePUforPUWeightBXP1, "nTruePUforPUWeightBXP1/I"); 
  tree->Branch("nTruePUafterPUWeightBXP1", &nTruePUafterPUWeightBXP1, "nTruePUafterPUWeightBXP1/D"); 
  tree->Branch("nTruePUforPUWeightBX0", &nTruePUforPUWeightBX0, "nTruePUforPUWeightBX0/I");
  tree->Branch("nTruePUafterPUWeightBX0", &nTruePUafterPUWeightBX0, "nTruePUafterPUWeightBX0/D");
  tree->Branch("Weight", &Weight, "Weight/D");
  tree->Branch("PUWeightTrue", &PUWeightTrue, "PUWeightTrue/D");

  nCandidates = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMu::endJob() 
{
  //  std::cout << "==> Number of candidates in the dataset : " << nCandidates << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
GammaGammaMuMu::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed;
  std::string triggerName_;
  
  changed = true;
  if (hltConfig_.init(iRun, iSetup, hltMenuLabel_, changed)) {
    if (changed) {
      // check if trigger name in (new) config
     // triggerName_ = "HLT_DoubleMu0";
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
        const UInt_t n = hltConfig_.size();
        const UInt_t triggerIndex = hltConfig_.triggerIndex(triggerName_);
        if (triggerIndex>=n) {
          std::cout << "GammaGammaMuMu::analyze:"
                    << " TriggerName " << triggerName_ 
                    << " not available in (new) config!"
		    << std::endl;
          //std::cout << "Available TriggerNames are: " << std::endl;
          //hltConfig_.dump("Triggers");
        }
      }
      //hltConfig_.dump("Streams");
      //hltConfig_.dump("Datasets");
      //hltConfig_.dump("PrescaleTable");
      //hltConfig_.dump("ProcessPSet");
    }
  }
  else {
    std::cout << "GammaGammaMuMu::beginRun:"
              << " config extraction failure with process name "
              << hltMenuLabel_
	      << std::endl;
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
GammaGammaMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GammaGammaMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GammaGammaMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GammaGammaMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/// Class PrimaryVertex

//
// constructors and destructor
//
PrimaryVertex::PrimaryVertex(std::vector<std::string>& _leptonsType, std::map<Int_t,TLorentzVector>& _muonsMomenta, std::map<Int_t,TLorentzVector>& _electronsMomenta) :
  nTracks(0),
  nMatchedTracks(0),
  nUnmatchedTracks(0),
  nMatchedMuons(0),
  nMatchedElectrons(0),
  FetchMuons(false),
  FetchElectrons(false)
{
  //LeptonsType = _leptonsType;
  for (i=0; i<_leptonsType.size(); i++) {
    if (_leptonsType[i]=="Muon") {
      FetchMuons = true;
    }
    else if (_leptonsType[i]=="Electron") {
      FetchElectrons = true;
    }
  }
  MuonMomenta = _muonsMomenta;
  ElectronMomenta = _electronsMomenta;
}

PrimaryVertex::~PrimaryVertex() {
}

void
PrimaryVertex::SetPositionMuMu(Double_t _x, Double_t _y, Double_t _z)
{
  Position.SetXYZ(_x, _y, _z);
#ifdef DEBUG
  std::cout << "[PrimaryVertex] SetPosition :: Vertex located at (" << Position.x() << ", " << Position.y() << ", " << Position.z() << ")" << std::endl;
#endif
}

/**
 * \brief Matches a track arising from the vertex with a lepton track from the
 *  internal collection
 */
Int_t
PrimaryVertex::AddTrackMuMu(const reco::TrackRef& _track, TString& _leptonType)
{
  nTracks++; // total number of tracks matched with the vertex
  std::map<Int_t,TLorentzVector>::iterator lep;
  for (lep=MuonMomenta.begin(); lep!=MuonMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "muon";
    MatchedMuons.push_back(lep->first);
    nMatchedMuons++;
    nMatchedTracks++;
    return lep->first;
  }

  nUnmatchedTracks++;
  return -1;
}

Double_t
PrimaryVertex::dZ(TVector3 _vmu, Int_t _muind)
{
  TLorentzVector m(MuonMomenta[_muind]);
  return (_vmu.Z()-Position.Z())-((_vmu.X()-Position.X())*m.Px()+(_vmu.Y()-Position.Y())*m.Py())/m.Pt()*m.Pz()/m.Pt();
}

HLTmatches::HLTmatches(std::vector<std::string> _HLTlist)
{
  for (i=0; i<_HLTlist.size(); i++) {
    HLTnames.push_back(_HLTlist[i].substr(0, _HLTlist[i].find_first_of("*")));
  }
#ifdef DEBUG
  for (i=0; i<HLTnames.size(); i++) {
    std::cout << i << " ==> " << HLTnames[i] << std::endl;
  }
#endif
}

HLTmatches::~HLTmatches()
{
}

Int_t
HLTmatches::TriggerNum(std::string _trigName)
{
  for (i=0; i<HLTnames.size(); i++) {
    if (_trigName.find(HLTnames[i])!=std::string::npos) return i;
  }
  return -1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaMuMu);
