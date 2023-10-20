// -*- C++ -*-
//   
//         Author: Sayan Chatterjee
//         Class: SayanCMW 
//         Checked version 24/06/2020
//         Last Modified on 19/10/2023
//         miniAOD:: Data
//         

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "Analyzers/SayanCMW/interface/SayanCMW.h"
#include <string>

//  constructors and destructor

SayanCMW::SayanCMW(const edm::ParameterSet& iConfig) :  //Parametrized Constructor
  //******TRACKED PARAMETER********
  
  //tracks & vertex
  trackTags_(consumes< edm::View< pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTagsgen_(consumes< edm::View< pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("tracksgen"))),
  chi2Map_( consumes< edm::ValueMap< float > >( iConfig.getParameter< edm::InputTag >( "trackschi2" ) ) ),
  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),

  //centrality bin
  cent_bin_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralitybin"))),
  cent_reco_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentReco"))),
  binTable(iConfig.getUntrackedParameter< std::vector< double > >("binTable")),
  binTable_up(iConfig.getUntrackedParameter< std::vector< double > >("binTable_up")),
  binTable_down(iConfig.getUntrackedParameter< std::vector< double > >("binTable_down")),
  //******UNTRACKED PARAMETER****************
  
  //Event classifier
  centmin_(iConfig.getUntrackedParameter<int>("centmin")),
  centmax_(iConfig.getUntrackedParameter<int>("centmax")),
  
  //vertex selection
  zminVtx_(iConfig.getUntrackedParameter<double>("zminVtx")),
  zmaxVtx_(iConfig.getUntrackedParameter<double>("zmaxVtx")),

  //EFF Correction
  fpt_(iConfig.getUntrackedParameter<edm::InputTag>("fpt")),
  fmb_(iConfig.getUntrackedParameter<edm::InputTag>("fmb")),
  fplus_(iConfig.getUntrackedParameter<edm::InputTag>("fplus")),
  fminus_(iConfig.getUntrackedParameter<edm::InputTag>("fminus")),
  fpix_(iConfig.getUntrackedParameter<edm::InputTag>("fpix")),

  //track selection
  pTmin_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk")),
  pTmax_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk")),
  pTmin_trg_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk_trg")),
  pTmax_trg_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk_trg")),
  pTmin_ass_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk_ass")),
  pTmax_ass_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk_ass")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),
  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),

  //trigger                                                                                                                                                                         
  pt_trigg_min(iConfig.getUntrackedParameter<double>("ptmin_trigg")),
  pt_trigg_max(iConfig.getUntrackedParameter<double>("ptmax_trigg")),
  eta_trigg_min(iConfig.getUntrackedParameter<double>("etamin_trigg")),
  eta_trigg_max(iConfig.getUntrackedParameter<double>("etamax_trigg")),
  //associate                                                    
  pt_ass_min(iConfig.getUntrackedParameter<double>("ptmin_ass")),
  pt_ass_max(iConfig.getUntrackedParameter<double>("ptmax_ass")),
  eta_ass_min(iConfig.getUntrackedParameter<double>("etamin_ass")),
  eta_ass_max(iConfig.getUntrackedParameter<double>("etamax_ass")),

  //binning
  trg_ptbin(iConfig.getUntrackedParameter< std::vector < double > >("trigger_ptbin")),
  ass_ptbin(iConfig.getUntrackedParameter< std::vector < double > >("associate_ptbin")),
  etabining(iConfig.getUntrackedParameter< std::vector < double > >("variable_etabin")),
//etabining5(iConfig.getUntrackedParameter< std::vector < double > >("etabinwidth_p5")),
  etabining4(iConfig.getUntrackedParameter< std::vector < double > >("etabinwidth_p4")),
  etabining3(iConfig.getUntrackedParameter< std::vector < double > >("etabinwidth_p3")),
  etabining2(iConfig.getUntrackedParameter< std::vector < double > >("etabinwidth_p2")),
  pTbining(iConfig.getUntrackedParameter< std::vector < double > >("variable_pTbin")),
  phibining(iConfig.getUntrackedParameter< std::vector < double > >("variable_phibin")),
  centbining(iConfig.getUntrackedParameter< std::vector < double > >("variable_centbin")),
  zvtxbining3(iConfig.getUntrackedParameter< std::vector < double > >("zvtxbinwidth_3")),
  zvtxbining2p5(iConfig.getUntrackedParameter< std::vector < double > >("zvtxbinwidth_2p5")),
  zvtxbining2(iConfig.getUntrackedParameter< std::vector < double > >("zvtxbinwidth_2")),
  zvtxbining1p5(iConfig.getUntrackedParameter< std::vector < double > >("zvtxbinwidth_1p5")),
  zvtxbining1(iConfig.getUntrackedParameter< std::vector < double > >("zvtxbinwidth_1")),
  bkgFactor(iConfig.getUntrackedParameter<unsigned int>("bkgFactor")),
  ifMcreco_(iConfig.getUntrackedParameter<bool>("ifMcreco")),
  cent_nom(iConfig.getUntrackedParameter<bool>("cent_nom")),
  cent_up(iConfig.getUntrackedParameter<bool>("cent_up")),
  cent_down(iConfig.getUntrackedParameter<bool>("cent_down")),
  cent_2023_run3(iConfig.getUntrackedParameter<bool>("cent_run3"))
{
  //---Dihedron header file----                                                          
  evt_ = new DiHadronCorrelationEvt(pTmin_.size(), pTmin_trg_.size(), pTmin_ass_.size());

  //*****Defining Histograms & Profile Histograms********************                  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();
  TProfile2D::SetDefaultSumw2();

  //**************************For efficiency correction ****************************************************** 
  TString f_PT(fpt_.label().c_str());
  edm::FileInPath f1(Form("Analyzers/SayanCMW/data/EFF/effHpT/%s",f_PT.Data()));
 
  TString f_MB(fmb_.label().c_str());
  edm::FileInPath f2(Form("Analyzers/SayanCMW/data/EFF/effMB/%s",f_MB.Data()));

  TString f_Plus(fplus_.label().c_str());
  edm::FileInPath f3(Form("Analyzers/SayanCMW/data/EFF/plus/%s",f_Plus.Data()));

  TString f_Minus(fminus_.label().c_str());
  edm::FileInPath f4(Form("Analyzers/SayanCMW/data/EFF/minus/%s",f_Minus.Data()));

  TString f_Pix(fpix_.label().c_str());
  edm::FileInPath f5(Form("Analyzers/SayanCMW/data/EFF/effPix/%s",f_Pix.Data()));

   
  TrkEff = new TrkEff2018PbPb(  "general", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff1 = new TrkEff2018PbPb(  "generalMB+", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff2 = new TrkEff2018PbPb(  "generalMB-", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());   

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  int ptbin_all = (ptmax_ - ptmin_)/0.1;
  int ptbin_trg = (pt_trigg_max - pt_trigg_min)/0.1;
  int ptbin_ass = (pt_ass_max - pt_ass_min)/0.1;

  TFileDirectory fGlobalHist = fs->mkdir("QAplots");
  hpt          = fGlobalHist.make<TH1F>("hpt", "", ptbin_all, ptmin_, ptmax_ );
  hptP         = fGlobalHist.make<TH1F>("hptP", "", ptbin_all, ptmin_, ptmax_ );
  hptN         = fGlobalHist.make<TH1F>("hptN", "", ptbin_all, ptmin_, ptmax_ );
  heta         = fGlobalHist.make<TH1F>("heta", "", 48, -2.4, 2.4 );
  hetaP        = fGlobalHist.make<TH1F>("hetaP", "", 48, -2.4, 2.4 );
  hetaN        = fGlobalHist.make<TH1F>("hetaN", "", 48, -2.4, 2.4 );
  h_ptbin      = fGlobalHist.make<TH1F>("h_ptbin", "", (pTbining.size() -1), &pTbining[0] );
  heta_nbin    = fGlobalHist.make<TH1F>("heta_nbin", "", (etabining.size() -1), &etabining[0] );
  hetaP_nbin   = fGlobalHist.make<TH1F>("hetaP_nbin", "", (etabining.size() -1), &etabining[0] );
  hetaN_nbin   = fGlobalHist.make<TH1F>("hetaN_nbin", "", (etabining.size() -1), &etabining[0] );

  hphi_nbin    = fGlobalHist.make<TH1F>("hphi_nbin", "", (phibining.size() -1), &phibining[0]);
  hphiP_nbin   = fGlobalHist.make<TH1F>("hphiP_nbin", "", (phibining.size() -1), &phibining[0]);
  hphiN_nbin   = fGlobalHist.make<TH1F>("hphiN_nbin", "", (phibining.size() -1), &phibining[0]);

  hpt_trigg          = fGlobalHist.make<TH1F>("hpt_trigg", "", ptbin_trg, pt_trigg_min, pt_trigg_max );
  heta_nbin_trigg    = fGlobalHist.make<TH1F>("heta_nbin_trigg", "", (etabining.size() -1), &etabining[0] );
  hphi_nbin_trigg    = fGlobalHist.make<TH1F>("hphi_nbin_trigg", "", (phibining.size() -1), &phibining[0]);

  hpt_ass          = fGlobalHist.make<TH1F>("hpt_ass", "", ptbin_ass, pt_ass_min, pt_ass_max );
  heta_nbin_ass    = fGlobalHist.make<TH1F>("heta_nbin_ass", "", (etabining.size() -1), &etabining[0] );
  hphi_nbin_ass    = fGlobalHist.make<TH1F>("hphi_nbin_ass", "", (phibining.size() -1), &phibining[0]);
  
  hpt_trigg_check          = fGlobalHist.make<TH1F>("hpt_trigg_check", "", ptbin_trg, pt_trigg_min, pt_trigg_max );
  hpt_ass_check          = fGlobalHist.make<TH1F>("hpt_ass_check", "", ptbin_trg, pt_ass_min, pt_ass_max );


  h_nHits      = fGlobalHist.make<TH1D>("h_nHits", "", 100, 0., 100.);
  h_pterr      = fGlobalHist.make<TH1D>("h_pterr", "", 1000, -5., 5.);
  h_ptreso     = fGlobalHist.make<TH1D>("h_ptreso", "", 1000, -5., 5.);
  h_chi2       = fGlobalHist.make<TH1D>("h_chi2", "", 2000, -10., 10.);
  h_DCAZ       = fGlobalHist.make<TH1D>("h_DCAZ", "", 200, -10., 10.);
  h_DCAXY       = fGlobalHist.make<TH1D>("h_DCAXY", "", 200, -10., 10.);

  hZBestVtx    = fGlobalHist.make<TH1F>("hZvtx", "", 600, -30.0, 30.0);
  hcent_bin    = fGlobalHist.make<TH1F>("hcent_bin", "", 200, 0.0, 200.0);
  hcentbin_array    = fGlobalHist.make<TH1F>("hcentbin_array", "", (centbining.size() -1), &centbining[0]);
  
  hzvtxbin3_array    = fGlobalHist.make<TH1F>("hzvtxbin3_array", "", (zvtxbining3.size() -1), &zvtxbining3[0]);
  hzvtxbin2p5_array    = fGlobalHist.make<TH1F>("hzvtxbin2p5_array", "", (zvtxbining2p5.size() -1), &zvtxbining2p5[0]);
  hzvtxbin2_array    = fGlobalHist.make<TH1F>("hzvtxbin2_array", "", (zvtxbining2.size() -1), &zvtxbining2[0]);
  hzvtxbin1p5_array    = fGlobalHist.make<TH1F>("hzvtxbin1p5_array", "", (zvtxbining1p5.size() -1), &zvtxbining1p5[0]);
  hzvtxbin1_array    = fGlobalHist.make<TH1F>("hzvtxbin1_array", "", (zvtxbining1.size() -1), &zvtxbining1[0]);
  
  heta_p5binwidth    = fGlobalHist.make<TH1F>("heta_p5binwidth", "", 1, -2.4, 2.4 );
  //heta_p5binwidth    = fGlobalHist.make<TH1F>("heta_p5binwidth", "", 10, -2.5, 2.5 );
  //heta_p5binwidth    = fGlobalHist.make<TH1F>("heta_p5binwidth", "", (etabining5.size() -1), &etabining5[0] );
  heta_p4binwidth    = fGlobalHist.make<TH1F>("heta_p4binwidth", "", (etabining4.size() -1), &etabining4[0] );
  heta_p3binwidth    = fGlobalHist.make<TH1F>("heta_p3binwidth", "", (etabining3.size() -1), &etabining3[0] );
  heta_p2binwidth    = fGlobalHist.make<TH1F>("heta_p2binwidth", "", (etabining2.size() -1), &etabining2[0] );

  hpt_trg    = fGlobalHist.make<TH1F>("hpt_trg", "", (trg_ptbin.size() -1), &trg_ptbin[0] );
  hpt_asso    = fGlobalHist.make<TH1F>("hpt_asso", "", (ass_ptbin.size() -1), &ass_ptbin[0] );

  //*************************************************************************************************************************************                                                                
  tp1d_mpteta_nbin   = fGlobalHist.make<TProfile>("tp1d_mpt_nbin", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);
  tp1d_mptetaP_nbin   = fGlobalHist.make<TProfile>("tp1d_mptP_nbin", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);
  tp1d_mptetaN_nbin   = fGlobalHist.make<TProfile>("tp1d_mptN_nbin", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);

  //*************************************************************************************************************************************                                                                
  TFileDirectory f_effw = fs->mkdir("Effeciency_weight");
  hpt_w         = f_effw.make<TH1F>("hpt_w", "", 25, 0.5, 3.0);
  hptP_w         = f_effw.make<TH1F>("hptP_w", "", 25, 0.5, 3.0);
  hptN_w         = f_effw.make<TH1F>("hptN_w", "", 25, 0.5, 3.0);
  heta_w        = f_effw.make<TH1F>("heta_w", "", 48, -2.4, 2.4 );
  hetaP_w        = f_effw.make<TH1F>("hetaP_w", "", 48, -2.4, 2.4 );
  hetaN_w        = f_effw.make<TH1F>("hetaN_w", "", 48, -2.4, 2.4 );
  heta_nbin_w    = f_effw.make<TH1F>("heta_nbin_w", "", (etabining.size() -1), &etabining[0] );
  hetaP_nbin_w    = f_effw.make<TH1F>("hetaP_nbin_w", "", (etabining.size() -1), &etabining[0] );
  hetaN_nbin_w    = f_effw.make<TH1F>("hetaN_nbin_w", "", (etabining.size() -1), &etabining[0] );

  //*************************************************************************************************************************************
  tp1d_mpteta_nbin_w   = f_effw.make<TProfile>("tp1d_mpt_nbin_w", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);
  tp1d_mptetaP_nbin_w   = f_effw.make<TProfile>("tp1d_mptP_nbin_w", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);
  tp1d_mptetaN_nbin_w   = f_effw.make<TProfile>("tp1d_mptN_nbin_w", "", (etabining.size() -1), &etabining[0], -1e10, 1e10);

  //default
  int nEtaBins_ =32;
  int nPhiBins_ =32;
  
  //int nEtaBins_ =64;
  //int nPhiBins_ =64;

  double etamax_trg_ = 2.4;             //default 2.4
  double etamin_trg_ = -2.4;            //default -2.4
  double etamax_ass_ = 2.4;             //default 2.4
  double etamin_ass_ = -2.4;            //default -2.4

  double etaW = (etamax_trg_ - etamin_ass_ - etamin_trg_ + etamax_ass_) / nEtaBins_;
  double phiW = 2.0*(TMath::Pi())/nPhiBins_;
  double minEta = etamin_trg_ - etamax_ass_ - etaW/2;
  double maxEta = etamax_trg_ - etamin_ass_ + etaW/2.;
  double minPhi = -(TMath::Pi() - phiW)/2.0;
  double maxPhi = (TMath::Pi() * 3.0 - phiW)/2.0;

  TFileDirectory f_dpt = fs->mkdir("Jet_study");
  
  //zvtxbin
  for(unsigned int i = 0; i< trg_ptbin.size() -1; i++)
    {
      //for(unsigned int j = 0; j< ass_ptbin.size()-1; j++)
      for(unsigned int j = 0; j<= i; j++)
	{
	  hsignal_c2_zvtx[i][j]   = f_dpt.make<TH2D>(Form("hsignal_c2_pttrg_%d_ptass_%d",i,j), "", nEtaBins_ + 1, minEta, maxEta, nPhiBins_ - 1, minPhi, maxPhi);
	  hsignal_c2_zvtx_2eff[i][j]   = f_dpt.make<TH2D>(Form("hsignal_c2_2eff_pttrg_%d_ptass_%d",i,j), "", nEtaBins_ + 1, minEta, maxEta, nPhiBins_ - 1, minPhi, maxPhi);
	  
	  //mixing
	  hsignal_c2_zvtx_mix[i][j]   = f_dpt.make<TH2D>(Form("hsignal_c2_mix_pttrg_%d_ptass_%d",i,j), "", nEtaBins_ + 1, minEta, maxEta, nPhiBins_ - 1, minPhi, maxPhi);
	  hsignal_c2_zvtx_mix_2eff[i][j]   = f_dpt.make<TH2D>(Form("hsignal_c2_mix_2eff_pttrg_%d_ptass_%d",i,j), "", nEtaBins_ + 1, minEta, maxEta, nPhiBins_ - 1, minPhi, maxPhi);
	  
	}
      
      hntrg_addbincontent[i]  = f_dpt.make<TH1D>(Form("hntrg_addbincontent_pttrg_%d",i), "", 3, 0, 3);
      hntrg_corr_addbincontent[i]  = f_dpt.make<TH1D>(Form("hntrg_corr_addbincontent_pttrg_%d",i), "", 3, 0, 3);
    }
}


SayanCMW::~SayanCMW() // Destructor 
{
  delete evt_;
}


//---------------method called for each event-------------------------------------------------------

void
SayanCMW::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  LoopCMWVertices(iEvent, iSetup);

}

//----------------method called once each job just before starting event loop---------------------------

void
SayanCMW::beginJob()
{
}

//-------------method called once each job just before ending the event loop--------------------------------------------

void
SayanCMW::endJob()
{

  std::cout<< "Start sorting the events!" << std::endl;
  std::sort(evtVec_.begin(),evtVec_.end());
  std::cout<< "Finish sorting the events!" << std::endl;

  std::cout<< "Total of " << evtVec_.size() << " events are selected! > endJob" << std::endl;

  //second event loop for analysis
  for( unsigned int ievt =0; ievt < evtVec_.size(); ievt++)
    {
      if( ievt % 100 ==0 ) std::cout<< " Processing " << ievt << "th event"<< "second event loop " <<std::endl;
      
      float cent_ = evtVec_[ievt].cent;
      float zvtx_ = evtVec_[ievt].zvtx;
      unsigned int itrack = 0;
      unsigned int ntrack = evtVec_[ievt].pVect[itrack].size();
      
      // ---- trigger track loop:: track loop 1
      for(unsigned int n = 0; n < ntrack; n++)
	{	  
	  TVector3 pvector_n = (evtVec_[ievt].pVect[itrack])[n];
	  int chg_n = (evtVec_[ievt].chgVect[itrack])[n];
	  double w_n = (evtVec_[ievt].weightVect[itrack])[n];
	  
	  double px_n = pvector_n.Px();
	  double py_n = pvector_n.Py();
	  double pz_n = pvector_n.Pz();
	  double pt_n = pvector_n.Pt();
	  double eta_n = pvector_n.Eta();
	  double phi_n = pvector_n.Phi();

	  double p_n = std::sqrt(px_n*px_n + py_n*py_n + pz_n*pz_n);
	  
	  //trigger pt cuts                                        
          if (pt_n < pt_trigg_min || pt_n >= pt_trigg_max) continue;
	  if (eta_n < eta_trigg_min || eta_n >= eta_trigg_max) continue;

	  int kpt_trg = (hpt_trg->FindBin(pt_n)) - 1;

	  // ---- associated track loop:: track loop 2
	  for(unsigned int m = 0; m < ntrack; m++)
	    {	  
	      TVector3 pvector_m = (evtVec_[ievt].pVect[itrack])[m];
	      int chg_m = (evtVec_[ievt].chgVect[itrack])[m];
	      double w_m = (evtVec_[ievt].weightVect[itrack])[m];
	      	      
	      double px_m = pvector_m.Px();
	      double py_m = pvector_m.Py();
	      double pz_m = pvector_m.Pz();
	      double pt_m = pvector_m.Pt();
	      double eta_m = pvector_m.Eta();
	      double phi_m = pvector_m.Phi();
	      
	      double p_m = std::sqrt(px_m*px_m + py_m*py_m + pz_m*pz_m);
	      
	      //associate pt cuts                                                                           
              if (pt_m < pt_ass_min || pt_m >= pt_ass_max) continue;
	      if (eta_m < eta_ass_min || eta_m >= eta_ass_max) continue;
              if (chg_n == chg_m  &&  phi_n == phi_m && eta_n == eta_m && pt_n == pt_m) continue;
	      
	      int kpt_ass = (hpt_asso->FindBin(pt_m)) - 1;
	      if (kpt_trg < kpt_ass) continue;
	      
	      //Track Splitting and Track Merging                                                          
              double cos_cut = 0.99999;
              double dpt_cut = 0.015;
	      Double_t cosa = TMath::Abs(px_m*px_n + py_m*py_n + pz_m*pz_n)/(p_m*p_n);
              Double_t deltapt = TMath::Abs(pt_m - pt_n);
              //if((cosa > cos_cut) && (deltapt < dpt_cut))continue;
	      	     
	      double deltaPhi = GetDeltaPhi(phi_n, phi_m);
	      double deltaPhi2 = GetDeltaPhi(phi_m, phi_n);
	      double deltaEta = GetDeltaEta(eta_n, eta_m);
	      
	      double w_nm = w_n*w_m;

	      //zvtxbin and ptbin
	      hsignal_c2_zvtx[kpt_trg][kpt_ass]->Fill( deltaEta, deltaPhi, 1.0 );
	      hsignal_c2_zvtx_2eff[kpt_trg][kpt_ass]->Fill( deltaEta, deltaPhi, 1.0*w_nm );

	      /*
		//zvtxbin and ptbin
		hsignal_c2_zvtx[kzvtx][keta5_n]->Fill( fabs(deltaEta), deltaPhi, 1.0/4.0 );
		hsignal_c2_zvtx[kzvtx][keta5_n]->Fill(-fabs(deltaEta), deltaPhi, 1.0/4.0 );
		hsignal_c2_zvtx[kzvtx][keta5_n]->Fill( fabs(deltaEta), deltaPhi2, 1.0/4.0 );
		hsignal_c2_zvtx[kzvtx][keta5_n]->Fill(-fabs(deltaEta), deltaPhi2, 1.0/4.0 );
	      
		hsignal_c2_zvtx_2eff[kzvtx][keta5_n]->Fill( fabs(deltaEta), deltaPhi, 1.0*w_nm/4.0 );
		hsignal_c2_zvtx_2eff[kzvtx][keta5_n]->Fill(-fabs(deltaEta), deltaPhi, 1.0*w_nm/4.0 );
		hsignal_c2_zvtx_2eff[kzvtx][keta5_n]->Fill( fabs(deltaEta), deltaPhi2, 1.0*w_nm/4.0 );
		hsignal_c2_zvtx_2eff[kzvtx][keta5_n]->Fill(-fabs(deltaEta), deltaPhi2, 1.0*w_nm/4.0 );
		
	      */
	    }//end of associated track loop
  
	}//end of trigger track loop



      //*********************************************************************************************************************                                                        
      //added by sayan:: mixing event background                                                                                       
      unsigned int mixstart = ievt - bkgFactor/2;
      unsigned int mixend   = ievt + bkgFactor/2 + 1;                                                                                             
      
      /*    
      if(ievt < bkgFactor)                                                                                                                         
	{                                                   
	  mixstart = 0;                                                                                                                              
	  mixend   = 2*bkgFactor + 1;                                                                                                              
	}
      
      else if(ievt > evtVec_.size() - bkgFactor - 1)                                     
	{                                                                                                                                            
	  mixstart = evtVec_.size() - 2*bkgFactor - 1;                                                                                              
	  mixend   = evtVec_.size();                                                                                                               
	}                                                                                                                                   

      if( mixend > evtVec_.size() )                                                                                                                
	mixend = evtVec_.size();       
      */      

      if(ievt < bkgFactor/2)                                                                                                                         
	{                                                   
	  mixstart = 0;                                                                                                                              
	  mixend   = bkgFactor + 1;                                                                                                              
	}
      
      else if(ievt > evtVec_.size() - bkgFactor/2 - 1)                                     
	{                                                                                                                                            
	  mixstart = evtVec_.size() - bkgFactor - 1;                                                                                              
	  mixend   = evtVec_.size();                                                                                                               
	}                                                                                                                                   

      if( mixend > evtVec_.size() )                                                                                                                
	mixend = evtVec_.size();       
      




      for( unsigned int jevt = mixstart; jevt < mixend; jevt++ )
	{
	  if(ievt == jevt) continue;

	  if( evtVec_[ievt].run ==  evtVec_[jevt].run && evtVec_[ievt].event == evtVec_[jevt].event ) continue;

	  double deltazvtx = evtVec_[ievt].zvtx-evtVec_[jevt].zvtx;
	  if(fabs(deltazvtx) > 2.0) continue;
	    
	  unsigned int nsize_ievt = evtVec_[ievt].pVect[itrack].size();
	  unsigned int nsize_jevt = evtVec_[jevt].pVect[itrack].size();
	  
	  if (nsize_ievt ==0 || nsize_jevt ==0) continue;
	  
	  //differential track loop for background, ******** event: ievt *******************                 
                                                                        
	  for(unsigned int mi = 0; mi < nsize_ievt; mi++)
	    {
	      TVector3 pvector_mi = (evtVec_[ievt].pVect[itrack])[mi];
	      int chg_mi = (evtVec_[ievt].chgVect[itrack])[mi];
	      double w_mi = (evtVec_[ievt].weightVect[itrack])[mi];
	      
	      double px_mi = pvector_mi.Px();
	      double py_mi = pvector_mi.Py();
	      double pz_mi = pvector_mi.Pz();
	      double pt_mi = pvector_mi.Pt();
	      double eta_mi = pvector_mi.Eta();
	      double phi_mi = pvector_mi.Phi();
	      
	      double p_mi = std::sqrt(px_mi*px_mi + py_mi*py_mi + pz_mi*pz_mi);

	      //trigger pt cuts                                           
	      if (pt_mi < pt_trigg_min || pt_mi >= pt_trigg_max) continue;
	      if (eta_mi < eta_trigg_min || eta_mi >= eta_trigg_max) continue;
		
	      int kpt_trg = (hpt_trg->FindBin(pt_mi)) - 1;

	      //differential track loop for background, ******** event: jevt ******************* 
	      for(unsigned int mj = 0; mj < nsize_jevt; mj++)
		{
		  TVector3 pvector_mj = (evtVec_[jevt].pVect[itrack])[mj];
		  int chg_mj = (evtVec_[jevt].chgVect[itrack])[mj];
		  double w_mj = (evtVec_[jevt].weightVect[itrack])[mj];
		  
		  double px_mj = pvector_mj.Px();
		  double py_mj = pvector_mj.Py();
		  double pz_mj = pvector_mj.Pz();
		  double pt_mj = pvector_mj.Pt();
		  double eta_mj = pvector_mj.Eta();
		  double phi_mj = pvector_mj.Phi();
	
		  double p_mj = std::sqrt(px_mj*px_mj + py_mj*py_mj + pz_mj*pz_mj);
	  
		  //associate pt cuts                                     
		  if (pt_mj < pt_ass_min || pt_mj >= pt_ass_max) continue;
		  if (eta_mj < eta_ass_min || eta_mj >= eta_ass_max) continue;
		  if (pt_mi == pt_mj && chg_mi == chg_mj && eta_mi == eta_mj && phi_mi == phi_mj) continue;

		  int kpt_ass = (hpt_asso->FindBin(pt_mj)) - 1;
		  if (kpt_trg < kpt_ass) continue;

		  //Track Splitting and Track Merging                                                          
		  double cos_cut_mix = 0.99999;
		  double dpt_cut_mix = 0.015;
		  Double_t cosa_mix = TMath::Abs(px_mj*px_mi + py_mj*py_mi + pz_mj*pz_mi)/(p_mj*p_mi);
		  Double_t deltapt_mix = TMath::Abs(pt_mj - pt_mi);
		  //if((cosa_mix > cos_cut_mix) && (deltapt_mix < dpt_cut_mix))continue;
		  
		  double deltaPhi = GetDeltaPhi( phi_mi, phi_mj );
		  double deltaPhi2 = GetDeltaPhi( phi_mj, phi_mi );
		  double deltaEta = GetDeltaEta( eta_mi, eta_mj );
		  
		  double w_mix = w_mi*w_mj;
		  

		  //zvtx
		  hsignal_c2_zvtx_mix[kpt_trg][kpt_ass]->Fill( deltaEta, deltaPhi, 1.0 );
		  hsignal_c2_zvtx_mix_2eff[kpt_trg][kpt_ass]->Fill( deltaEta, deltaPhi, 1.0*w_mix );
		  
		  /*		  
		  //zvtx
		  hsignal_c2_zvtx_mix[kzvtx1][keta5_m]->Fill( fabs(deltaEta), deltaPhi, 1.0/4.0 );
		  hsignal_c2_zvtx_mix[kzvtx1][keta5_m]->Fill(-fabs(deltaEta), deltaPhi, 1.0/4.0 );
		  hsignal_c2_zvtx_mix[kzvtx1][keta5_m]->Fill( fabs(deltaEta), deltaPhi2, 1.0/4.0 );
		  hsignal_c2_zvtx_mix[kzvtx1][keta5_m]->Fill(-fabs(deltaEta), deltaPhi2, 1.0/4.0 );
		  
		  hsignal_c2_zvtx_mix_2eff[kzvtx1][keta5_m]->Fill( fabs(deltaEta), deltaPhi, 1.0*w_mix/4.0 );
		  hsignal_c2_zvtx_mix_2eff[kzvtx1][keta5_m]->Fill(-fabs(deltaEta), deltaPhi, 1.0*w_mix/4.0 );
		  hsignal_c2_zvtx_mix_2eff[kzvtx1][keta5_m]->Fill( fabs(deltaEta), deltaPhi2, 1.0*w_mix/4.0 );
		  hsignal_c2_zvtx_mix_2eff[kzvtx1][keta5_m]->Fill(-fabs(deltaEta), deltaPhi2, 1.0*w_mix/4.0 );
		  */

		}// mj loop
		
	      }//mi loop
	    
	  }// end of mixing event loop

	//}// bool mixing
      
    } //end of the 2nd event loop
  
}

//==============================================================================================

void
SayanCMW::fillDescriptions(edm::ConfigurationDescriptions&  descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//=============================================================================

void
SayanCMW::LoopCMWVertices( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  //track collection
  auto trks = iEvent.getHandle( trackTags_ );

  auto trksgen = iEvent.getHandle( trackTagsgen_ );

  //access tracks chi2/ndf
  auto chi2Map = iEvent.getHandle( chi2Map_ ); 

  //vtx collection
  auto pvs = iEvent.getHandle( vtxTags_ );

  //best vertex
  double bestvzError;
  math::XYZPoint bestvtx;
  math::Error<3>::type vtx_cov;
  if ( !pvs->empty() ) {
    const reco::Vertex& vtx = (*pvs)[0];
    bestvzError = vtx.zError();
    bestvtx = vtx.position();
    vtx_cov = vtx.covariance();      
  }else { 
    return; 
  }

  xBestVtx_ = bestvtx.x();
  yBestVtx_ = bestvtx.y();
  zBestVtx_ = bestvtx.z();

  if ( zBestVtx_ < zminVtx_ || zBestVtx_ >= zmaxVtx_ ) return; 
  //if ( fabs(zBestVtx_) < zminVtx_ || fabs(zBestVtx_) > zmaxVtx_ ) return; 


  // ----------------- centrality selection -------------------------------

  //access centrality bins
  auto cbin = iEvent.getHandle( cent_bin_ );
  //centBin = ( float ) (*cbin);
  
  //access reco centrality
  auto recocent = iEvent.getHandle(cent_reco_);
  const Double_t hiHF = recocent->EtHFtowerSum();
  
  // nominal, up, down

  if (cent_2023_run3)
    {
      centBin = ( float ) (*cbin);
      //if (cent_nom) centBin = getHiBinFromhiHF(hiHF, true, false, false);
    }
  else {
  
  if (cent_nom) centBin = getHiBinFromhiHF(hiHF, true, false, false);
  if (cent_up) centBin = getHiBinFromhiHF(hiHF, false, true, false);
  if (cent_down) centBin = getHiBinFromhiHF(hiHF, false, false, true);
  }
  //centBin = getHiBinFromhiHF(hiHF);
  if (centBin < centmin_*2 || centBin >= centmax_*2) return;

  hZBestVtx -> Fill(zBestVtx_); 
  //hcent_bin -> Fill(centBin/2.);
  hcent_bin -> Fill(centBin);
  hcentbin_array -> Fill(centBin);
  
  run_no = iEvent.id().run();
  event_no = iEvent.id().event();

  evt_->run = run_no;
  evt_->event = event_no;
  evt_->zvtx = zBestVtx_;
  //evt_->cent = (centBin/2.);
  evt_->cent = (centBin);


  //********* start track loop *********
  
  if(ifMcreco_){
  int trkIndx = -1;
  // Loop over tracks
  for (auto const& trk : *trks)
    {
      trkIndx++;
      
      if ( !trk.hasTrackDetails() ) continue;
      auto iter_tk = trk.pseudoTrack();
      
      double dzvtx = iter_tk.dz( bestvtx );
      double dxyvtx = iter_tk.dxy( bestvtx );
      double dzerror = std::hypot( iter_tk.dzError(), bestvzError );
      double dxyerror = iter_tk.dxyError( bestvtx, vtx_cov );
      double pterror = iter_tk.ptError();
      
      // Get eta, pt, and charge of the track
      double pt = iter_tk.pt();
      double eta = iter_tk.eta();
      int charge = iter_tk.charge();
      double phi = iter_tk.phi();
      
      auto hit_pattern = iter_tk.hitPattern();
      
      //HI specific cuts
      double chi2ndof = ( double ) ( *chi2Map )[ trks->ptrAt( trkIndx ) ];
      double dcaxy = (dxyvtx / dxyerror);
      double dcaz = (dzvtx / dzerror);
      double ptreso = (fabs(pterror) / pt);
      int nHits = iter_tk.numberOfValidHits();
      double chi2n = ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() );

      //selected tracks
      if (cent_2023_run3)
	{
	  if( charge == 0 ) continue;
	}
      else
	{
	  if( charge == 0 ) continue;
	  if( fabs(pterror) / pt >= 0.1 ) continue;                                           //Default cuts
	  if( fabs(dzvtx / dzerror) >= 3.0 ) continue;                                        //Default cuts
	  if( fabs(dxyvtx / dxyerror) >= 3.0  ) continue;                                     //Default cuts
	  if( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.18 ) continue;   //Default cuts
	  
	  if ( iter_tk.numberOfValidHits() < 11 ) continue;
	}
      
      /*
	if( fabs(pterror) / pt >= 0.15 ) continue;                                           //Loose cuts
	if( fabs(dzvtx / dzerror) >= 5.0 ) continue;                                        //Loose cuts
	if( fabs(dxyvtx / dxyerror) >= 5.0  ) continue;                                     //Loose cuts
	if ( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.18 ) continue;  //Loose cuts
	
      
	if( fabs(pterror) / pt >= 0.05 ) continue;                                           //Tight cuts
	if( fabs(dzvtx / dzerror) >= 2.0 ) continue;                                        //Tight cuts
	if( fabs(dxyvtx / dxyerror) >= 2.0  ) continue;                                     //Tight cuts
	if ( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.15 ) continue;  //Tight cuts
      */
      
      if( pt <= ptmin_ || pt > ptmax_ ) continue;
      //if( pt <= ptmin_ ) continue;                                                        //ptmin = 0.5
      if( eta <= etamin_ || eta >= etamax_ ) continue;                                    //|eta|<2.4
      //if ( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.18 ) continue;  //Default cuts
      //if ( iter_tk.numberOfValidHits() < 11 ) continue;
      
      int index = GetpTbin(pt);
      if(index == -1) continue;
      //int index = 0;
      
      float wgt = TrkEff->getCorrection(pt, eta, centBin);
      //float wgt = TrkEff->getCorrection(pt, eta, zBestVtx_, centBin);
      
      float weight = -999.0;
      
      
      //weight(pt,eta,centbin)
      if (charge > 0)
	{
	  weight = TrkEff1->getCorrection(pt, eta, centBin);
	  if(weight == -999.0) continue;
	}
      else
	{
	  weight = TrkEff2->getCorrection(pt, eta, centBin);
	  if(weight == -999.0) continue;
	}
      
      /*
      //weight(pt,eta,zvtx,centbin)
      if (charge > 0)
	{
	  weight = TrkEff1->getCorrection(pt, eta, zBestVtx_, centBin);
	  if(weight == -999.0) continue;
	}
      else
	{
	  weight = TrkEff2->getCorrection(pt, eta, zBestVtx_, centBin);
	  if(weight == -999.0) continue;
	}
      */

      if (cent_2023_run3) weight = 1.0;
      AssignpTbins(pt, eta, phi, wgt, weight, charge, index);
      
      //=========Filling histograms======================
      
      hpt->Fill(pt);
      heta->Fill(eta);
      heta_nbin->Fill(eta);
      hphi_nbin->Fill(phi);
      h_nHits->Fill(nHits);
      h_pterr->Fill(fabs(pterror));
      h_ptreso->Fill(ptreso);
      h_chi2->Fill(chi2n);
      h_DCAZ->Fill(dcaz);
      h_DCAXY->Fill(dcaxy);
      
      h_ptbin->Fill(pt);
      
      hpt_w ->Fill(pt, weight);
      heta_w ->Fill(eta, weight);
      heta_nbin_w ->Fill(eta, weight);

      int kpt_trg = (hpt_trg->FindBin(pt)) - 1;
      
      if(pt >= pt_trigg_min && pt < pt_trigg_max && eta >= eta_trigg_min && eta < eta_trigg_max )
        {
	  hpt_trg->Fill(pt);
          hpt_trigg->Fill(pt);
          heta_nbin_trigg->Fill(eta);
          hphi_nbin_trigg->Fill(phi);
	  
	  hntrg_addbincontent[kpt_trg] ->AddBinContent(1, 1);
          hntrg_corr_addbincontent[kpt_trg] ->AddBinContent(1, weight);

	  if(charge > 0)
	    {
	      hntrg_addbincontent[kpt_trg] ->AddBinContent(2, 1);
              hntrg_corr_addbincontent[kpt_trg] ->AddBinContent(2, weight);
	    }
	  
	  if(charge < 0)
	    {
	      hntrg_addbincontent[kpt_trg] ->AddBinContent(3, 1);
              hntrg_corr_addbincontent[kpt_trg] ->AddBinContent(3, weight);
	    }
	  
	}

      if(pt >= pt_ass_min && pt < pt_ass_max && eta >= eta_ass_min && eta < eta_ass_max   )
        {
	  hpt_asso->Fill(pt);
          hpt_ass->Fill(pt);
          heta_nbin_ass->Fill(eta);
          hphi_nbin_ass->Fill(phi);
        }


      
      if (charge > 0)
	{
	  hptP ->Fill(pt);
	  hetaP ->Fill(eta);
	  hetaP_nbin ->Fill(eta);
	  hphiP_nbin->Fill(phi);

	  hptP_w ->Fill(pt, weight);
          hetaP_w ->Fill(eta, weight);
          hetaP_nbin_w ->Fill(eta, weight);
          
	  tp1d_mptetaP_nbin ->Fill(eta, pt);
	  tp1d_mptetaP_nbin_w ->Fill(eta, pt, weight);
	}
      
      if (charge < 0)
	{
	  hptN ->Fill(pt);
	  hetaN ->Fill(eta);
	  hetaN_nbin ->Fill(eta);
	  hphiN_nbin->Fill(phi);

	  hptN_w ->Fill(pt, weight);
          hetaN_w ->Fill(eta, weight);
          hetaN_nbin_w ->Fill(eta, weight);
          
	  tp1d_mptetaN_nbin ->Fill(eta, pt);
	  tp1d_mptetaN_nbin_w ->Fill(eta, pt, weight);
	}
    } //end of Track loop
  }
  else
    {
      for (auto const& iter_tk : *trksgen)
        {
          if(iter_tk.status() != 1) continue;

          // Get eta, pt, and charge of the track                                                                                                                     
          double pt = iter_tk.pt();
          double eta = iter_tk.eta();
          int charge = iter_tk.charge();
          double phi = iter_tk.phi();

          //selected tracks                                                                                                                                           
          if( charge == 0 ) continue;
          if( pt <= ptmin_ || pt > ptmax_ ) continue;
          if( eta <= etamin_ || eta >= etamax_ ) continue;
	  
	  
          int index = GetpTbin(pt);
          if(index == -1) continue;
	  
          float weight = 1.0, wgt =1.0;
	  
          AssignpTbins(pt, eta, phi, wgt, weight, charge, index);
	  
          //=========Filling histograms======================                                                                                               
          hpt->Fill(pt);
          heta->Fill(eta);
          heta_nbin->Fill(eta);
          hphi_nbin->Fill(phi);
          h_ptbin->Fill(pt);
	  
          hpt_w ->Fill(pt, weight);
          heta_w ->Fill(eta, weight);
          heta_nbin_w ->Fill(eta, weight);

          if(pt >= pt_trigg_min && pt < pt_trigg_max && eta >= eta_trigg_min && eta < eta_trigg_max )
            {
              hpt_trigg->Fill(pt);
              heta_nbin_trigg->Fill(eta);
              hphi_nbin_trigg->Fill(phi);
            }

          if(pt >= pt_ass_min && pt < pt_ass_max && eta >= eta_ass_min && eta < eta_ass_max   )
            {
              hpt_ass->Fill(pt);
              heta_nbin_ass->Fill(eta);
              hphi_nbin_ass->Fill(phi);
            }

          if (charge > 0)
            {
              hptP ->Fill(pt);
              hetaP ->Fill(eta);
              hetaP_nbin ->Fill(eta);
              hphiP_nbin->Fill(phi);

              hptP_w ->Fill(pt, weight);
              hetaP_w ->Fill(eta, weight);
              hetaP_nbin_w ->Fill(eta, weight);

              tp1d_mptetaP_nbin ->Fill(eta, pt);
              tp1d_mptetaP_nbin_w ->Fill(eta, pt, weight);
            }

          if (charge < 0)
            {
              hptN ->Fill(pt);
              hetaN ->Fill(eta);
              hetaN_nbin ->Fill(eta);
              hphiN_nbin->Fill(phi);

              hptN_w ->Fill(pt, weight);
              hetaN_w ->Fill(eta, weight);
              hetaN_nbin_w ->Fill(eta, weight);

              tp1d_mptetaN_nbin ->Fill(eta, pt);
              tp1d_mptetaN_nbin_w ->Fill(eta, pt, weight);
            }
        }//end of track loop                                                                                                                                          
    }// end of if loop                       

  evtVec_.push_back(*evt_);

  //reset evt container                                                                                                                                                                                   
  evt_->reset();
  
}//end of LoopCMWVertices


//=============================== Getptbin================================
int SayanCMW::GetpTbin(double pt)
{
  int idx = -1;
  
  for(unsigned int i = 0; i<pTmin_.size(); ++i)
    {
      if( pt >= pTmin_[i] && pt <= pTmax_[i]) idx = i;
    }
  
  return idx;
}

//=========================== AssignpTbins =================================
  
void
SayanCMW::AssignpTbins(double pt, double eta, double phi, float wgt, float weight, int charge, int idx)
{
  TVector3 pvector;
  pvector.SetPtEtaPhi(pt, eta, phi);
    
  (evt_->pVect[idx]).push_back(pvector);
  (evt_->chgVect[idx]).push_back(charge);
  (evt_->wgtallVect[idx]).push_back(wgt);
    
  if ( charge > 0)
    {
      (evt_->pVect_trg[idx]).push_back(pvector);
      (evt_->chgVect_trg[idx]).push_back(charge);
      (evt_->weightVect[idx]).push_back(weight);
    }
  else
    {
      (evt_->pVect_ass[idx]).push_back(pvector);
      (evt_->chgVect_ass[idx]).push_back(charge);
      (evt_->weightVect[idx]).push_back(weight);
    }
}

//=================================================================                                                                                                                  
double
SayanCMW::GetDeltaEta(double eta_trg, double eta_ass)
{
  double deltaEta = eta_ass - eta_trg;                                                                                                                                             
  //double deltaEta = eta_trg - eta_ass;
  return deltaEta;
}
//=================================================================                                                                                                                  
double
SayanCMW::GetDeltaPhi(double phi_trg, double phi_ass)
{
  //double deltaPhi = phi_trg - phi_ass;
  double deltaPhi = phi_ass - phi_trg;                                                                                                                                             

  if(deltaPhi > 1.5*TMath::Pi())
    deltaPhi = deltaPhi - 2.0*TMath::Pi();

  else if(deltaPhi < -1.0*TMath::Pi() / 2.0)
    deltaPhi = deltaPhi + 2.0*TMath::Pi();

  return deltaPhi;

}

Int_t 
//SayanCMW::getHiBinFromhiHF(const Double_t hiHF, Double_t *binTable)
SayanCMW::getHiBinFromhiHF(const Double_t hiHF, bool nom, bool up, bool down)
{
  Int_t binPos = -1;
  for(int i = 0; i < ncBins; ++i){
    if(nom){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
    }
    }
    else if (up) {
      if(hiHF >= binTable_up[i] && hiHF < binTable_up[i+1]){
	binPos = i;
	break;
      }
    }
    else if (down) {
      if(hiHF >= binTable_down[i] && hiHF < binTable_down[i+1]){
	binPos = i;
	break;
      }
    }
  }
  
  binPos = ncBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)ncBins));
  
}

DEFINE_FWK_MODULE(SayanCMW);
