// -*- Header -*-
//
// Package:       Analyzers/SayanCMW
// Class:         SayanCMW
//
//
// Author:        Sayan Chatterjee
// Last modified: 29/06/2022




// system include files

#include <memory>

//CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb.h"
//#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb_newEffv1.h"
//#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb_newEffv1_pTintegratedweight_w1.h"    //weight(pt,eta,centbin)
//#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb_thnsparse4D_etaptzvtxcent_newEffv1.h"  //weight(pt,eta,zvtx,centbin)
#include "Analyzers/SayanCMW/interface/DiHadronCorrelationEvt.h"


//user include files

#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TVector3.h"
#include <string>
#include "TProfile.h"
#include "TProfile2D.h"
#include <vector>
using std::vector;

class SayanCMW : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit SayanCMW(const edm::ParameterSet&);
  ~SayanCMW();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  void LoopCMWVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void AssignpTbins(double pt, double eta, double phi, float wgt,float weight, int charge, int idx);
  int GetpTbin(double pt);
  double GetDeltaEta(double eta_trg, double eta_ass);
  double GetDeltaPhi(double phi_trg, double phi_ass);
  //Int_t getHiBinFromhiHF(const Double_t hiHF,  Double_t *array);
  Int_t getHiBinFromhiHF(const Double_t hiHF, bool nom, bool up, bool down);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //****************member data*****************************************

  //## tracks ##
  edm::EDGetTokenT< edm::View < pat::PackedCandidate > > trackTags_;
  edm::EDGetTokenT< edm::View < pat::PackedGenParticle > > trackTagsgen_;
  edm::EDGetTokenT< edm::ValueMap < float > > chi2Map_;
  
  //## vertex ##     
  
  edm::EDGetTokenT < reco::VertexCollection > vtxTags_;
  
  //## centrality ##
  // used to access centrality bins
  edm::EDGetTokenT<int> cent_bin_;
  edm::EDGetTokenT < reco::Centrality > cent_reco_;
  std::vector< double > binTable;
  std::vector< double > binTable_up;
  std::vector< double > binTable_down;
  
  // ## Event selection ##
  int centmin_;
  int centmax_;
  float centBin;
  int run_no, event_no;
  
  //Event Id
  int ev_id;
  int Ntracks;
  
  // ## Vertex variables ##
  double xBestVtx_;
  double yBestVtx_;
  double zBestVtx_;
  double rhoBestVtx_;
  double xBestVtxError_;
  double yBestVtxError_;
  double zBestVtxError_;
  double zminVtx_; 
  double zmaxVtx_;
  double rhomaxVtx_;
  double nvtx_;
  
  //efficiency
  edm::InputTag fpt_;
  edm::InputTag fmb_;
  edm::InputTag fplus_;
  edm::InputTag fminus_;
  edm::InputTag fpix_;

  // ## track variables ##
  std::vector< double > pTmin_;
  std::vector< double > pTmax_;
  std::vector< double > pTmin_trg_;
  std::vector< double > pTmax_trg_;
  std::vector< double > pTmin_ass_;
  std::vector< double > pTmax_ass_;
  
  double ptmin_;
  double ptmax_;
  double etamin_;
  double etamax_;

  double pt_trigg_min;
  double pt_trigg_max;
  double eta_trigg_min;
  double eta_trigg_max;

  double pt_ass_min;
  double pt_ass_max;
  double eta_ass_min;
  double eta_ass_max;

  std::vector<double> trg_ptbin;
  std::vector<double> ass_ptbin;
  std::vector<double> etabining;
  std::vector<double> etabining5;
  std::vector<double> etabining4;
  std::vector<double> etabining3;
  std::vector<double> etabining2;
  std::vector<double> pTbining;
  std::vector<double> phibining;
  std::vector<double> centbining;
  std::vector<double> zvtxbining3;
  std::vector<double> zvtxbining2p5;
  std::vector<double> zvtxbining2;
  std::vector<double> zvtxbining1p5;
  std::vector<double> zvtxbining1;
  //std::vector<double> zvtxfull_;
  unsigned int bkgFactor;
  bool ifMcreco_;
  bool cent_nom;
  bool cent_up;
  bool cent_down;
  bool cent_2023_run3;

  // ## Dihadron corr events ##
  DiHadronCorrelationEvt* evt_;
  std::vector< DiHadronCorrelationEvt > evtVec_;
  
  TrkEff2018PbPb* TrkEff;
  TrkEff2018PbPb* TrkEff1;
  TrkEff2018PbPb* TrkEff2;
  
  // ## histograms ##
  // QA_plots
  
  TH1F* hpt;
  TH1F* hptP;
  TH1F* hptN;
  TH1F* heta;
  TH1F* hetaP;
  TH1F* hetaN;
  TH1F* h_ptbin;
  TH1F* heta_nbin;
  TH1F* hetaP_nbin;
  TH1F* hetaN_nbin;

  TH1F* hphi_nbin;
  TH1F* hphiP_nbin;
  TH1F* hphiN_nbin;

  TH1F* hpt_trigg;
  TH1F* heta_nbin_trigg;
  TH1F* hphi_nbin_trigg;

  TH1F* hpt_ass;
  TH1F* heta_nbin_ass;
  TH1F* hphi_nbin_ass;

  TH1F* hpt_trigg_check;
  TH1F* hpt_ass_check;

  TH1D* h_nHits;
  TH1D* h_pterr;
  TH1D* h_ptreso;
  TH1D* h_chi2;
  TH1D* h_DCAZ;
  TH1D* h_DCAXY;
  
  TH1F* hZBestVtx;
  TH1F* hcent_bin;
  TH1F* hcentbin_array;
  TH1F* hzvtxbin3_array;
  TH1F* hzvtxbin2p5_array;
  TH1F* hzvtxbin2_array;
  TH1F* hzvtxbin1p5_array;
  TH1F* hzvtxbin1_array;
  
  TH1F* heta_p5binwidth;
  TH1F* heta_p4binwidth;
  TH1F* heta_p3binwidth;
  TH1F* heta_p2binwidth;

  //TH1I* hntrk_pos;
  //TH1I* hntrk_neg;
  TH1F* hpt_trg;
  TH1F* hpt_asso;

  TProfile* tp1d_mpteta_nbin;
  TProfile* tp1d_mptetaP_nbin;
  TProfile* tp1d_mptetaN_nbin;
  
  //efficiency weight
  TH1F* hpt_w;
  TH1F* hptP_w;
  TH1F* hptN_w;
  TH1F* heta_w;
  TH1F* hetaP_w;
  TH1F* hetaN_w;
  TH1F* heta_nbin_w;
  TH1F* hetaP_nbin_w;
  TH1F* hetaN_nbin_w;

  TProfile* tp1d_mpteta_nbin_w;
  TProfile* tp1d_mptetaP_nbin_w;
  TProfile* tp1d_mptetaN_nbin_w;

  //trigger ptbin  
  static const int ptbn_trg =2;
  static const int ptbn =2;

  TH2D* hsignal_c2[ptbn];
  TH2D* hsignal_c2_2eff[ptbn];
  
  TH2D* hsignal_c2_mix[ptbn];
  TH2D* hsignal_c2_mix_2eff[ptbn];

  //PP
  TH2D* hsignal_c2_PP[ptbn];
  TH2D* hsignal_c2_2eff_PP[ptbn];
  
  TH2D* hsignal_c2_mix_PP[ptbn];
  TH2D* hsignal_c2_mix_2eff_PP[ptbn];

  //NN
  TH2D* hsignal_c2_NN[ptbn];
  TH2D* hsignal_c2_2eff_NN[ptbn];
  
  TH2D* hsignal_c2_mix_NN[ptbn];
  TH2D* hsignal_c2_mix_2eff_NN[ptbn];

  //PN
  TH2D* hsignal_c2_PN[ptbn];
  TH2D* hsignal_c2_2eff_PN[ptbn];
  
  TH2D* hsignal_c2_mix_PN[ptbn];
  TH2D* hsignal_c2_mix_2eff_PN[ptbn];

  //NP
  TH2D* hsignal_c2_NP[ptbn];
  TH2D* hsignal_c2_2eff_NP[ptbn];
  
  TH2D* hsignal_c2_mix_NP[ptbn];
  TH2D* hsignal_c2_mix_2eff_NP[ptbn];


  
  TH1F* h_ntrg_trigg[ptbn];
  TH1F* h_ntrg_ass;
  TH1F* h_ntrg_trigg_eff[ptbn];
  TH1F* h_ntrg_ass_eff;
    
  TH1F* h_ntrg_trigg_P[ptbn];
  TH1F* h_ntrg_trigg_eff_P[ptbn];

  TH1F* h_ntrg_trigg_N[ptbn];
  TH1F* h_ntrg_trigg_eff_N[ptbn];

  TH1D * hntrg_addbincontent[ptbn_trg];
  TH1D * hntrg_corr_addbincontent[ptbn_trg];
  
  TH1D * hntrg_addbincontent_3p0trg6p0;
  TH1D * hntrg_corr_addbincontent_3p0trg6p0;

  
  //zvtxbin
  TH2D* hsignal_c2_zvtx[ptbn_trg][ptbn];
  TH2D* hsignal_c2_zvtx_2eff[ptbn_trg][ptbn];
  
  TH2D* hsignal_c2_zvtx_mix[ptbn_trg][ptbn];
  TH2D* hsignal_c2_zvtx_mix_2eff[ptbn_trg][ptbn];

  //TH1D * hntrg_addbincontent_zvtx[3][ptbn];
  //TH1D * hntrg_corr_addbincontent_zvtx[3][ptbn];

  const Int_t ncBins = 200;
};
