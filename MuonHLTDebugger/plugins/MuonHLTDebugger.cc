// -*- C++ -*-
//
// Package:    MuonHLT/MuonHLTDebugger
// Class:      MuonHLTDebugger
// 
/**\class MuonHLTDebugger MuonHLTDebugger.cc MuonHLT/MuonHLTDebugger/plugins/MuonHLTDebugger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Santiago Folgueras
//         Created:  Thu, 22 Sep 2016 13:30:13 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"


//
// class declaration
//
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"

const std::string EFFICIENCY_SUFFIXES[2] = {"denom", "numer"};

double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};

using namespace edm;
using namespace std;


class MuonHLTDebugger : public edm::EDAnalyzer  {
public:
  explicit MuonHLTDebugger(const edm::ParameterSet&);
  ~MuonHLTDebugger();
  
  
private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) ;
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup) ;
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  // Extra Methods
  
  // Member Variables
  HLTConfigProvider hltConfig_;

  //  edm::InputTag vertexTag_;
  //  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  
  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

  edm::InputTag l3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_; 
  edm::InputTag l2candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l2candToken_; 
  edm::InputTag l1candTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; 

  //  edm::InputTag muonTag_;
  //  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  
  // Trigger process
  edm::InputTag triggerResultsTag_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultsToken_;
  edm::InputTag triggerSummaryTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken_;
  
  // Trigger process
  std::string triggerProcess_;
  std::string triggerName_;
  std::string l1filterLabel_;
  std::string l2filterLabel_;
  std::string l3filterLabel_;

  
  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  
  edm::Service<TFileService> outfile_;
  std::map<std::string, TH1*> hists_;     
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonHLTDebugger::MuonHLTDebugger(const edm::ParameterSet& cfg):
  genTag_                 (cfg.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
  genToken_               (consumes<reco::GenParticleCollection>(genTag_)),
  l3candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L3Candidates")),
  l3candToken_            (consumes<reco::RecoChargedCandidateCollection>(l3candTag_)),
  l2candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L2Candidates")),
  l2candToken_            (consumes<reco::RecoChargedCandidateCollection>(l2candTag_)),
  l1candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L1Candidates")),
  l1candToken_            (consumes<l1t::MuonBxCollection>(l1candTag_)),
  triggerResultsTag_      (cfg.getUntrackedParameter<edm::InputTag>("triggerResults")), 
  triggerResultsToken_    (consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerSummaryTag_      (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")), 
  triggerSummaryToken_    (consumes<trigger::TriggerEvent>(triggerSummaryTag_)),
  triggerProcess_         (cfg.getParameter<std::string>("triggerProcess")), 
  triggerName_            (cfg.getParameter<std::string>("triggerName")), 
  l1filterLabel_          (cfg.getParameter<std::string>("l1filterLabel")), 
  l2filterLabel_          (cfg.getParameter<std::string>("l2filterLabel")), 
  l3filterLabel_          (cfg.getParameter<std::string>("l3filterLabel"))
{
  
}


MuonHLTDebugger::~MuonHLTDebugger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called for each event  ------------
void
MuonHLTDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // capable of reading edm::Views
  //  if (DEBUGLEVEL==1) std::cout << std::endl << "========== Analysing EVENT: " << iEvent.id() << "====" << std::endl;


  /// READING THE OBJECTS
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genToken_, genParticles);
    
  //########################### Trigger Info ###########################
  // Get objects from the event.  
  Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByToken(triggerSummaryToken_, triggerSummary);

  if(!triggerSummary.isValid()) 
    {
      LogError("MuonHLTDebugger")<<"Missing triggerSummary collection" << endl;
      return;
    }

  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);

  if(!triggerResults.isValid()) 
    {
      LogError("MuonHLTDebugger")<<"Missing triggerResults collection" << endl;
      return;
    }


  //Get the Reconstructed object collections:
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  iEvent.getByToken(l1candToken_,l1Muons);
  edm::Handle <reco::RecoChargedCandidateCollection> l2Muons;
  iEvent.getByToken(l2candToken_,l2Muons);
  edm::Handle <reco::RecoChargedCandidateCollection> l3Muons;
  iEvent.getByToken(l3candToken_, l3Muons);

  //Get filter objects, these are the names of the paths for the Mu50 path:
  size_t L1MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l1filterLabel_,"",triggerProcess_));//The L1 Filter
  size_t L2MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l2filterLabel_,"",triggerProcess_));//The L2 Filter
  size_t L3MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l3filterLabel_,"",triggerProcess_));//The L3 Filter

  trigger::TriggerObjectCollection L1MuonTrigObjects;
  trigger::TriggerObjectCollection L2MuonTrigObjects;
  trigger::TriggerObjectCollection L3MuonTrigObjects;
  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();

  if (L1MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L1MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L1MuonTrigObjects.push_back(foundObject);
    }
  }
  if (L2MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L2MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L2MuonTrigObjects.push_back(foundObject);
    }
  }
  if (L3MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L3MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L3MuonTrigObjects.push_back(foundObject);
    }
  }

  std::cout << "Number of L1s passing filter = " << L1MuonTrigObjects.size() << std::endl;
  std::cout << "Number of L2s passing filter = " << L2MuonTrigObjects.size() << std::endl;
  std::cout << "Number of L3s passing filter = " << L3MuonTrigObjects.size() << std::endl;

  std::cout << "Num L1: " << l1Muons->size() << std::endl;
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonHLTDebugger::beginJob()
{
  TH1::SetDefaultSumw2();
  
  // global quantities 
  hists_["hlt_pt"]  = outfile_->make<TH1F>("hlt_pt",  "HLT p_{T}; p_{T} of HLT object", 11,  pt_bins );
  hists_["hlt_eta"] = outfile_->make<TH1F>("hlt_eta", "HLT #eta; #eta of HLT object", 15, eta_bins );
  hists_["hlt_phi"] = outfile_->make<TH1F>("hlt_phi", "HLT #phi;#phi of HLT object", 15, -3.2, 3.2);
  hists_["nVtx"]    = outfile_->make<TH1F>("nvtx_event","nvtx_event"       ,   60,  0,   60 );
  
  hists_["resolutionEta"] = outfile_->make<TH1F>("resolutionEta", "Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.01,   0.01);
  hists_["resolutionPhi"] = outfile_->make<TH1F>("resolutionPhi", "Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.01,   0.01);
  hists_["resolutionPt"]  = outfile_->make<TH1F>("resolutionPt", "Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.30,   0.30);
 
  hists_["efficiencyEta_num" ] = outfile_->make<TH1F>("efficiencyEta_num", "Efficiency;#eta",  15, eta_bins );
  hists_["efficiencyPhi_num" ] = outfile_->make<TH1F>("efficiencyPhi_num", "Efficiency;#phi",  15, -3.2, 3.2);
  hists_["efficiencyPt_num"  ] = outfile_->make<TH1F>("efficiencyPt_num",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["efficiencyNVtx_num"] = outfile_->make<TH1F>("efficiencyNVtx_num","Efficiency;NVtx",  60,  0,   60.);    

  hists_["efficiencyEta_den" ] = outfile_->make<TH1F>("efficiencyEta_den", "Efficiency;#eta",  15, eta_bins );
  hists_["efficiencyPhi_den" ] = outfile_->make<TH1F>("efficiencyPhi_den", "Efficiency;#phi",  15, -3.2, 3.2);
  hists_["efficiencyPt_den"  ] = outfile_->make<TH1F>("efficiencyPt_den",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["efficiencyNVtx_den"] = outfile_->make<TH1F>("efficiencyNVtx_den","Efficiency;NVtx",  60,  0,   60.);    
}

void 
MuonHLTDebugger::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) 
{
  // Initialize hltConfig
  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    // Now crash
    assert(false);
  }
  
  triggerIndex_ = -1; 
  
  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.find(triggerName_) != std::string::npos) {
      triggerIndex_ = int(iHltPath);
    }
    
    if( triggerIndex_>-1) break; 
  } // end for each path
  
  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }
}
// ------------ method called once each job just after ending the event loop  ------------
void MuonHLTDebugger::endJob() {}
void MuonHLTDebugger::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTDebugger);
