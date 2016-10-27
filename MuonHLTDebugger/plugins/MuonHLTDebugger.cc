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
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"


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
  
  unsigned int debuglevel_;

  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  
  edm::Service<TFileService> outfile_;
  std::map<std::string, TH1*> hists_;     
  
  edm::EDGetTokenT<TrajectorySeedCollection>  theSeedsIter0Token_;
  edm::EDGetTokenT<TrajectorySeedCollection>  theSeedsIter2Token_;
  edm::EDGetTokenT<reco::TrackCollection>           theTracksIter0Token_;
  edm::EDGetTokenT<reco::TrackCollection>           theTracksIter2Token_;
  edm::EDGetTokenT<reco::TrackCollection>           theTracksToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection>  theIOLinksToken_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > >  puSummaryInfo_;
  edm::EDGetTokenT<reco::BeamSpot>                     theBeamSpotToken_;
};

//
// constants, enums and typedefs
//
double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};

using namespace edm;
using namespace std;

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
  l3filterLabel_          (cfg.getParameter<std::string>("l3filterLabel")),
  debuglevel_             (cfg.getUntrackedParameter<unsigned int>("debuglevel")),
  theSeedsIter0Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltNewIter0HighPtTkMuPixelSeedsFromPixelTracks","","SFHLT"))),
  theSeedsIter2Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltNewIter2HighPtTkMuPixelSeeds","","SFHLT"))),
  //  theTracksNoHPIter0Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter0HighPtTkMuCtfWithMaterialTracks","","SFHLT"))),
  //  theTracksNoHPIter2Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuCtfWithMaterialTracks","","SFHLT"))),
  theTracksIter0Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter0HighPtTkMuTrackSelectionHighPurity","","SFHLT"))),
  theTracksIter2Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuTrackSelectionHighPurity","","SFHLT"))),
  theTracksToken_         (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuMerged","","SFHLT"))),
  theIOLinksToken_        (mayConsume<reco::MuonTrackLinksCollection>(edm::InputTag("NewIterL3Muons","","SFHLT"))),
  puSummaryInfo_          (consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("addPileupInfo"))),
  theBeamSpotToken_       (consumes<reco::BeamSpot>(edm::InputTag("hltOnlineBeamSpot")))
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
    

  edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
  iEvent.getByToken(puSummaryInfo_,puInfo);
  int nbMCvtx = -1;
  if (puInfo.isValid()) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing()==0){
	nbMCvtx = PVI->getPU_NumInteractions();
	break;
      }
    }
  }

  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(theBeamSpotToken_,recoBeamSpotHandle);
  const reco::BeamSpot& beamSpot = *recoBeamSpotHandle;

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
  if (debuglevel_ > 1) {
    std::cout << "Number of L1s passing filter = " << L1MuonTrigObjects.size() << std::endl;
    std::cout << "Number of L2s passing filter = " << L2MuonTrigObjects.size() << std::endl;
    std::cout << "Number of L3s passing filter = " << L3MuonTrigObjects.size() << std::endl;
  }

  // Check by-module... 
  //  const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
  //  for (uint32_t i = 0; i < moduleLabels.size(); ++i){    
  //    size_t filterIndex = (*triggerSummary).filterIndex(edm::InputTag( moduleLabels[i].c_str(), "", triggerProcess_));
  //    
  //    if (filterIndex < (*triggerSummary).sizeFilters()) {
  //      const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
  //      for (size_t j = 0; j < keys.size(); j++ ){
  //	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
  //	if (foundObject.pt()>0) hists_["hlt_counter"]->Fill(i+1);
  //      }
  //    }
  //  }
  
  double NumL1MatchedToGenInEvent=0;
  double NumL2MatchedToGenInEvent=0;
  double NumL3MatchedToGenInEvent=0;

  std::vector<const trigger::TriggerObject*> foundL3Trig;

  std::vector<const reco::GenParticle*> L1FoundGens;
  std::vector<const reco::GenParticle*> L2FoundGens;
  std::vector<const reco::GenParticle*> L3FoundGens;
  std::vector<const reco::GenParticle*> L3FoundWrtL1Gens;

  // Loop over muons and fill histograms: 
  int numGenPerEvent=0;
  for (unsigned int g(0); g < genParticles->size(); ++g){
    const reco::GenParticle* gen = &genParticles->at(g);
    if (fabs(gen->pdgId())!=13) continue;
    if (gen->pt()<10)           continue;
    if (gen->status()!=1)       continue;
    if (fabs(gen->eta())>2.4)   continue;
    ++numGenPerEvent;
    
    hists_["gen_pt"] ->Fill(gen->pt());
    hists_["gen_eta"]->Fill(gen->eta());
    hists_["gen_phi"]->Fill(gen->phi());

    if (debuglevel_ > 1) 
      std::cout << "gen muon found: pt: " << gen->pt() 
		<< " eta: " << gen->eta() 
		<< " phi: " << gen->phi() << std::endl;
    
    double NumL1MatchedToGen=0;
    double NumL2MatchedToGen=0;
    double NumL3MatchedToGen=0;
    
    int numL2Found=0;
    int numL3Found=0;
    int numL3NotFound=0;
    
    /// DRAWING RESOLUTION PLOTS: 
    for (int ibx = l1Muons->getFirstBX(); ibx <= l1Muons->getLastBX(); ++ibx) {
      if (ibx != 0) continue; 
      for (auto it = l1Muons->begin(ibx); it != l1Muons->end(ibx); it++){
	l1t::MuonRef l1muon(l1Muons, distance(l1Muons->begin(l1Muons->getFirstBX()),it) );
	if (gen->charge() > 0 && l1muon->charge()>0) {
	  hists_["hltL1p_resEta"]->Fill(l1muon->eta()/gen->eta()); 
	  hists_["hltL1p_resPhi"]->Fill(l1muon->phi()/gen->phi());
	  hists_["hltL1p_resPt"] ->Fill(l1muon->pt() /gen->pt() );

	  if (fabs(gen->eta()) < 0.9) {
	    hists_["hltL1p_resEta_barrel"]->Fill(l1muon->eta()/gen->eta());
	    hists_["hltL1p_resPhi_barrel"]->Fill(l1muon->phi()/gen->phi());
	  }
	  if (fabs(gen->eta()) > 0.9) {
	    hists_["hltL1p_resEta_endcap"]->Fill(l1muon->eta()/gen->eta());
	    hists_["hltL1p_resPhi_endcap"]->Fill(l1muon->phi()/gen->phi());
	  }
	}
	else if (gen->charge() < 0 && l1muon->charge()<0)  {
	  hists_["hltL1m_resEta"]->Fill(l1muon->eta()/gen->eta()); 
	  hists_["hltL1m_resPhi"]->Fill(l1muon->phi()/gen->phi());
	  hists_["hltL1m_resPt"] ->Fill(l1muon->pt() /gen->pt() );

	  if (fabs(gen->eta()) < 0.9) {
	    hists_["hltL1m_resEta_barrel"]->Fill(l1muon->eta()/gen->eta());
	    hists_["hltL1m_resPhi_barrel"]->Fill(l1muon->phi()/gen->phi());
	  }
	  if (fabs(gen->eta()) > 0.9) {
	    hists_["hltL1m_resEta_endcap"]->Fill(l1muon->eta()/gen->eta());
	    hists_["hltL1m_resPhi_endcap"]->Fill(l1muon->phi()/gen->phi());
	  }
	}

      }
    }
    
    
    //L1 Match:
    bool foundL1=false;
    for (unsigned int t(0); t < L1MuonTrigObjects.size(); ++t){
      trigger::TriggerObject* l1mu = &L1MuonTrigObjects.at(t);
      if (debuglevel_ > 1) std::cout << "\tL1["<<t<<"]: deltaR(*gen,*l1mu): " << deltaR(*gen,*l1mu) << std::endl;
      hists_["hltL1_DeltaR"]->Fill(deltaR(*gen,*l1mu));
      hists_["hltL1_pt"]    ->Fill(l1mu->pt());
      hists_["hltL1_eta"]   ->Fill(l1mu->eta());
      hists_["hltL1_phi"]   ->Fill(l1mu->phi());
  
      if (deltaR(*gen,*l1mu)<0.4 && !foundL1){
	if (debuglevel_ > 1)  std::cout << "\t\tL1 found: pt: " << l1mu->pt() << " eta: " << l1mu->eta() << std::endl;
	hists_["effL3L1Eta_den" ] ->Fill(gen->eta());
	hists_["effL3L1Phi_den" ] ->Fill(gen->phi());
	hists_["effL3L1Pt_den"  ] ->Fill(gen->pt());
	hists_["effL3L1NPU_den" ] ->Fill(nbMCvtx);
	//	hists_["effL3L1NVtx_den"] ->Fill(
	++NumL1MatchedToGen;
	++NumL1MatchedToGenInEvent;
	foundL1=true;
	L1FoundGens.push_back(gen);
      }
    }
    
    //L2 Match:
    bool foundL2=false;
    for (unsigned int t(0); t < L2MuonTrigObjects.size(); ++t){
      trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(t);
      if (debuglevel_ > 1)  std::cout << "\tL2["<<t<<"] deltaR(*gen,*l2mu): " << deltaR(*gen,*l2mu) << std::endl;

      if (gen->charge() > 0) {
	hists_["hltL2p_resEta"]->Fill(l2mu->eta()-gen->eta()); 
	hists_["hltL2p_resPhi"]->Fill(l2mu->phi()-gen->phi());
	hists_["hltL2p_resPt"] ->Fill(l2mu->pt() -gen->pt() );
	
	if (fabs(gen->eta()) < 0.9) {
	  hists_["hltL2p_resEta_barrel"]->Fill(l2mu->eta()-gen->eta());
	  hists_["hltL2p_resPhi_barrel"]->Fill(l2mu->phi()-gen->phi());
	}
	if (fabs(gen->eta()) > 0.9) {
	  hists_["hltL2p_resEta_endcap"]->Fill(l2mu->eta()-gen->eta());
	  hists_["hltL2p_resPhi_endcap"]->Fill(l2mu->phi()-gen->phi());
	}
      }
      else if (gen->charge() < 0)  {
	hists_["hltL2m_resEta"]->Fill(l2mu->eta()-gen->eta()); 
	hists_["hltL2m_resPhi"]->Fill(l2mu->phi()-gen->phi());
	hists_["hltL2m_resPt"] ->Fill(l2mu->pt() -gen->pt() );
	
	if (fabs(gen->eta()) < 0.9) {
	  hists_["hltL2m_resEta_barrel"]->Fill(l2mu->eta()-gen->eta());
	  hists_["hltL2m_resPhi_barrel"]->Fill(l2mu->phi()-gen->phi());
	}
	if (fabs(gen->eta()) > 0.9) {
	  hists_["hltL2m_resEta_endcap"]->Fill(l2mu->eta()-gen->eta());
	  hists_["hltL2m_resPhi_endcap"]->Fill(l2mu->phi()-gen->phi());
	}
      }
      
      hists_["hltL2_DeltaR"]->Fill(deltaR(*gen,*l2mu));
      hists_["hltL2_pt"]    ->Fill(l2mu->pt());
      hists_["hltL2_eta"]   ->Fill(l2mu->eta());
      hists_["hltL2_phi"]   ->Fill(l2mu->phi());

      if (deltaR(*gen,*l2mu)<0.2 && !foundL2){
	if (debuglevel_ > 1)   std::cout << "\t\tL2 found: pt: " << l2mu->pt() << " eta: " << l2mu->eta() << std::endl;
	hists_["effL3L2Eta_den" ] ->Fill(gen->eta());
	hists_["effL3L2Phi_den" ] ->Fill(gen->phi());
	hists_["effL3L2Pt_den"  ] ->Fill(gen->pt());
	hists_["effL3L2NPU_den" ] ->Fill(nbMCvtx);
	foundL2=true;
	++numL2Found;	  
	++NumL2MatchedToGen;
	++NumL2MatchedToGenInEvent;
	L2FoundGens.push_back(gen);
      }
    }

    if (foundL2){
      //L3 Match:
      bool foundL3=false;
      std::vector<const trigger::TriggerObject*> foundL3Trig;    
      for (unsigned int t(0); t < L3MuonTrigObjects.size(); ++t){
	trigger::TriggerObject* l3mu = &L3MuonTrigObjects.at(t);
	hists_["hltL3_DeltaR"]->Fill(deltaR(*gen,*l3mu));
	hists_["hltL3_pt"]    ->Fill(l3mu->pt());
	hists_["hltL3_eta"]   ->Fill(l3mu->eta());
	hists_["hltL3_phi"]   ->Fill(l3mu->phi());
	hists_["hltL3_resEta"]->Fill(l3mu->eta()-gen->eta());
	hists_["hltL3_resPhi"]->Fill(l3mu->phi()-gen->phi());
	hists_["hltL3_resPt"] ->Fill(l3mu->pt() -gen->pt() );
	if (debuglevel_ > 1)  std::cout << "\tL3["<<t<<"]: deltaR(*gen,*l3mu): " << deltaR(*gen,*l3mu) << std::endl;
	if (deltaR(*gen,*l3mu)<0.01 && !foundL3){
	  if (std::find(foundL3Trig.begin(), foundL3Trig.end(), l3mu)!=foundL3Trig.end()) std::cout << "THIS L3 WAS ALREADY FOUND!" << std::endl;
	  foundL3Trig.push_back(l3mu);
	  if (debuglevel_ > 1)  std::cout << "\t\tL3 found: pt: " << l3mu->pt() << " eta: " << l3mu->eta() << std::endl;
	  hists_["effL3Eta_num" ] ->Fill(gen->eta());
	  hists_["effL3Phi_num" ] ->Fill(gen->phi());
	  hists_["effL3Pt_num"  ] ->Fill(gen->pt());
	  hists_["effL3NPU_num" ] ->Fill(nbMCvtx);
	  ++NumL3MatchedToGen;
	  ++NumL3MatchedToGenInEvent;
	  L3FoundGens.push_back(gen);
	  foundL3=true;
	  ++numL3Found;
	}
	else {
	  ++numL3NotFound;
	}
      }
    }

    hists_["hlt_NumL1Match" ]->Fill(NumL1MatchedToGen);
    hists_["hlt_NumL2Match" ]->Fill(NumL2MatchedToGen);
    hists_["hlt_NumL3Match" ]->Fill(NumL3MatchedToGen);
    hists_["hlt_NumL2" ]     ->Fill(numL2Found);
    hists_["hlt_NumL3" ]     ->Fill(numL3Found);
    hists_["hlt_NumNoL3" ]   ->Fill(numL3NotFound);
    
    if (debuglevel_ > 1)   std::cout << "numL2Found: " << numL2Found
				     << " numL3Found: " << numL3Found
				     << " numL3NotFound: " << numL3NotFound
				     << std::endl;
  } //genMuon
    
  /// NOW FOR DIMUON EVENTS CHECK ALSO THE PER-EVENT-EFFICIENCY
  if (L2FoundGens.size()==2){
    if (debuglevel_ > 1) std::cout << "Yes there are 2 L2FoundGens: " << L2FoundGens.size() << std::endl;
    for (unsigned int t(0); t<L2FoundGens.size(); ++t){
      const reco::GenParticle* gen = L2FoundGens.at(t);
      hists_["eff2L3L2Eta_den"]->Fill(gen->eta());
      hists_["eff2L3L2Phi_den"]->Fill(gen->phi());
      hists_["eff2L3L2Pt_den"] ->Fill(gen->pt() );
      hists_["eff2L3L2NPU_den"]->Fill(nbMCvtx);
      if (L3FoundGens.size()==2){
	hists_["eff2L3L2Eta_num"]->Fill(gen->eta());
	hists_["eff2L3L2Phi_num"]->Fill(gen->phi());
	hists_["eff2L3L2Pt_num"] ->Fill(gen->pt() );
	hists_["eff2L3L2NPU_num"]->Fill(nbMCvtx);
      }
    }
  }
  if (L1FoundGens.size()==2){
    if (debuglevel_ > 1)  std::cout << "Yes there are 2 L1FoundGens: " << L1FoundGens.size() << std::endl;
    for (unsigned int t(0); t<L1FoundGens.size();++t){
      const reco::GenParticle* gen = L1FoundGens.at(t);
      hists_["eff2L3L1Eta_den"]->Fill(gen->eta());
      hists_["eff2L3L1Phi_den"]->Fill(gen->phi());
      hists_["eff2L3L1Pt_den"] ->Fill(gen->pt() );
      hists_["eff2L3L1NPU_den"]->Fill(nbMCvtx);
      if (L3FoundGens.size()==2){
	hists_["eff2L3L1Eta_num"]->Fill(gen->eta());
	hists_["eff2L3L1Phi_num"]->Fill(gen->phi());
	hists_["eff2L3L1Pt_num"] ->Fill(gen->pt() );
	hists_["eff2L3L1NPU_num"]->Fill(nbMCvtx);
      }
    }
  }

  if (debuglevel_ > 1) std::cout << "============================================================" << std::endl
				 << "============================================================" << std::endl;

  
  try {
    /// NOW do other stuff::
    if (L2MuonTrigObjects.size() > 0) {
      if (L3MuonTrigObjects.size()>0) {
	// get the number of seeds... 
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter0;
	iEvent.getByToken(theSeedsIter0Token_, hltL3TrajSeedIter0);
	hists_["hlt_numSeedsIter0"]->Fill(hltL3TrajSeedIter0->size());
	hists_["hlt_numSeedsIter0_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter0->size());
	if (debuglevel_ > 1)    std::cout << "# of hltL3TrajSeedIter0: " << hltL3TrajSeedIter0->size() << std::endl;
	
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter2;
	iEvent.getByToken(theSeedsIter2Token_, hltL3TrajSeedIter2);
	hists_["hlt_numSeedsIter2"]->Fill(hltL3TrajSeedIter2->size());
	hists_["hlt_numSeedsIter2_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter2->size());
	if (debuglevel_ > 1)    std::cout << "# of hltL3TrajSeedIter2: " << hltL3TrajSeedIter2->size() << std::endl;
	
	hists_["hlt_numSeeds"]->Fill(hltL3TrajSeedIter0->size()+hltL3TrajSeedIter2->size());
	hists_["hlt_numSeeds_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter0->size()+hltL3TrajSeedIter2->size());
	
	// now for the number of tracks
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter0;
	iEvent.getByToken(theTracksIter0Token_, hltL3TkTracksIter0);
	hists_["hlt_numTracksIter0"]->Fill(hltL3TkTracksIter0->size());
	if (debuglevel_ > 1)     std::cout << "# of hltL3TkTracksIter0 = " << hltL3TkTracksIter0->size() << std::endl;    
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter0->begin(); t!=hltL3TkTracksIter0->end(); t++) {
	  hists_["hlt_numPixelHitsIter0"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
	
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter2;
	iEvent.getByToken(theTracksIter2Token_, hltL3TkTracksIter2);
	hists_["hlt_numTracksIter2"]->Fill(hltL3TkTracksIter2->size());
	if (debuglevel_ > 1)     std::cout << "# of hltL3TkTracksIter2 = " << hltL3TkTracksIter2->size() << std::endl;
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter2->begin(); t!=hltL3TkTracksIter2->end(); t++) {
	  hists_["hlt_numPixelHitsIter2"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
	
	edm::Handle<reco::TrackCollection> hltL3TkTracks;
	iEvent.getByToken(theTracksToken_, hltL3TkTracks);    
	hists_["hlt_numTracks"]->Fill(hltL3TkTracks->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracks->begin(); t!=hltL3TkTracks->end(); t++) {
	  hists_["hlt_numPixelHits"]->Fill(t->hitPattern().numberOfValidPixelHits());
	  hists_["hlt_normalizedChi2"]->Fill(t->normalizedChi2()); 
	  hists_["hlt_numValidHits"]->Fill(t->numberOfValidHits());	
	  hists_["hlt_numValidMuonHits"]->Fill(t->hitPattern().numberOfValidMuonHits());
	}
	edm::Handle<reco::MuonTrackLinksCollection> IOlink;
	iEvent.getByToken(theIOLinksToken_, IOlink);
	
	for(reco::RecoChargedCandidateCollection::const_iterator cand=l3Muons->begin(); cand!=l3Muons->end(); cand++){
	  for(unsigned int l(0); l <IOlink->size(); ++l){
	    const reco::MuonTrackLinks* link = &IOlink->at(l);
	    bool useThisLink=false;
	    reco::TrackRef tk = cand->track();
	    reco::TrackRef trkTrack = link->trackerTrack();
	    
	    // Using the same method that was used to create the links
	    // ToDo: there should be a better way than dR,dPt matching
	    const reco::Track& globalTrack = *link->globalTrack();
	    float dR2 = deltaR2(tk->eta(),tk->phi(),globalTrack.eta(),globalTrack.phi());
	    float dPt = std::abs(tk->pt() - globalTrack.pt())/tk->pt();
	    
	    hists_["hltL3mu_dR2withLink"]->Fill(dR2);
	    hists_["hltL3mu_dPtwithLink"]->Fill(dPt);
	    
	    if (dR2 < 0.02*0.02 and dPt < 0.001) {
	      useThisLink=true;
	    }
	    
	    if (useThisLink) {
	      const reco::TrackRef staTrack = link->standAloneTrack();
	      hists_["hltL3mu_pt"]->Fill(cand->pt());
	      hists_["hltL3mu_eta"]->Fill(cand->eta());
	      hists_["hltL3mu_numValidHits"]->Fill(tk->numberOfValidHits());
	      hists_["hltL3mu_dr"]->Fill(std::abs( (- (cand->vx()-beamSpot.x0()) * cand->py() + (cand->vy()-beamSpot.y0()) * cand->px() ) / cand->pt() ));
	      hists_["hltL3mu_dz"]->Fill(std::abs((cand->vz()-beamSpot.z0()) - ((cand->vx()-beamSpot.x0())*cand->px()+(cand->vy()-beamSpot.y0())*cand->py())/cand->pt() * cand->pz()/cand->pt()));
	      hists_["hltL3mu_dxySig"]->Fill(std::abs(tk->dxy(beamSpot.position())/tk->dxyError()));
	      hists_["hltL3mu_normalizedChi2"]->Fill(tk->normalizedChi2());
	      hists_["hltL3mu_dxy"]->Fill(std::abs(tk->dxy(beamSpot.position())));
	      hists_["hltL3mu_numValidMuonHits"]->Fill(tk->hitPattern().numberOfValidMuonHits());
	    } //end of useThisLink
	  } //end of muons in links collection
	} //end of RecoCand collection

      }   // L3 Objects
      else { 
	//      cout << "[FAILING EVENT] Run:Event" << iEvent.id().run() << " : " << iEvent.id().event() << endl;
	/// Store some information of such events failing the L3 reconstruction: 
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter0;
	iEvent.getByToken(theSeedsIter0Token_, hltL3TrajSeedIter0);
	hists_["hlt_noL3_numSeedsIter0"]->Fill(hltL3TrajSeedIter0->size());
	if (debuglevel_ > 1)    std::cout << "# of hltL3TrajSeedIter0: " << hltL3TrajSeedIter0->size() << std::endl;
      
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter2;
	iEvent.getByToken(theSeedsIter2Token_, hltL3TrajSeedIter2);
	hists_["hlt_noL3_numSeedsIter2"]->Fill(hltL3TrajSeedIter2->size());
	if (debuglevel_ > 1)    std::cout << "# of hltL3TrajSeedIter2: " << hltL3TrajSeedIter2->size() << std::endl;
      
	// now for the number of tracks
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter0;
	iEvent.getByToken(theTracksIter0Token_, hltL3TkTracksIter0);
	hists_["hlt_noL3_numTracksIter0"]->Fill(hltL3TkTracksIter0->size());
	if (debuglevel_ > 1)     std::cout << "# of hltL3TkTracksIter0 = " << hltL3TkTracksIter0->size() << std::endl;    
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter0->begin(); t!=hltL3TkTracksIter0->end(); t++) {
	  hists_["hlt_noL3_numPixelHitsIter0"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
      
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter2;
	iEvent.getByToken(theTracksIter2Token_, hltL3TkTracksIter2);
	hists_["hlt_noL3_numTracksIter2"]->Fill(hltL3TkTracksIter2->size());
	if (debuglevel_ > 1)     std::cout << "# of hltL3TkTracksIter2 = " << hltL3TkTracksIter2->size() << std::endl;
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter2->begin(); t!=hltL3TkTracksIter2->end(); t++) {
	  hists_["hlt_noL3_numPixelHitsIter2"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
      
	edm::Handle<reco::TrackCollection> hltL3TkTracks;
	iEvent.getByToken(theTracksToken_, hltL3TkTracks);    
	hists_["hlt_noL3_numTracks"]->Fill(hltL3TkTracks->size());
	trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(0);  
	for (reco::TrackCollection::const_iterator t=hltL3TkTracks->begin(); t!=hltL3TkTracks->end(); t++) {
	  hists_["hlt_noL3_numPixelHits"]    ->Fill(t->hitPattern().numberOfValidPixelHits());
	  hists_["hlt_noL3_normalizedChi2"]  ->Fill(t->normalizedChi2()); 
	  hists_["hlt_noL3_numValidHits"]    ->Fill(t->numberOfValidHits());
	  hists_["hlt_noL3_numValidMuonHits"]->Fill(t->hitPattern().numberOfValidMuonHits());
	  hists_["hlt_noL3_trackpt"]         ->Fill(t->pt());
	  hists_["hlt_noL3_tracketa"]        ->Fill(t->eta());
	  hists_["hlt_noL3_trackphi"]        ->Fill(t->phi());
	  hists_["hlt_noL3_DeltaR"]          ->Fill(deltaR(*t,*l2mu));
	  hists_["hlt_noL3_DeltaEta"]        ->Fill(t->eta()-l2mu->eta());
	  hists_["hlt_noL3_DeltaPhi"]        ->Fill(t->phi()-l2mu->phi());
	  hists_["hlt_noL3_DeltaPt"]         ->Fill((t->pt()-l2mu->pt())/l2mu->pt());
	}
	/*
	edm::Handle<reco::MuonTrackLinksCollection> IOlink;
	iEvent.getByToken(theIOLinksToken_, IOlink);
	
	for(reco::RecoChargedCandidateCollection::const_iterator cand=l3Muons->begin(); cand!=l3Muons->end(); cand++){
	  for(unsigned int l(0); l <IOlink->size(); ++l){
	    const reco::MuonTrackLinks* link = &IOlink->at(l);
	    bool useThisLink=false;
	    reco::TrackRef tk = cand->track();
	    reco::TrackRef trkTrack = link->trackerTrack();
	    
	    // Using the same method that was used to create the links
	    // ToDo: there should be a better way than dR,dPt matching
	    const reco::Track& globalTrack = *link->globalTrack();
	    float dR2 = deltaR2(tk->eta(),tk->phi(),globalTrack.eta(),globalTrack.phi());
	    float dPt = std::abs(tk->pt() - globalTrack.pt())/tk->pt();
	    	    
	    if (dR2 < 0.02*0.02 and dPt < 0.001) {
	      useThisLink=true;
	    }
	    
	    //	    if (useThisLink) {
	    //	      const reco::TrackRef staTrack = link->standAloneTrack();
//	      cout << "Run:Event " << iEvent.id().run() << ":" << iEvent.id().event() << endl;
//	      cout << "      pT: " << cand->pt()  << " , eta: " << cand->eta()  << endl; 
//	      cout << "      nHits: " << tk->numberOfValidHits() << " , NMuonHits: " << tk->hitPattern().numberOfValidMuonHits() << endl;
//	      cout << "      dr: " << std::abs( (- (cand->vx()-beamSpot.x0()) * cand->py() + (cand->vy()-beamSpot.y0()) * cand->px() ) / cand->pt() ) 
		   << " , dz " << std::abs((cand->vz()-beamSpot.z0()) - ((cand->vx()-beamSpot.x0())*cand->px()+(cand->vy()-beamSpot.y0())*cand->py())/cand->pt() * cand->pz()/cand->pt()) 
		   << " , dxy " << tk->dxy(beamSpot.position()) << " , dxyError " << tk->dxyError() 
		   << " , dxy(sig) " << std::abs(tk->dxy(beamSpot.position())/tk->dxyError()) << endl;
	    } //end of useThisLink
	  } //end of muons in links collection
	} //end of RecoCand collection
	*/
      } // NO L3 Objects
    } //L2Objects
  }
  catch (...) {
  }
}
 

// ------------ method called once each job just before starting event loop  ------------
void 
MuonHLTDebugger::beginJob()
{
  hists_["gen_pt"]  = outfile_->make<TH1F>("gen_pt",  "Gen Muon p_{T}; gen #mu p_{T}", 30,  pt_bins[0], pt_bins[11]);
  hists_["gen_eta"] = outfile_->make<TH1F>("gen_eta", "Gen Muon #eta; gen #mu #eta", 30, eta_bins[0], eta_bins[15]);
  hists_["gen_phi"] = outfile_->make<TH1F>("gen_phi", "Gen Muon #phi; gen #mu #phi", 30, -3.3, 3.3);
 
  // L1,L2,L3 values and efficiencies: 
  hists_["hltL1_pt"]     = outfile_->make<TH1F>("hltL1_pt",  "HLT (L1) p_{T}; p_{T} of L1 object", 11,  pt_bins );
  hists_["hltL1_eta"]    = outfile_->make<TH1F>("hltL1_eta", "HLT (L1) #eta; #eta of L1 object", 15, eta_bins );
  hists_["hltL1_phi"]    = outfile_->make<TH1F>("hltL1_phi", "HLT (L1) #phi;#phi of L1 object", 15, -3.3, 3.3);
  hists_["hltL1_DeltaR"] = outfile_->make<TH1F>("hltL1_DeltaR", "HLT (L1) #Delta R; #Delta wrt L1 object", 15, 0., 1.);
  hists_["hltL1p_resEta"] = outfile_->make<TH1F>("hltL1p_resEta", "L1 Resolution (+);#eta^{reco}-#eta^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1p_resPhi"] = outfile_->make<TH1F>("hltL1p_resPhi", "L1 Resolution (+);#phi^{reco}-#phi^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1m_resEta"] = outfile_->make<TH1F>("hltL1m_resEta", "L1 Resolution (-);#eta^{reco}-#eta^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1m_resPhi"] = outfile_->make<TH1F>("hltL1m_resPhi", "L1 Resolution (-);#phi^{reco}-#phi^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1p_resEta_barrel"] = outfile_->make<TH1F>("hltL1p_resEta_barrel", "L1 Resolution (+);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resEta_endcap"] = outfile_->make<TH1F>("hltL1p_resEta_endcap", "L1 Resolution (+);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPhi_barrel"] = outfile_->make<TH1F>("hltL1p_resPhi_barrel", "L1 Resolution (+);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPhi_endcap"] = outfile_->make<TH1F>("hltL1p_resPhi_endcap", "L1 Resolution (+);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resEta_barrel"] = outfile_->make<TH1F>("hltL1m_resEta_barrel", "L1 Resolution (-);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resEta_endcap"] = outfile_->make<TH1F>("hltL1m_resEta_endcap", "L1 Resolution (-);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resPhi_barrel"] = outfile_->make<TH1F>("hltL1m_resPhi_barrel", "L1 Resolution (-);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resPhi_endcap"] = outfile_->make<TH1F>("hltL1m_resPhi_endcap", "L1 Resolution (-);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPt"]  = outfile_->make<TH1F>("hltL1p_resPt",  "L1 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);
  hists_["hltL1m_resPt"]  = outfile_->make<TH1F>("hltL1p_resPt",  "L1 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);

  hists_["hltL2_pt"]            = outfile_->make<TH1F>("hltL2_pt",  "HLT (L2) p_{T}; p_{T} of L2 object", 11,  pt_bins );
  hists_["hltL2_eta"]           = outfile_->make<TH1F>("hltL2_eta", "HLT (L2) #eta; #eta of L2 object", 15, eta_bins );
  hists_["hltL2_phi"]           = outfile_->make<TH1F>("hltL2_phi", "HLT (L2) #phi;#phi of L2 object", 15, -3.3, 3.3);
  hists_["hltL2_DeltaR"]        = outfile_->make<TH1F>("hltL2_DeltaR", "HLT (L2) #Delta R; #Delta wrt L2 object", 15, 0., 1.);
  hists_["hltL2p_resEta"]        = outfile_->make<TH1F>("hltL2p_resEta", "L2 Resolution (+);#eta^{reco}-#eta^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2p_resPhi"]        = outfile_->make<TH1F>("hltL2p_resPhi", "L2 Resolution (+);#phi^{reco}-#phi^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2m_resEta"]        = outfile_->make<TH1F>("hltL2m_resEta", "L2 Resolution (-);#eta^{reco}-#eta^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2m_resPhi"]        = outfile_->make<TH1F>("hltL2m_resPhi", "L2 Resolution (-);#phi^{reco}-#phi^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2p_resEta_barrel"] = outfile_->make<TH1F>("hltL2p_resEta_barrel", "L2 Resolution (+);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resEta_endcap"] = outfile_->make<TH1F>("hltL2p_resEta_endcap", "L2 Resolution (+);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPhi_barrel"] = outfile_->make<TH1F>("hltL2p_resPhi_barrel", "L2 Resolution (+);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPhi_endcap"] = outfile_->make<TH1F>("hltL2p_resPhi_endcap", "L2 Resolution (+);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resEta_barrel"] = outfile_->make<TH1F>("hltL2m_resEta_barrel", "L2 Resolution (-);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resEta_endcap"] = outfile_->make<TH1F>("hltL2m_resEta_endcap", "L2 Resolution (-);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resPhi_barrel"] = outfile_->make<TH1F>("hltL2m_resPhi_barrel", "L2 Resolution (-);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resPhi_endcap"] = outfile_->make<TH1F>("hltL2m_resPhi_endcap", "L2 Resolution (-);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPt"]         = outfile_->make<TH1F>("hltL2p_resPt",         "L2 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);
  hists_["hltL2m_resPt"]         = outfile_->make<TH1F>("hltL2m_resPt",         "L2 Resolution (-);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);

  hists_["hltL3_pt"]     = outfile_->make<TH1F>("hltL3_pt",  "HLT (L3) p_{T}; p_{T} of L3 object", 11,  pt_bins );
  hists_["hltL3_eta"]    = outfile_->make<TH1F>("hltL3_eta", "HLT (L3) #eta; #eta of L3 object", 15, eta_bins );
  hists_["hltL3_phi"]    = outfile_->make<TH1F>("hltL3_phi", "HLT (L3) #phi;#phi of L3 object", 15, -3.3, 3.3);
  hists_["hltL3_DeltaR"] = outfile_->make<TH1F>("hltL3_DeltaR", "HLT (L3) #Delta R; #Delta wrt L3 object", 15, 0., 1.);
  hists_["hltL3_resEta"] = outfile_->make<TH1F>("hltL3_resEta", "L3 Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPhi"] = outfile_->make<TH1F>("hltL3_resPhi", "L3 Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPt"]  = outfile_->make<TH1F>("hltL3_resPt",  "L3 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.30,   0.30);

  hists_["effL3Eta_num" ] = outfile_->make<TH1F>("effL3Eta_num", "Efficiency;#eta",  15, eta_bins );
  hists_["effL3Phi_num" ] = outfile_->make<TH1F>("effL3Phi_num", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["effL3Pt_num"  ] = outfile_->make<TH1F>("effL3Pt_num",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["effL3NPU_num" ] = outfile_->make<TH1F>("effL3NPU_num", "Efficiency;NPU",   75,  0,   75.);

  hists_["effL3L2Eta_den" ] = outfile_->make<TH1F>("effL3L2Eta_den", "Efficiency;#eta",  15, eta_bins );
  hists_["effL3L2Phi_den" ] = outfile_->make<TH1F>("effL3L2Phi_den", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["effL3L2Pt_den"  ] = outfile_->make<TH1F>("effL3L2Pt_den",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["effL3L2NPU_den" ] = outfile_->make<TH1F>("effL3L2NPU_den", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["effL3L2NVtx_den"] = outfile_->make<TH1F>("effL3L2NVtx_den","Efficiency;NVtx",  60,  0,   60.);    

  hists_["effL3L1Eta_den" ] = outfile_->make<TH1F>("effL3L1Eta_den", "Efficiency;#eta",  15, eta_bins );
  hists_["effL3L1Phi_den" ] = outfile_->make<TH1F>("effL3L1Phi_den", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["effL3L1Pt_den"  ] = outfile_->make<TH1F>("effL3L1Pt_den",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["effL3L1NPU_den" ] = outfile_->make<TH1F>("effL3L1NPU_den", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["effL3L1NVtx_den"] = outfile_->make<TH1F>("effL3L1NVtx_den","Efficiency;NVtx",  60,  0,   60.);    

  //// DOUBLE EFFICIENCY
  hists_["eff2L3L2Eta_num" ] = outfile_->make<TH1F>("eff2L3L2Eta_num", "Efficiency;#eta",  15, eta_bins );
  hists_["eff2L3L2Phi_num" ] = outfile_->make<TH1F>("eff2L3L2Phi_num", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["eff2L3L2Pt_num"  ] = outfile_->make<TH1F>("eff2L3L2Pt_num",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["eff2L3L2NPU_num" ] = outfile_->make<TH1F>("eff2L3L2NPU_num", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["eff2L3L2NVtx_num"] = outfile_->make<TH1F>("eff2L3L2NVtx_num","Efficiency;NVtx",  60,  0,   60.);    

  hists_["eff2L3L1Eta_num" ] = outfile_->make<TH1F>("eff2L3L1Eta_num", "Efficiency;#eta",  15, eta_bins );
  hists_["eff2L3L1Phi_num" ] = outfile_->make<TH1F>("eff2L3L1Phi_num", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["eff2L3L1Pt_num"  ] = outfile_->make<TH1F>("eff2L3L1Pt_num",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["eff2L3L1NPU_num" ] = outfile_->make<TH1F>("eff2L3L1NPU_num", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["eff2L3NVtx_num"] = outfile_->make<TH1F>("eff2L3NVtx_num","Efficiency;NVtx",  60,  0,   60.);    

  hists_["eff2L3L2Eta_den" ] = outfile_->make<TH1F>("eff2L3L2Eta_den", "Efficiency;#eta",  15, eta_bins );
  hists_["eff2L3L2Phi_den" ] = outfile_->make<TH1F>("eff2L3L2Phi_den", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["eff2L3L2Pt_den"  ] = outfile_->make<TH1F>("eff2L3L2Pt_den",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["eff2L3L2NPU_den" ] = outfile_->make<TH1F>("eff2L3L2NPU_den", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["eff2L3L2NVtx_den"] = outfile_->make<TH1F>("eff2L3L2NVtx_den","Efficiency;NVtx",  60,  0,   60.);    

  hists_["eff2L3L1Eta_den" ] = outfile_->make<TH1F>("eff2L3L1Eta_den", "Efficiency;#eta",  15, eta_bins );
  hists_["eff2L3L1Phi_den" ] = outfile_->make<TH1F>("eff2L3L1Phi_den", "Efficiency;#phi",  15, -3.3, 3.3);
  hists_["eff2L3L1Pt_den"  ] = outfile_->make<TH1F>("eff2L3L1Pt_den",  "Efficiency;p_{T}", 11,  pt_bins );
  hists_["eff2L3L1NPU_den" ] = outfile_->make<TH1F>("eff2L3L1NPU_den", "Efficiency;NPU",   75,  0,   75.);
  //  hists_["eff2L3L1NVtx_den"] = outfile_->make<TH1F>("eff2L3L1NVtx_den","Efficiency;NVtx",  60,  0,   60.);    

  //// COUNTERS
  hists_["hlt_NumL1Match" ] = outfile_->make<TH1F>("hlt_NumL1Match","Number of L1 Gen Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL2Match" ] = outfile_->make<TH1F>("hlt_NumL2Match","Number of L2 Gen Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL3Match" ] = outfile_->make<TH1F>("hlt_NumL3Match","Number of L3 Gen Matched", 5, -0.5, 4.5);

  hists_["hlt_NumL2" ] = outfile_->make<TH1F>("hlt_NumL2","Number of L2 Found", 5, -0.5, 4.5);
  hists_["hlt_NumL3" ] = outfile_->make<TH1F>("hlt_NumL3","Number of L3 Found", 5, -0.5, 4.5);
  hists_["hlt_NumNoL3" ] = outfile_->make<TH1F>("hlt_NumNoL3","Number of L3 Not Found", 5, -0.5, 4.5);

  /// OTHER CHECKS:   
  hists_["hlt_numSeeds"] = outfile_->make<TH1F>("hlt_numSeeds","Number of Seeds (Iter0+Iter2)", 50, -0.5, 49.5);  
  hists_["hlt_numSeedsIter0"] = outfile_->make<TH1F>("hlt_numSeedsIter0","Number of Seeds (Iter0)", 50, -0.5, 49.5);
  hists_["hlt_numSeedsIter2"] = outfile_->make<TH1F>("hlt_numSeedsIter2","Number of Seeds (Iter2)", 50, -0.5, 49.5);

  hists_["hlt_numSeeds_PU"] = outfile_->make<TH2F>("hlt_numSeeds_PU","Number of Seeds (Iter0+Iter2) vs NPU"    , 75,  0,   75., 50, -0.5, 49.5);  
  hists_["hlt_numSeedsIter0_PU"] = outfile_->make<TH2F>("hlt_numSeedsIter0_PU","Number of Seeds (Iter0) vs NPU", 75,  0,   75., 50, -0.5, 49.5);
  hists_["hlt_numSeedsIter2_PU"] = outfile_->make<TH2F>("hlt_numSeedsIter2_PU","Number of Seeds (Iter2) vs NPU", 75,  0,   75., 50, -0.5, 49.5);

  hists_["hlt_numTracks"] = outfile_->make<TH1F>("hlt_numTracks","Number of Tracks (Iter0+Iter2)", 15, -0.5, 14.5);
  hists_["hlt_numTracksIter0"] = outfile_->make<TH1F>("hlt_numTracksIter0","Number of Tracks (Iter0)", 15, -0.5, 14.5);
  hists_["hlt_numTracksIter2"] = outfile_->make<TH1F>("hlt_numTracksIter2","Number of Tracks (Iter2)", 15, -0.5, 14.5);

  hists_["hlt_numPixelHits"]      = outfile_->make<TH1F>("hlt_numPixelHits"     ,"Number of PixelHits (Iter0+Iter2)", 50, -0.5, 49.5);
  hists_["hlt_numPixelHitsIter0"] = outfile_->make<TH1F>("hlt_numPixelHitsIter0","Number of PixelHits (Iter0)", 50, -0.5, 49.5);
  hists_["hlt_numPixelHitsIter2"] = outfile_->make<TH1F>("hlt_numPixelHitsIter2","Number of PixelHits (Iter2)", 50, -0.5, 49.5);
  hists_["hlt_numValidHits"]       = outfile_->make<TH1F>("hlt_numValidHits", "N Valid Hits Iter0+Iter2 ", 70,  -0.5, 69.5 ); 
  hists_["hlt_normalizedChi2"]     = outfile_->make<TH1F>("hlt_normalizedChi2", "Normalised Chi2 of Iter0+Iter2 ", 100, 0., 10.); 
  hists_["hlt_numValidMuonHits"]   = outfile_->make<TH1F>("hlt_numValidMuonHits", "N Valid Muon Hits of Iter0+Iter2", 70,  -0.5, 69.5 );

  /// for failing L3 
  hists_["hlt_noL3_numSeedsIter0"]     = outfile_->make<TH1F>("hlt_noL3_numSeedsIter0"    ,"Number of Seeds of Failed L3 (Iter0)",           50, -0.5, 49.5);
  hists_["hlt_noL3_numSeedsIter2"]     = outfile_->make<TH1F>("hlt_noL3_numSeedsIter2"    ,"Number of Seeds of Failed L3 (Iter2)",           50, -0.5, 49.5);

  hists_["hlt_noL3_numTracks"]         = outfile_->make<TH1F>("hlt_noL3_numTracks"        ,"Number of Tracks of Failed L3  (Iter0+Iter2)",   15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksIter0"]    = outfile_->make<TH1F>("hlt_noL3_numTracksIter0"   ,"Number of Tracks of Failed L3 (Iter0)",          15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksIter2"]    = outfile_->make<TH1F>("hlt_noL3_numTracksIter2"   ,"Number of Tracks of Failed L3 (Iter2)",          15, -0.5, 14.5);

  hists_["hlt_noL3_numPixelHits"]      = outfile_->make<TH1F>("hlt_noL3_numPixelHits"     ,"Number of PixelHits of Failed L3 (Iter0+Iter2)", 50, -0.5, 49.5);
  hists_["hlt_noL3_numPixelHitsIter0"] = outfile_->make<TH1F>("hlt_noL3_numPixelHitsIter0","Number of PixelHits of Failed L3 (Iter0)",       50, -0.5, 49.5);
  hists_["hlt_noL3_numPixelHitsIter2"] = outfile_->make<TH1F>("hlt_noL3_numPixelHitsIter2","Number of PixelHits of Failed L3 (Iter2)",       50, -0.5, 49.5);
  hists_["hlt_noL3_numValidHits"]      = outfile_->make<TH1F>("hlt_noL3_numValidHits"     ,"Number of Valid Hits of Failed L3 ",             70, -0.5, 69.5); 
  hists_["hlt_noL3_normalizedChi2"]    = outfile_->make<TH1F>("hlt_noL3_normalizedChi2"   ,"Normalised Chi2 of Failed L3 Muon",             100,  0. , 10. ); 
  hists_["hlt_noL3_numValidMuonHits"]  = outfile_->make<TH1F>("hlt_noL3_numValidMuonHits" ,"Number of Valid Muon Hits of Failed L3 Muon",    70, -0.5, 69.5);

  // other track properties...
  hists_["hlt_noL3_trackpt"]       = outfile_->make<TH1F>("hlt_noL3_trackpt"  ,"Track (L3) p_{T}; p_{T} of track",  11,  pt_bins );
  hists_["hlt_noL3_tracketa"]      = outfile_->make<TH1F>("hlt_noL3_tracketa" ,"Track (L3) #eta; #eta of L3 track",  15, eta_bins );
  hists_["hlt_noL3_trackphi"]      = outfile_->make<TH1F>("hlt_noL3_trackphi" ,"Track (L3) #phi; #phi of L3 track",  15, -3.3, 3.3);
  hists_["hlt_noL3_DeltaR"]        = outfile_->make<TH1F>("hlt_noL3_DeltaR",   "Iter L3 #Delta R; #Delta R (tk,L2)", 15, 0., 1.);
  hists_["hlt_noL3_DeltaEta"]      = outfile_->make<TH1F>("hlt_noL3_DeltaEta", "Iter L3 #Delta #eta;#eta^{tk}-#eta^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPhi"]      = outfile_->make<TH1F>("hlt_noL3_DeltaPhi", "Iter L3 #Delta #phi;#phi^{tk}-#phi^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPt"]       = outfile_->make<TH1F>("hlt_noL3_DeltaPt",  "Iter L3 #Delta p_{T};p_{T}^{tk}-p_{T}^{L2}", 60,  -0.30,   0.30);
  
  hists_["hltL3mu_dR2withLink"]      = outfile_->make<TH1F>("hltL3mu_dR2withLink"     , "#Delta R2 with Link", 100, 0, 0.02);
  hists_["hltL3mu_dPtwithLink"]      = outfile_->make<TH1F>("hltL3mu_dPtwithLink"     , "#Delta pt with Link", 100, 0, 0.05);
  hists_["hltL3mu_pt"]               = outfile_->make<TH1F>("hltL3mu_pt"              , "#p_{T} of the L3 Muon; p_{T} L3 muon", 11,  pt_bins );
  hists_["hltL3mu_eta"]		     = outfile_->make<TH1F>("hltL3mu_eta"             , "#eta of the L3 Muon; #eta L3 muons", 15, eta_bins ); 
  hists_["hltL3mu_numValidHits"]     = outfile_->make<TH1F>("hltL3mu_numValidHits"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_dr"]		     = outfile_->make<TH1F>("hltL3mu_dr"              , "Dr w.r.t. BS", 50, 0., 1. ); 
  hists_["hltL3mu_dz"]		     = outfile_->make<TH1F>("hltL3mu_dz"              , "Dz w.r.t. BS", 50, 0., 10.); 
  hists_["hltL3mu_dxySig"]	     = outfile_->make<TH1F>("hltL3mu_dxySig"          , "Dxy Significance ", 50, 0., 15.); 
  hists_["hltL3mu_dxy"]		     = outfile_->make<TH1F>("hltL3mu_dxy"             , "Dxy w.r.t. BS", 50, 0., 1.); 
  hists_["hltL3mu_normalizedChi2"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2"  , "Normalised Chi2 of L3 Muon", 100, 0., 10.); 
  hists_["hltL3mu_numValidMuonHits"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );
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
