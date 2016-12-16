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
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"
#include "TPRegexp.h"
#include "TString.h"
//class TrajectoryStateOnSurface;
//class TrajectorySeed;
class MuonServiceProxy;

const double NOMATCH = 999.;
const double NOMATCHITS =  0.;
const std::string EFFICIENCY_SUFFIXES[2] = {"Den", "num"};

class MuonHLTDebugger : public edm::EDAnalyzer  {
public:
  explicit MuonHLTDebugger(const edm::ParameterSet&);
  ~MuonHLTDebugger();
  
  template <class T1, class T2> std::vector<size_t> matchByDeltaR(const std::vector<T1> &, const std::vector<T2> &, const double maxDeltaR = NOMATCH);
  std::vector<size_t> matchBySharedHits(const reco::MuonCollection& Muons, trigger::TriggerObjectCollection& hltl3muons, 
					const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac=NOMATCHITS) ; 
  
  int sharedHits(const reco::Track& track1, const reco::Track& track2) const;
  
  reco::MuonCollection selectedMuons(const reco::MuonCollection &, const StringCutObjectSelector<reco::Muon> &);

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
  edm::InputTag muonTag_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  //  edm::EDGetTokenT<reco::TrackCollection>  TrackCollectionToken_;
  
  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

  edm::InputTag l3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_; 
  edm::InputTag l2candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l2candToken_; 
  edm::InputTag l1candTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; 

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
  bool isMC_;

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
  edm::EDGetTokenT<reco::MuonTrackLinksCollection>  theOILinksToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection>  theLinksToken_;
  edm::EDGetTokenT<reco::TrackCollection>  theMuonsWithHitsToken_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > >  puSummaryInfo_;
  edm::EDGetTokenT<reco::BeamSpot>                     theBeamSpotToken_;
  edm::EDGetTokenT<reco::VertexCollection>  theVertexToken_;
  edm::EDGetTokenT<TrajectorySeedCollection> theSeedsOIToken_;

  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> theB;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometry_;
  
  std::string ttrhbuilder_ ;
  
  StringCutObjectSelector<reco::Muon> targetMuonSelector_;
  MuonServiceProxy *theService;
};

//
// constants, enums and typedefs
//
double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double phi_bins[14] = {-M_PI, -M_PI*(11/12.), -M_PI*(9/12.), -M_PI*(7/12.), -M_PI*(5/12.), -M_PI*(3/12.), -M_PI*(1/12.), M_PI*(1/12.), M_PI*(3/12.), M_PI*(5/12.), M_PI*(7/12.), M_PI*(9/12.), M_PI*(11/12.), M_PI};
using namespace edm;
using namespace std;

//
// static data member definitions
//

//
// constructors and destructor
//
MuonHLTDebugger::MuonHLTDebugger(const edm::ParameterSet& cfg):
  muonTag_                (cfg.getUntrackedParameter<edm::InputTag>("muonTag")),
  muonToken_              (consumes<std::vector<reco::Muon> >(muonTag_)),
  //  TrackCollectionToken_   (consumes<reco::TrackExtraCollection>(edm::InputTag("globalMuons"))),
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
  isMC_                   (cfg.getUntrackedParameter<bool>("isMC")),
  theSeedsIter0Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltNewIter0HighPtTkMuPixelSeedsFromPixelTracks","","SFHLT"))),
  theSeedsIter2Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltNewIter2HighPtTkMuPixelSeeds","","SFHLT"))),
  //  theTracksNoHPIter0Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter0HighPtTkMuCtfWithMaterialTracks","","SFHLT"))),
  //  theTracksNoHPIter2Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuCtfWithMaterialTracks","","SFHLT"))),
  theTracksIter0Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter0HighPtTkMuTrackSelectionHighPurity","","SFHLT"))),
  theTracksIter2Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuTrackSelectionHighPurity","","SFHLT"))),
  theTracksToken_         (mayConsume<reco::TrackCollection>(edm::InputTag("hltNewIter2HighPtTkMuMerged","","SFHLT"))),
  theIOLinksToken_        (mayConsume<reco::MuonTrackLinksCollection>(edm::InputTag("NewIterL3Muons","","SFHLT"))),
  theOILinksToken_        (mayConsume<reco::MuonTrackLinksCollection>(edm::InputTag("hltL3MuonsNewOI","","SFHLT"))),
  theLinksToken_          (consumes<reco::MuonTrackLinksCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonLinksTag"))),
  theMuonsWithHitsToken_  (consumes<reco::TrackCollection>(edm::InputTag("globalMuons"))),
  puSummaryInfo_          (consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("addPileupInfo"))),
  theBeamSpotToken_       (consumes<reco::BeamSpot>(edm::InputTag("hltOnlineBeamSpot"))),
  theVertexToken_         (consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"))),
  theSeedsOIToken_        (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltNewOISeedsFromL2Muons","","SFHLT"))),
  targetMuonSelector_     ("isGlobalMuon && abs(eta) < 2.4 && pt > 10")
{
  theService = new MuonServiceProxy(cfg.getParameter<ParameterSet>("ServiceParameters"));
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
  theService->update(iSetup);
  
  /// READING THE OBJECTS
  edm::Handle<reco::GenParticleCollection> genParticles;  
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  if (isMC_) {
    iEvent.getByToken(genToken_, genParticles);
  }
  
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
  
  edm::Handle<reco::VertexCollection> pvHandle; 
  iEvent.getByToken(theVertexToken_, pvHandle);
  const reco::VertexCollection & vertices = *pvHandle.product();
  unsigned int nGoodVtx = 0; 
  for(reco::VertexCollection::const_iterator it=vertices.begin(); it!=vertices.end(); ++it) {
    if( it->ndof()>4                     && 
	(fabs(it->z())<=24.)             && 
	(fabs(it->position().rho())<=2.)   ) 
      nGoodVtx++;
  }
  if( nGoodVtx==0 ) return;
  const reco::Vertex & pv = vertices[0];
  
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(theBeamSpotToken_,recoBeamSpotHandle);
  //  const reco::BeamSpot& beamSpot = *recoBeamSpotHandle;

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

  
  
  // Loop over muons and fill histograms: 
  int numGenPerEvent=0;
  if (isMC_) { 
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


      /// DRAWING RESOLUTION PLOTS: 
      for (int ibx = l1Muons->getFirstBX(); ibx <= l1Muons->getLastBX(); ++ibx) {
	if (ibx != 0) continue; 
	for (auto it = l1Muons->begin(ibx); it != l1Muons->end(ibx); it++){
	  l1t::MuonRef l1muon(l1Muons, distance(l1Muons->begin(l1Muons->getFirstBX()),it) );
	  if (gen->charge() > 0 && l1muon->charge()>0) {
	    hists_["hltL1p_resEta"]->Fill(l1muon->eta()-gen->eta()); 
	    hists_["hltL1p_resPhi"]->Fill(l1muon->phi()-gen->phi());
	    hists_["hltL1p_resPt"] ->Fill(l1muon->pt() -gen->pt() );
	    
	    if (fabs(gen->eta()) < 0.9) {
	      hists_["hltL1p_resEta_barrel"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1p_resPhi_barrel"]->Fill(l1muon->phi()-gen->phi());
	    }
	    if (fabs(gen->eta()) > 0.9) {
	      hists_["hltL1p_resEta_endcap"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1p_resPhi_endcap"]->Fill(l1muon->phi()-gen->phi());
	    }
	  }
	  else if (gen->charge() < 0 && l1muon->charge()<0)  {
	    hists_["hltL1m_resEta"]->Fill(l1muon->eta()-gen->eta()); 
	    hists_["hltL1m_resPhi"]->Fill(l1muon->phi()-gen->phi());
	    hists_["hltL1m_resPt"] ->Fill(l1muon->pt() -gen->pt() );
	    
	    if (fabs(gen->eta()) < 0.9) {
	      hists_["hltL1m_resEta_barrel"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1m_resPhi_barrel"]->Fill(l1muon->phi()-gen->phi());
	    }
	    if (fabs(gen->eta()) > 0.9) {
	      hists_["hltL1m_resEta_endcap"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1m_resPhi_endcap"]->Fill(l1muon->phi()-gen->phi());
	    }
	  }

	} //l1Muons->begin(ibx)
      } //L1 l1Muons->getFirstBX()
      
      for (unsigned int t(0); t < L2MuonTrigObjects.size(); ++t){
	trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(t);
	
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
      } // L2
    } //genParticles
  }//isMC

  if (isMC_ && numGenPerEvent==0) return; //if no st1 muon skip the event.

  reco::MuonCollection targetMuons = selectedMuons(*muons, targetMuonSelector_);

  edm::Handle<reco::MuonTrackLinksCollection> links;
  iEvent.getByToken(theLinksToken_, links);
  
  vector<size_t> matchesL1 = matchByDeltaR(targetMuons,L1MuonTrigObjects,0.4); 
  vector<size_t> matchesL2 = matchByDeltaR(targetMuons,L2MuonTrigObjects,0.3); 
  vector<size_t> matchesL3 = matchByDeltaR(targetMuons,L3MuonTrigObjects,0.1); 
  try { 
    vector<size_t> matchesL3_hits = matchBySharedHits(targetMuons,L3MuonTrigObjects,links, 0.1); // fill with shared hits...
  }
  catch(...){
  }
  bool pairalreadyconsidered = false;
  
  //  cout << "Run:Event --> " << iEvent.id().run() << " : " << iEvent.id().event() << endl;
  for (size_t i = 0; i < targetMuons.size(); i++) {
    reco::Muon & recoMu = targetMuons[i];
    
    //    cout << matchesL2[i] << " , " << matchesL3[i] << " < " << targetMuons.size() << " ?" << endl;
    if (matchesL2[i] < targetMuons.size()) { 
      trigger::TriggerObject & l2mu = L2MuonTrigObjects[matchesL2[i]];      
      hists_["hltL2_DeltaR"]->Fill(deltaR(recoMu,l2mu));
      hists_["hltL2_pt"]    ->Fill(l2mu.pt());
      hists_["hltL2_eta"]   ->Fill(l2mu.eta());
      hists_["hltL2_phi"]   ->Fill(l2mu.phi());
      hists_["hltL2_resEta"]->Fill( recoMu.eta()-l2mu.eta());
      hists_["hltL2_resPhi"]->Fill( recoMu.phi()-l2mu.phi());
      hists_["hltL2_resPt"] ->Fill((recoMu.pt() -l2mu.pt())/recoMu.pt());
    }

    if (matchesL3[i] < targetMuons.size()) {
      trigger::TriggerObject & l3mu = L3MuonTrigObjects[matchesL3[i]];
      hists_["hltL3_DeltaR"]->Fill(deltaR(recoMu,l3mu));
      hists_["hltL3_pt"]    ->Fill(l3mu.pt());
      hists_["hltL3_eta"]   ->Fill(l3mu.eta());
      hists_["hltL3_phi"]   ->Fill(l3mu.phi());
      hists_["hltL3_resEta"]->Fill( recoMu.eta()-l3mu.eta());
      hists_["hltL3_resPhi"]->Fill( recoMu.phi()-l3mu.phi());
      hists_["hltL3_resPt"] ->Fill((recoMu.pt() -l3mu.pt())/recoMu.pt());
    }    
      
    for (size_t j = 0; j < 2; j++) {
      string suffix = EFFICIENCY_SUFFIXES[j];
      
      // L3 EFFICIENCY OVER RECO MUO...   
      // If no match was found, then the numerator plots don't get filled.
      if (suffix == "num" && matchesL3[i] >= targetMuons.size()) continue;
      hists_["effL3Eta_"+suffix ] ->Fill(recoMu.eta());
      hists_["effL3Phi_"+suffix ] ->Fill(recoMu.phi());
      hists_["effL3Pt_"+suffix  ] ->Fill(recoMu.pt());
      hists_["effL3NPU_"+suffix ] ->Fill(nbMCvtx);
      hists_["effL3NPV_"+suffix ] ->Fill(nGoodVtx);
      hists_["effL3dxy_"+suffix ] ->Fill(std::abs(recoMu.muonBestTrack()->dxy(pv.position())));
      hists_["effL3dz_"+suffix ]  ->Fill(std::abs(recoMu.muonBestTrack()->dz(pv.position())));
      hists_["effL3Chi2_"+suffix ]->Fill(recoMu.muonBestTrack()->normalizedChi2());
      
      // NOW try to compute L3/L2:  
      // If no L2 match was found, skip both denominator and numer
      if (suffix =="den" && matchesL2[i] >= targetMuons.size()) continue;
      if (matchesL3[i] != matchesL2[i]) continue;

      hists_["effL3L2Eta_"+suffix ] ->Fill(recoMu.eta());
      hists_["effL3L2Phi_"+suffix ] ->Fill(recoMu.phi());
      hists_["effL3L2Pt_"+suffix  ] ->Fill(recoMu.pt());
      hists_["effL3L2NPU_"+suffix ] ->Fill(nbMCvtx);
      hists_["effL3L2NPV_"+suffix ] ->Fill(nGoodVtx);
      hists_["effL3L2dxy_"+suffix ] ->Fill(std::abs(recoMu.muonBestTrack()->dxy(pv.position())));
      hists_["effL3L2dz_"+suffix ]  ->Fill(std::abs(recoMu.muonBestTrack()->dz(pv.position())));
      hists_["effL3L2Chi2_"+suffix ]->Fill(recoMu.muonBestTrack()->normalizedChi2());
    }

    for (size_t j = 0; j < 2; j++) {
      string suffix = EFFICIENCY_SUFFIXES[j];

      // L3 EFFICIENCY OVER RECO MUO...   
      // If no match was found, then the numerator plots don't get filled.
      if (suffix == "num" && matchesL3[i] >= targetMuons.size()) continue;
      hists_["effL3Eta_"+suffix ] ->Fill(recoMu.eta());
      hists_["effL3Phi_"+suffix ] ->Fill(recoMu.phi());
      hists_["effL3Pt_"+suffix  ] ->Fill(recoMu.pt());
      hists_["effL3NPU_"+suffix ] ->Fill(nbMCvtx);
      hists_["effL3NPV_"+suffix ] ->Fill(nGoodVtx);
      hists_["effL3dxy_"+suffix ] ->Fill(std::abs(recoMu.muonBestTrack()->dxy(pv.position())));
      hists_["effL3dz_"+suffix ]  ->Fill(std::abs(recoMu.muonBestTrack()->dz(pv.position())));
      hists_["effL3Chi2_"+suffix ]->Fill(recoMu.muonBestTrack()->normalizedChi2());
      
      // NOW try to compute L3/L2:  
      // If no L2 match was found, skip both denominator and numer
      if (suffix =="den" && matchesL2[i] >= targetMuons.size()) continue;
      if (matchesL3[i] != matchesL2[i]) continue;

      hists_["effL3L2Eta_"+suffix ] ->Fill(recoMu.eta());
      hists_["effL3L2Phi_"+suffix ] ->Fill(recoMu.phi());
      hists_["effL3L2Pt_"+suffix  ] ->Fill(recoMu.pt());
      hists_["effL3L2NPU_"+suffix ] ->Fill(nbMCvtx);
      hists_["effL3L2NPV_"+suffix ] ->Fill(nGoodVtx);
      hists_["effL3L2dxy_"+suffix ] ->Fill(std::abs(recoMu.muonBestTrack()->dxy(pv.position())));
      hists_["effL3L2dz_"+suffix ]  ->Fill(std::abs(recoMu.muonBestTrack()->dz(pv.position())));
      hists_["effL3L2Chi2_"+suffix ]->Fill(recoMu.muonBestTrack()->normalizedChi2());
    }

 
    
    // NOW for dilepton events: 
    // Muon cannot be a tag because doesn't match an hlt muon
    for (size_t k = 0; k < targetMuons.size(); k++) {
      reco::Muon & theProbe = targetMuons[k];
      if (pairalreadyconsidered) continue;
      hists_["eff2L3Eta_den"]->Fill(theProbe.eta());
      hists_["eff2L3Phi_den"]->Fill(theProbe.phi());
      hists_["eff2L3Pt_den"] ->Fill(theProbe.pt() );
      hists_["eff2L3NPU_den"]->Fill(nbMCvtx);
      hists_["eff2L3NPV_den"]->Fill(nGoodVtx);
      if (matchesL2[i] < targetMuons.size()){ 
	hists_["eff2L3L2Eta_den"]->Fill(theProbe.eta());
	hists_["eff2L3L2Phi_den"]->Fill(theProbe.phi());
	hists_["eff2L3L2Pt_den"] ->Fill(theProbe.pt() );
	hists_["eff2L3L2NPU_den"]->Fill(nbMCvtx);
	hists_["eff2L3L2NPV_den"]->Fill(nGoodVtx);
	if (matchesL3[k] < targetMuons.size()) {
	  hists_["eff2L3L2Eta_num"]->Fill(theProbe.eta());
	  hists_["eff2L3L2Phi_num"]->Fill(theProbe.phi());
	  hists_["eff2L3L2Pt_num"] ->Fill(theProbe.pt() );
	  hists_["eff2L3L2NPU_num"]->Fill(nbMCvtx);
	  hists_["eff2L3L2NPV_num"]->Fill(nGoodVtx);
	  
	  hists_["eff2L3Eta_num"]->Fill(theProbe.eta());
	  hists_["eff2L3Phi_num"]->Fill(theProbe.phi());
	  hists_["eff2L3Pt_num"] ->Fill(theProbe.pt() );
	  hists_["eff2L3NPU_num"]->Fill(nbMCvtx);
	  hists_["eff2L3NPV_num"]->Fill(nGoodVtx);
	  pairalreadyconsidered = true;
	}
      }
    }
  }
  hists_["hlt_NumL1" ]     ->Fill(L1MuonTrigObjects.size());
  hists_["hlt_NumL2" ]     ->Fill(L2MuonTrigObjects.size());
  hists_["hlt_NumL3" ]     ->Fill(L3MuonTrigObjects.size());
  hists_["hlt_NumL1Match" ]->Fill(matchesL1.size());
  hists_["hlt_NumL2Match" ]->Fill(matchesL2.size());
  hists_["hlt_NumL3Match" ]->Fill(matchesL3.size());
  
//  if (L2MuonTrigObjects.size() != L3MuonTrigObjects.size())  
//    cout << "[FAILING EVENT] Run:Event" << iEvent.id().run() << " : " << iEvent.id().event() << endl;    
//  if (matchesL3.size() != matchesL2.size()) 
//    cout << "[FAILING EVENT - No Match] Run:Event" << iEvent.id().run() << " : " << iEvent.id().event() << endl;    

//  /// NOW FOR DIMUON EVENTS CHECK ALSO THE PER-EVENT-EFFICIENCY
//  if (L2Found.size()==2){
//    for (unsigned int t(0); t<L2Found.size(); ++t){
//      const reco::Muon* mu = L2Found.at(t);
//      if (L3Found.size()==2){
//	hists_["eff2L3L2Eta_num"]->Fill(mu->eta());
//	hists_["eff2L3L2Phi_num"]->Fill(mu->phi());
//	hists_["eff2L3L2Pt_num"] ->Fill(mu->pt() );
//	hists_["eff2L3L2NPU_num"]->Fill(nbMCvtx);
//	hists_["eff2L3L2NPV_num"]->Fill(nGoodVtx);
//      }
//    }
//  }
//  if (L1Found.size()==2){
//    for (unsigned int t(0); t<L1Found.size();++t){
//      const reco::Muon* mu = L1Found.at(t);
//      hists_["eff2L3L1Eta_den"]->Fill(mu->eta());
//      hists_["eff2L3L1Phi_den"]->Fill(mu->phi());
//      hists_["eff2L3L1Pt_den"] ->Fill(mu->pt() );
//      hists_["eff2L3L1NPU_den"]->Fill(nbMCvtx);
//      hists_["eff2L3L1NPV_den"]->Fill(nGoodVtx);
//      if (L3Found.size()==2){
//	hists_["eff2L3L1Eta_num"]->Fill(mu->eta());
//	hists_["eff2L3L1Phi_num"]->Fill(mu->phi());
//	hists_["eff2L3L1Pt_num"] ->Fill(mu->pt() );
//	hists_["eff2L3L1NPU_num"]->Fill(nbMCvtx);
//	hists_["eff2L3L1NPV_num"]->Fill(nGoodVtx);
//      }
//    }

/*  try {
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
	cout << "[FAILING EVENT] Run:Event" << iEvent.id().run() << " : " << iEvent.id().event() << endl;
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
      } // NO L3 Objects
    } //L2Objects
  }
  catch (...) {
  }
  /// DEBUGGING THE OUTSIDE-IN COMPONENT 
  try {
    if (L2MuonTrigObjects.size() > 0 && L3MuonTrigObjects.size() > 0) {
      edm::Handle<TrajectorySeedCollection> hltL3TrajSeedOI;
      iEvent.getByToken(theSeedsOIToken_, hltL3TrajSeedOI);
      hists_["hlt_OI_NumSeeds"]->Fill(hltL3TrajSeedOI->size()/L2MuonTrigObjects.size());
      
      for (unsigned int t(0); t < L2MuonTrigObjects.size(); ++t){
	trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(t);
	hists_["hlt_OI_NumSeedsVsEta"]->Fill(l2mu->eta(),hltL3TrajSeedOI->size()/L2MuonTrigObjects.size());
      }
      
      edm::Handle<reco::TrackExtraCollection> TrackCollection;
      iEvent.getByToken(TrackCollectionToken_, TrackCollection);
      for (reco::TrackExtraCollection::const_iterator trk = TrackCollection->begin(); trk!=TrackCollection->end(); ++trk){
	reco::TrackExtra track = (*trk);
	
	for(TrajectorySeedCollection::const_iterator seed = hltL3TrajSeedOI->begin(); seed != hltL3TrajSeedOI->end(); ++seed){
	  TrajectorySeed::range seedHits = seed->recHits();
	  hists_["hlt_OI_NumberOfRecHitsPerSeed"]->Fill(seed->nHits());	    
	  
	  PTrajectoryStateOnDet pTSOD = seed->startingState();
	  DetId seedDetId(pTSOD.detId());
	  const GeomDet* gdet = theService->trackingGeometry()->idToDet( seedDetId );
	  TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(pTSOD, &(gdet->surface()), &*(theService)->magneticField());
	  AlgebraicSymMatrix66 errors = seedTSOS.cartesianError().matrix();
	  double partialPterror = errors(3,3)*pow(seedTSOS.globalMomentum().x(),2) + errors(4,4)*pow(seedTSOS.globalMomentum().y(),2);

	  float seedeta = seedTSOS.globalMomentum().eta(); 
	  hists_["hlt_OI_seedPtErrVsEta"]->Fill(seedeta,sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPtErr"]     ->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	  if (std::abs(seedeta) < 0.9) {
	    hists_["hlt_OI_seedPtErr_barrel"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	    hists_["hlt_OI_seedPhiErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	    hists_["hlt_OI_seedEtaErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	  }
	  else { 
	    hists_["hlt_OI_seedPtErr_endcap"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	    hists_["hlt_OI_seedPhiErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	    hists_["hlt_OI_seedEtaErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	  }	    
	  
	  for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
	    if (!(*iseed).isValid()) continue;
	    hists_["hlt_OI_hitPt"]->Fill((*iseed).globalPosition().perp());
	    hists_["hlt_OI_hitPhi"]->Fill((*iseed).globalPosition().phi());
	    hists_["hlt_OI_hitEta"]->Fill((*iseed).globalPosition().eta());
	    hists_["hlt_OI_hitx"]->Fill((*iseed).globalPosition().x());
	    hists_["hlt_OI_hity"]->Fill((*iseed).globalPosition().y());
	    hists_["hlt_OI_hitz"]->Fill((*iseed).globalPosition().z());
	    
	    float deltax(9999.); // mindx(9999.),minTOBdx(9999.),minTECdx(9999.);
	    float deltay(9999.); // mindy(9999.),minTOBdy(9999.),minTECdy(9999.);
	    
	    DetId detidSeed = (*iseed).geographicalId();
	    int subDetSeed = detidSeed.subdetId();
	    
	    for (trackingRecHit_iterator hitIt = track.recHitsBegin(); hitIt != track.recHitsEnd(); ++hitIt) {
	      if (not (*hitIt)->isValid()) continue;
	      const TrackingRecHit* rechit = (*hitIt)->hit();
	      DetId detid = rechit->geographicalId();
	      int subDet = detid.subdetId();
	      
	      // Make sure both the hit from the seed and the tracker hit are in the same detector... 
	      if (detid.det() != detidSeed.det()) continue;
	      if (subDetSeed != subDet) continue;
	      if (detid.rawId() != detidSeed.rawId()) continue; // SAME layer, rod, module (TOB) or wheel, petal, ring, module (TEC) detector.
	      
	      deltax = (*iseed).localPosition().x() - rechit->localPosition().x();
	      deltay = (*iseed).localPosition().y() - rechit->localPosition().y();
	      //  		deltaz = (*iseed).localPosition().z() - rechit->localPosition().z();
	      
	      hists_["hlt_OI_hitdx"]->Fill(deltax);
	      hists_["hlt_OI_hitdy"]->Fill(deltay);
	      if (subDet == StripSubdetector::TOB) {
		TOBDetId myDet(detid.rawId());
		hists_["hlt_OI_TOBhitdx"]->Fill(deltax);
		hists_["hlt_OI_TOBhitdy"]->Fill(deltay);
		
		if (myDet.isRPhi()) { 
		  hists_["hlt_OI_TOBmonohitdx"]->Fill(deltax);
		  hists_["hlt_OI_TOBmonohitdy"]->Fill(deltay);
		}
		if (myDet.isStereo()){
		  hists_["hlt_OI_TOBstereohitdx"]->Fill(deltax);
		  hists_["hlt_OI_TOBstereohitdy"]->Fill(deltay);
		}
	      }
	      if (subDet == StripSubdetector::TEC) {
		TECDetId myDet(detid.rawId());
		hists_["hlt_OI_TEChitdx"]->Fill(deltax);
		hists_["hlt_OI_TEChitdy"]->Fill(deltay);
		if (myDet.isRPhi()) { 
		  hists_["hlt_OI_TECmonohitdx"]->Fill(deltax);
		  hists_["hlt_OI_TECmonohitdy"]->Fill(deltay);
		}
		if (myDet.isStereo()){
		  hists_["hlt_OI_TECstereohitdx"]->Fill(deltax);
		  hists_["hlt_OI_TECstereohitdy"]->Fill(deltay);
		}
	      }
	      
		if (abs(deltax) < mindx) mindx = deltax;
		if (abs(deltay) < mindy) mindy = deltay;
		if (abs(deltaz) < mindz) mindz = deltaz;
		
		// enum SubDetector { TIB=3,TID=4,TOB=5,TEC=6 };
		if (subDet == StripSubdetector::TOB) { 
		if (abs(deltax) < minTOBdx) minTOBdx = deltax;
		if (abs(deltay) < minTOBdy) minTOBdy = deltay;
		if (abs(deltaz) < minTOBdz) minTOBdz = deltaz;
		}
		if (subDetSeed == StripSubdetector::TEC) {
		if (abs(deltax) < minTECdx) minTECdx = deltax;
		if (abs(deltay) < minTECdy) minTECdy = deltay;
		if (abs(deltaz) < minTECdz) minTECdz = deltaz;		  
		}
		}
	    }	//rechits
	  } // seeds
	} // iterator over HLT Trajectory Seed collection
      } // iterator  over muon collection 
    } // L2TriggerObject > 0 && L3TriggerObject > 0
  }
  catch (...) {
  }
*/
}
reco::MuonCollection MuonHLTDebugger::selectedMuons(const reco::MuonCollection & allMuons,  const StringCutObjectSelector<reco::Muon> &selector){ 
  
  reco::MuonCollection reducedMuons;
  for (auto const& mu : allMuons){
    const reco::Track * track = 0;
    if (mu.isTrackerMuon())         track = & * mu.innerTrack();
    else if (mu.isStandAloneMuon()) track = & * mu.outerTrack();
    if (track && selector(mu))
      reducedMuons.push_back(mu);
  }
  return reducedMuons;
}

int MuonHLTDebugger::sharedHits(const reco::Track& track1, const reco::Track& track2) const {
  
  int match = 0;
  for (trackingRecHit_iterator hit1 = track1.recHitsBegin(); hit1 != track1.recHitsEnd(); ++hit1) {
    if ( !(*hit1)->isValid() ) continue;
    DetId id1 = (*hit1)->geographicalId();
    GlobalPoint pos1 = theService->trackingGeometry()->idToDet(id1)->surface().toGlobal((*hit1)->localPosition());
    for (trackingRecHit_iterator hit2 = track2.recHitsBegin(); hit2 != track2.recHitsEnd(); ++hit2) {
      if ( !(*hit2)->isValid() ) continue;
      DetId id2 = (*hit2)->geographicalId();
      if (id2.subdetId() != id1.subdetId()) continue;
      if (id2.rawId() != id1.rawId() ) continue;
      if (id2.det() == DetId::Muon ) continue;
      GlobalPoint pos2 = theService->trackingGeometry()->idToDet(id2)->surface().toGlobal((*hit2)->localPosition());
      double diff = ( pos1 - pos2 ).mag();
      cout << "GLOBAL: Pos1 - Pos2 = " << pos1 << " - " << pos2 << " = " << diff << endl;
      cout << "LOCAL:  Pos1 - Pos2 = " << (*hit1)->localPosition() <<" - "<< (*hit2)->localPosition() << " = " << ((*hit1)->localPosition()-(*hit2)->localPosition()).mag()<< endl;
      cout << "Shares input? " << (*hit1)->sharesInput(*hit2,TrackingRecHit::some) << endl;

      
      //hists_["hltL3mu_DiffHits"] ->Fill(diff);
      if ( diff < 1 ) match++;
    }
  }
  cout << "N Shared Hits: " << match << " / " << track1.numberOfValidHits() << " --> " << track1.pt() << " , " << track1.eta() << " , " << track1.phi() << endl;
  cout << "N Shared Hits: " << match << " / " << track2.numberOfValidHits() << " --> " << track2.pt() << " , " << track2.eta() << " , " << track2.phi() << endl;
  return match;
}

template <class T1, class T2> vector<size_t> MuonHLTDebugger::matchByDeltaR(const vector<T1> & collection1, const vector<T2> & collection2, const double maxDeltaR)  {
  
  const size_t n1 = collection1.size();
  const size_t n2 = collection2.size();
  
  vector<size_t> result(n1, -1);
  vector<vector<double> > deltaRMatrix(n1, vector<double>(n2, NOMATCH));
  
  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
      deltaRMatrix[i][j] = deltaR(collection1[i], collection2[j]);
    }
  
  // Run through the matrix n1 times to make sure we've found all matches.
  for (size_t k = 0; k < n1; k++) {
    size_t i_min = -1;
    size_t j_min = -1;
    double minDeltaR = maxDeltaR;
    // find the smallest deltaR
    for (size_t i = 0; i < n1; i++)
      for (size_t j = 0; j < n2; j++)
	if (deltaRMatrix[i][j] < minDeltaR) {
	  i_min = i;
	  j_min = j;
	  minDeltaR = deltaRMatrix[i][j];
	}
    
    if (minDeltaR < maxDeltaR) {
      result[i_min] = j_min;
      deltaRMatrix[i_min] = vector<double>(n2, NOMATCH);
      for (size_t i = 0; i < n1; i++)
	deltaRMatrix[i][j_min] = NOMATCH;
    }
  }
  return result;
}

vector<size_t> MuonHLTDebugger::matchBySharedHits(const reco::MuonCollection& muons, trigger::TriggerObjectCollection& hltl3muons,
						  const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac)  
{
  
  const size_t n1 = muons.size();
  const size_t n2 = hltl3muons.size();
  
  vector<size_t> result(n1, -1);
  vector<vector<double> > SharedHitMatrix(n1, vector<double>(n2, NOMATCHITS));

  for (size_t i = 0; i < n1; i++) {
    const reco::Muon  recoMu = muons[i];
    reco::TrackRef mutk = recoMu.globalTrack();
    
    for (size_t j = 0; j < n2; j++) {
      trigger::TriggerObject & l3mu = hltl3muons[j];
      
      for(unsigned int l(0); l <links->size(); ++l){
	const reco::MuonTrackLinks* link = &links->at(l);
	
	const reco::Track& globalTrack = *link->globalTrack();
	float dR2 = deltaR2(l3mu.eta(),l3mu.phi(),globalTrack.eta(),globalTrack.phi());
	float dPt = std::abs(l3mu.pt() - globalTrack.pt())/l3mu.pt();
	
	if (dR2 < 0.02*0.02 and dPt < 0.001) {
	  reco::TrackRef tk = link->globalTrack();
	  SharedHitMatrix[i][j] = sharedHits(*mutk, *tk)/recoMu.innerTrack()->numberOfValidHits();
	  hists_["hltL3mu_numValidHits"]->Fill(tk->numberOfValidHits());
	  hists_["hltL3mu_normalizedChi2"]->Fill(tk->normalizedChi2());
	  hists_["hltL3mu_numOfLostHits"]->Fill(tk->numberOfLostHits());
	  hists_["hltL3mu_numValidMuonHits"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  
	  if (std::abs(l3mu.eta())<0.9) {
	    hists_["hltL3mu_numValidHits_barrel"]->Fill(tk->numberOfValidHits());
	    hists_["hltL3mu_normalizedChi2_barrel"]->Fill(tk->normalizedChi2());
	    hists_["hltL3mu_numOfLostHits_barrel"]->Fill(tk->numberOfLostHits());
	    hists_["hltL3mu_numValidMuonHits_barrel"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  }
	  else {
	    hists_["hltL3mu_numValidHits_endcap"]->Fill(tk->numberOfValidHits());
	    hists_["hltL3mu_normalizedChi2_endcap"]->Fill(tk->normalizedChi2());
	    hists_["hltL3mu_numOfLostHits_endcap"]->Fill(tk->numberOfLostHits());
	    hists_["hltL3mu_numValidMuonHits_endcap"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  }
	  hists_["hltL3mu_NumSharedHits"]->Fill(SharedHitMatrix[i][j]);
	  hists_["hltL3mu_FracSharedHits"]->Fill(SharedHitMatrix[i][j]/recoMu.innerTrack()->numberOfValidHits());
	}
      }
    }
    // RECO MUO param
    hists_["RecoMu_numValidHits"]->Fill(recoMu.innerTrack()->numberOfValidHits());
    hists_["RecoMu_normalizedChi2"]->Fill(recoMu.innerTrack()->normalizedChi2());
    hists_["RecoMu_numOfLostHits"]->Fill(recoMu.innerTrack()->numberOfLostHits());
    hists_["RecoMu_numValidMuonHits"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    if (std::abs(recoMu.eta())<0.9) {
      hists_["RecoMu_numValidHits_barrel"]->Fill(recoMu.innerTrack()->numberOfValidHits());
      hists_["RecoMu_normalizedChi2_barrel"]->Fill(recoMu.innerTrack()->normalizedChi2());
      hists_["RecoMu_numOfLostHits_barrel"]->Fill(recoMu.innerTrack()->numberOfLostHits());
      hists_["RecoMu_numValidMuonHits_barrel"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    }
    else{ 
      hists_["RecoMu_numValidHits_endcap"]->Fill(recoMu.innerTrack()->numberOfValidHits());
      hists_["RecoMu_normalizedChi2_endcap"]->Fill(recoMu.innerTrack()->normalizedChi2());
      hists_["RecoMu_numOfLostHits_endcap"]->Fill(recoMu.innerTrack()->numberOfLostHits());
      hists_["RecoMu_numValidMuonHits_endcap"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    }
  }
  
  
  // Run through the matrix n1 times to make sure we've found all matches.
  for (size_t k = 0; k < n1; k++) {
    size_t i_min = -1;
    size_t j_min = -1;
    double maxSharedFrac = minSharedFrac;
    // find the smallest deltaR
    for (size_t i = 0; i < n1; i++)
      for (size_t j = 0; j < n2; j++)
	if (SharedHitMatrix[i][j] > maxSharedFrac) {
	  i_min = i;
	  j_min = j;
	  maxSharedFrac = SharedHitMatrix[i][j];
	}

    // If a match has been made, save it and make those candidates unavailable.
    if (maxSharedFrac > minSharedFrac) {
      result[i_min] = j_min;
      SharedHitMatrix[i_min] = vector<double>(n2, NOMATCHITS);
      for (size_t i = 0; i < n1; i++)
	SharedHitMatrix[i][j_min] = NOMATCHITS;
    }
  }
  return result;
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
  hists_["hltL1_phi"]    = outfile_->make<TH1F>("hltL1_phi", "HLT (L1) #phi;#phi of L1 object", 13, phi_bins);
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

  hists_["hltL2_pt"]            = outfile_->make<TH1F>("hltL2_pt",  "HLT (L2) p_{T}; p_{T} of L2 object", 11,  pt_bins );
  hists_["hltL2_eta"]           = outfile_->make<TH1F>("hltL2_eta", "HLT (L2) #eta; #eta of L2 object", 15, eta_bins );
  hists_["hltL2_phi"]           = outfile_->make<TH1F>("hltL2_phi", "HLT (L2) #phi;#phi of L2 object", 13, phi_bins);
  hists_["hltL2_DeltaR"] = outfile_->make<TH1F>("hltL2_DeltaR", "HLT (L2) #Delta R; #Delta wrt L2 object", 15, 0., 1.);
  hists_["hltL2_resEta"] = outfile_->make<TH1F>("hltL2_resEta", "L2 Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.05,   0.05);
  hists_["hltL2_resPhi"] = outfile_->make<TH1F>("hltL2_resPhi", "L2 Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.05,   0.05);
  hists_["hltL2_resPt"]  = outfile_->make<TH1F>("hltL2_resPt",  "L2 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.5,   0.5);

  hists_["hltL3_pt"]     = outfile_->make<TH1F>("hltL3_pt",  "HLT (L3) p_{T}; p_{T} of L3 object", 11,  pt_bins );
  hists_["hltL3_eta"]    = outfile_->make<TH1F>("hltL3_eta", "HLT (L3) #eta; #eta of L3 object", 15, eta_bins );
  hists_["hltL3_phi"]    = outfile_->make<TH1F>("hltL3_phi", "HLT (L3) #phi;#phi of L3 object", 13, phi_bins);
  hists_["hltL3_DeltaR"] = outfile_->make<TH1F>("hltL3_DeltaR", "HLT (L3) #Delta R; #Delta wrt L3 object", 15, 0., 1.);
  hists_["hltL3_resEta"] = outfile_->make<TH1F>("hltL3_resEta", "L3 Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPhi"] = outfile_->make<TH1F>("hltL3_resPhi", "L3 Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPt"]  = outfile_->make<TH1F>("hltL3_resPt",  "L3 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.30,   0.30);
  
  for (size_t j = 0; j < 2; j++) {
    string suffix = EFFICIENCY_SUFFIXES[j];
    hists_["effL3L2Eta_"+suffix ] = outfile_->make<TH1F>(("effL3L2Eta_"+suffix).c_str(), "Efficiency;#eta",  15, eta_bins );
    hists_["effL3L2Phi_"+suffix ] = outfile_->make<TH1F>(("effL3L2Phi_"+suffix).c_str(), "Efficiency;#phi",  13, phi_bins);
    hists_["effL3L2Pt_"+suffix  ] = outfile_->make<TH1F>(("effL3L2Pt_"+suffix).c_str(),  "Efficiency;p_{T}", 11,  pt_bins );
    hists_["effL3L2NPU_"+suffix ] = outfile_->make<TH1F>(("effL3L2NPU_"+suffix).c_str(), "Efficiency;NPU",   75,  0,   75.);
    hists_["effL3L2NPV_"+suffix ] = outfile_->make<TH1F>(("effL3L2NPV_"+suffix).c_str(), "Efficiency;NPV",   75,  0,   75.);
    hists_["effL3L2dxy_"+suffix ] = outfile_->make<TH1F>(("effL3L2dxy_"+suffix).c_str(), "Efficiency;dxy",   50,  0., 0.2);
    hists_["effL3L2dz_"+suffix ]  = outfile_->make<TH1F>(("effL3L2dz_"+suffix).c_str() , "Efficiency;dz" ,   50,  0., 0.5);
    hists_["effL3L2Chi2_"+suffix ]= outfile_->make<TH1F>(("effL3L2Chi2_"+suffix).c_str() , "Efficiency;#chi^2/n.d.o.f.",   50,  0., 20.);
    
    hists_["eff2L3L2Eta_"+suffix ] = outfile_->make<TH1F>(("eff2L3L2Eta_"+suffix).c_str(), "Efficiency;#eta",  15, eta_bins );
    hists_["eff2L3L2Phi_"+suffix ] = outfile_->make<TH1F>(("eff2L3L2Phi_"+suffix).c_str(), "Efficiency;#phi",  13, phi_bins);
    hists_["eff2L3L2Pt_"+suffix  ] = outfile_->make<TH1F>(("eff2L3L2Pt_"+suffix).c_str(),  "Efficiency;p_{T}", 11,  pt_bins );
    hists_["eff2L3L2NPU_"+suffix ] = outfile_->make<TH1F>(("eff2L3L2NPU_"+suffix).c_str(), "Efficiency;NPU",   75,  0,   75.);
    hists_["eff2L3L2NPV_"+suffix ] = outfile_->make<TH1F>(("eff2L3L2NPV_"+suffix).c_str(), "Efficiency;NPV",   75,  0,   75.);

    hists_["effL3Eta_"+suffix ] = outfile_->make<TH1F>(("effL3Eta_"+suffix).c_str(), "Efficiency;#eta",  15, eta_bins );
    hists_["effL3Phi_"+suffix ] = outfile_->make<TH1F>(("effL3Phi_"+suffix).c_str(), "Efficiency;#phi",  13, phi_bins);
    hists_["effL3Pt_"+suffix  ] = outfile_->make<TH1F>(("effL3Pt_"+suffix).c_str(),  "Efficiency;p_{T}", 11,  pt_bins );
    hists_["effL3NPU_"+suffix ] = outfile_->make<TH1F>(("effL3NPU_"+suffix).c_str(), "Efficiency;NPU",   75,  0,   75.);
    hists_["effL3NPV_"+suffix ] = outfile_->make<TH1F>(("effL3NPV_"+suffix).c_str(), "Efficiency;NPV",   75,  0,   75.);
    hists_["effL3dxy_"+suffix ] = outfile_->make<TH1F>(("effL3dxy_"+suffix).c_str(), "Efficiency;dxy",   50,  0., 0.2);
    hists_["effL3dz_"+suffix ]  = outfile_->make<TH1F>(("effL3dz_"+suffix).c_str(),  "Efficiency;dz" ,   50,  0., 0.5);
    hists_["effL3Chi2_"+suffix ]= outfile_->make<TH1F>(("effL3Chi2_"+suffix).c_str() , "Efficiency;#chi^2/n.d.o.f.",   50,  0., 20.);
    
    hists_["eff2L3Eta_"+suffix ] = outfile_->make<TH1F>(("eff2L3Eta_"+suffix).c_str(), "Efficiency;#eta",  15, eta_bins );
    hists_["eff2L3Phi_"+suffix ] = outfile_->make<TH1F>(("eff2L3Phi_"+suffix).c_str(), "Efficiency;#phi",  13, phi_bins);
    hists_["eff2L3Pt_"+suffix  ] = outfile_->make<TH1F>(("eff2L3Pt_"+suffix).c_str(),  "Efficiency;p_{T}", 11,  pt_bins );
    hists_["eff2L3NPU_"+suffix ] = outfile_->make<TH1F>(("eff2L3NPU_"+suffix).c_str(), "Efficiency;NPU",   75,  0,   75.);
    hists_["eff2L3NPV_"+suffix ] = outfile_->make<TH1F>(("eff2L3NPV_"+suffix).c_str(), "Efficiency;NPV",   75,  0,   75.);
    
    
    //// SHARED HITS... 
    hists_["effSharedL3L2Eta_"+suffix ] = outfile_->make<TH1F>(("effSharedL3L2Eta_"+suffix).c_str(), "Efficiency (Shared Hits);#eta",  15, eta_bins );
    hists_["effSharedL3L2Phi_"+suffix ] = outfile_->make<TH1F>(("effSharedL3L2Phi_"+suffix).c_str(), "Efficiency (Shared Hits);#phi",  13, phi_bins);
    hists_["effSharedL3L2Pt_"+suffix  ] = outfile_->make<TH1F>(("effSharedL3L2Pt_"+suffix).c_str(),  "Efficiency (Shared Hits);p_{T}", 11,  pt_bins );
    hists_["effSharedL3L2NPU_"+suffix ] = outfile_->make<TH1F>(("effSharedL3L2NPU_"+suffix).c_str(), "Efficiency (Shared Hits);NPU",   75,  0,   75.);
    hists_["effSharedL3L2NPV_"+suffix ] = outfile_->make<TH1F>(("effSharedL3L2NPV_"+suffix).c_str(), "Efficiency (Shared Hits);NPV",   75,  0,   75.);
    hists_["effSharedL3L2dxy_"+suffix ] = outfile_->make<TH1F>(("effSharedL3L2dxy_"+suffix).c_str(), "Efficiency (Shared Hits);dxy",   50,  0., 0.2);
    hists_["effSharedL3L2dz_"+suffix ]  = outfile_->make<TH1F>(("effSharedL3L2dz_"+suffix).c_str() , "Efficiency (Shared Hits);dz" ,   50,  0., 0.5);
    hists_["effSharedL3L2Chi2_"+suffix ]= outfile_->make<TH1F>(("effSharedL3L2Chi2_"+suffix).c_str() ,"Efficiency (Shared Hits);#chi^2/n.d.o.f.",   50,  0., 20.);
    
    hists_["effSharedL3Eta_"+suffix ] = outfile_->make<TH1F>(("effSharedL3Eta_"+suffix).c_str(), "Efficiency (Shared Hits);#eta",  15, eta_bins );
    hists_["effSharedL3Phi_"+suffix ] = outfile_->make<TH1F>(("effSharedL3Phi_"+suffix).c_str(), "Efficiency (Shared Hits);#phi",  13, phi_bins);
    hists_["effSharedL3Pt_"+suffix  ] = outfile_->make<TH1F>(("effSharedL3Pt_"+suffix).c_str(),  "Efficiency (Shared Hits);p_{T}", 11,  pt_bins );
    hists_["effSharedL3NPU_"+suffix ] = outfile_->make<TH1F>(("effSharedL3NPU_"+suffix).c_str(), "Efficiency (Shared Hits);NPU",   75,  0,   75.);
    hists_["effSharedL3NPV_"+suffix ] = outfile_->make<TH1F>(("effSharedL3NPV_"+suffix).c_str(), "Efficiency (Shared Hits);NPV",   75,  0,   75.);
    hists_["effSharedL3dxy_"+suffix ] = outfile_->make<TH1F>(("effSharedL3dxy_"+suffix).c_str(), "Efficiency (Shared Hits);dxy",   50,  0., 0.2);
    hists_["effSharedL3dz_"+suffix ]  = outfile_->make<TH1F>(("effSharedL3dz_"+suffix).c_str(),  "Efficiency (Shared Hits);dz" ,   50,  0., 0.5);
    hists_["effSharedL3Chi2_"+suffix ]= outfile_->make<TH1F>(("effSharedL3Chi2_"+suffix).c_str() , "Efficiency (Shared Hits);#chi^2/n.d.o.f.",   50,  0., 20.);    
  }
  

  //// COUNTERS
  hists_["hlt_NumL1Match" ] = outfile_->make<TH1F>("hlt_NumL1Match","Number of L1 Gen Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL2Match" ] = outfile_->make<TH1F>("hlt_NumL2Match","Number of L2 Gen Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL3Match" ] = outfile_->make<TH1F>("hlt_NumL3Match","Number of L3 Gen Matched", 5, -0.5, 4.5);

  hists_["hlt_NumL1" ] = outfile_->make<TH1F>("hlt_NumL1","Number of L1 Found", 5, -0.5, 4.5);
  hists_["hlt_NumL2" ] = outfile_->make<TH1F>("hlt_NumL2","Number of L2 Found", 5, -0.5, 4.5);
  hists_["hlt_NumL3" ] = outfile_->make<TH1F>("hlt_NumL3","Number of L3 Found", 5, -0.5, 4.5);

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
  hists_["hlt_normalizedChi2"]     = outfile_->make<TH1F>("hlt_normalizedChi2", "Normalised Chi2 of Iter0+Iter2 ", 100, 0., 20); 
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
  hists_["hlt_noL3_trackphi"]      = outfile_->make<TH1F>("hlt_noL3_trackphi" ,"Track (L3) #phi; #phi of L3 track",  13, phi_bins);
  hists_["hlt_noL3_DeltaR"]        = outfile_->make<TH1F>("hlt_noL3_DeltaR",   "Iter L3 #Delta R; #Delta R (tk,L2)", 15, 0., 0.2);
  hists_["hlt_noL3_DeltaEta"]      = outfile_->make<TH1F>("hlt_noL3_DeltaEta", "Iter L3 #Delta #eta;#eta^{tk}-#eta^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPhi"]      = outfile_->make<TH1F>("hlt_noL3_DeltaPhi", "Iter L3 #Delta #phi;#phi^{tk}-#phi^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPt"]       = outfile_->make<TH1F>("hlt_noL3_DeltaPt",  "Iter L3 #Delta p_{T};p_{T}^{tk}-p_{T}^{L2}", 60,  -0.30,   0.30);
  
  hists_["hltL3mu_dR2withLink"]      = outfile_->make<TH1F>("hltL3mu_dR2withLink"     , "#Delta R2 with Link", 100, 0, 0.02);
  hists_["hltL3mu_dPtwithLink"]      = outfile_->make<TH1F>("hltL3mu_dPtwithLink"     , "#Delta pt with Link", 100, 0, 0.05);
  hists_["hltL3mu_pt"]               = outfile_->make<TH1F>("hltL3mu_pt"              , "#p_{T} of the L3 Muon; p_{T} L3 muon", 11,  pt_bins );
  hists_["hltL3mu_eta"]		     = outfile_->make<TH1F>("hltL3mu_eta"             , "#eta of the L3 Muon; #eta L3 muons", 15, eta_bins ); 
  hists_["hltL3mu_dr"]		     = outfile_->make<TH1F>("hltL3mu_dr"              , "Dr w.r.t. BS", 50, 0., 0.2 ); 
  hists_["hltL3mu_dz"]		     = outfile_->make<TH1F>("hltL3mu_dz"              , "Dz w.r.t. BS", 50, 0., 0.5); 
  hists_["hltL3mu_dxySig"]	     = outfile_->make<TH1F>("hltL3mu_dxySig"          , "Dxy Significance ", 50, 0., 15.); 
  hists_["hltL3mu_dxy"]		     = outfile_->make<TH1F>("hltL3mu_dxy"             , "Dxy w.r.t. BS", 50, 0., 0.2); 
  hists_["hltL3mu_numValidHits"]     = outfile_->make<TH1F>("hltL3mu_numValidHits"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits"]     = outfile_->make<TH1F>("RecoMu_numValidHits"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_numValidHits_barrel"]     = outfile_->make<TH1F>("hltL3mu_numValidHits_barrel"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits_barrel"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits_barrel"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2_barrel"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2_barrel"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits_barrel"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits_barrel", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits_barrel"]     = outfile_->make<TH1F>("RecoMu_numValidHits_barrel"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits_barrel"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits_barrel"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2_barrel"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2_barrel"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits_barrel"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits_barrel", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_numValidHits_endcap"]     = outfile_->make<TH1F>("hltL3mu_numValidHits_endcap"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits_endcap"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits_endcap"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2_endcap"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2_endcap"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits_endcap"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits_endcap", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits_endcap"]     = outfile_->make<TH1F>("RecoMu_numValidHits_endcap"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits_endcap"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits_endcap"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2_endcap"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2_endcap"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits_endcap"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits_endcap", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_DiffHits"]               = outfile_->make<TH1F>("hltL3mu_DiffHits"    , "Diff Shared Hits L3/RECO Muon", 100,  0., 5. ); 
  hists_["hltL3mu_NumSharedHits"]          = outfile_->make<TH1F>("hltL3mu_NumSharedHits"    , "N Shared Hits L3/RECO Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_FracSharedHits"]         = outfile_->make<TH1F>("hltL3mu_FracValidHits"    , "Fraction Shared Hits L3/RECO Muon", 100,  0., 1.05 ); 
 
  hists_["hlt_OI_NumSeeds"]               = outfile_->make<TH1F>("hlt_OI_NumSeeds","Number of Seeds (OI)", 50, -0.5, 49.5);
  hists_["hlt_OI_NumSeedsVsEta"]          = outfile_->make<TH2F>("hlt_OI_NumSeedsVsEta","Number of Seeds (OI) vs #eta", 15, eta_bins, 50, -0.5, 49.5);
  hists_["hlt_OI_NumberOfRecHitsPerSeed"] = outfile_->make<TH1F>("hlt_OI_NumberOfRecHitsPerSeed","Number Hits per Seed (OI)", 50, -0.5, 49.5);
  hists_["hlt_OI_seedPtErr"]              = outfile_->make<TH1F>("hlt_OI_seedPtErr", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErr_barrel"]       = outfile_->make<TH1F>("hlt_OI_seedPtErr_barrel", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErr_endcap"]       = outfile_->make<TH1F>("hlt_OI_seedPtErr_endcap", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErrVsEta"]         = outfile_->make<TH2F>("hlt_OI_seedPrErrVsEta", "Seed Pt Error vs Eta (OI)", 15, eta_bins, 50, 0., 50. );
  hists_["hlt_OI_seedPhiErr"]             = outfile_->make<TH1F>("hlt_OI_seedPhiErr","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedPhiErr_barrel"]      = outfile_->make<TH1F>("hlt_OI_seedPhiErr_barrel","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedPhiErr_endcap"]      = outfile_->make<TH1F>("hlt_OI_seedPhiErr_endcap","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr"]             = outfile_->make<TH1F>("hlt_OI_seedEtaErr","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr_barrel"]      = outfile_->make<TH1F>("hlt_OI_seedEtaErr_barrel","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr_endcap"]      = outfile_->make<TH1F>("hlt_OI_seedEtaErr_endcap","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_hitPt"]                  = outfile_->make<TH1F>("hlt_OI_hitPt","Hit p_{T} (OI);p_{T} hit", 11,  pt_bins );
  hists_["hlt_OI_hitPhi"]                 = outfile_->make<TH1F>("hlt_OI_hitPhi","Hit #phi (OI);#phi hit", 13, phi_bins);
  hists_["hlt_OI_hitEta"]                 = outfile_->make<TH1F>("hlt_OI_hitEta","Hit #eta (OI);#eta hit", 15,  eta_bins ) ;
  hists_["hlt_OI_hitx"]                   = outfile_->make<TH1F>("hlt_OI_hitx","Hit x (OI);x hit", 200, -100, 100);
  hists_["hlt_OI_hity"]                   = outfile_->make<TH1F>("hlt_OI_hity","Hit y (OI);y hit", 200, -100, 100);
  hists_["hlt_OI_hitz"]                   = outfile_->make<TH1F>("hlt_OI_hitz","Hit z (OI);z hit", 600, -300, 300);
  hists_["hlt_OI_hitdx"]                  = outfile_->make<TH1F>("hlt_OI_hitdx","Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_hitdy"]                  = outfile_->make<TH1F>("hlt_OI_hitdy","Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);

  /// TOB hits
  hists_["hlt_OI_TOBhitdx"]               = outfile_->make<TH1F>("hlt_OI_TOBhitdx",      "TOB Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBhitdy"]               = outfile_->make<TH1F>("hlt_OI_TOBhitdy",      "TOB Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);
  hists_["hlt_OI_TOBmonohitdx"]           = outfile_->make<TH1F>("hlt_OI_TOBmonohitdx",  "TOB Mono Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBmonohitdy"]           = outfile_->make<TH1F>("hlt_OI_TOBmonohitdy",  "TOB Mono Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);
  hists_["hlt_OI_TOBstereohitdx"]         = outfile_->make<TH1F>("hlt_OI_TOBstereohitdx","TOB Stereo Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBstereohitdy"]         = outfile_->make<TH1F>("hlt_OI_TOBstereohitdy","TOB Stereo Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);

  /// TEC hits
  hists_["hlt_OI_TEChitdx"]               = outfile_->make<TH1F>("hlt_OI_TEChitdx",      "TEC Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TEChitdy"]               = outfile_->make<TH1F>("hlt_OI_TEChitdy",      "TEC Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);
  hists_["hlt_OI_TECmonohitdx"]           = outfile_->make<TH1F>("hlt_OI_TECmonohitdx",  "TEC Mono Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TECmonohitdy"]           = outfile_->make<TH1F>("hlt_OI_TECmonohitdy",  "TEC Mono Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);
  hists_["hlt_OI_TECstereohitdx"]         = outfile_->make<TH1F>("hlt_OI_TECstereohitdx","TEC Stereo Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TECstereohitdy"]         = outfile_->make<TH1F>("hlt_OI_TECstereohitdy","TEC Stereo Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);

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
