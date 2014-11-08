// -*- C++ -*-
//
// Package:    MetTesting/MetFilter
// Class:      MetFilter
// 
/**\class MetFilter MetFilter.cc MetTesting/MetFilter/plugins/MetFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Young Ho Shin
//         Created:  Wed, 15 Oct 2014 07:43:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// For creating/writing histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// DataFormats
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CorrMETData.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
// class declaration
//

class MetFilter : public edm::EDFilter {
   public:
      explicit MetFilter(const edm::ParameterSet&);
      ~MetFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
			std::string mOutputFileName;

			// Tokens for consumed products
			edm::EDGetTokenT<reco::PFJetCollection> mJetsToken, mChsJetsToken;
			edm::EDGetTokenT<reco::VertexCollection> mVerticesToken;
};

//
// constants, enums and typedefs
//
typedef std::vector<bool> boolVector;
typedef std::vector<int> intVector;
typedef std::vector<float> floatVector;
typedef std::vector<double> doubleVector;

//
// static data member definitions
//

//
// constructors and destructor
//
MetFilter::MetFilter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed

	// Register consumed products
	// Using InputTag from iConfig is proper way to do it.
	// For now just set sensible defaults.
	// edm::InputTag jets_     = iConfig.getParameter<edm::InputTag> ("jetCollection"    ) ;
	// edm::InputTag chsJets_  = iConfig.getParameter<edm::InputTag> ("chsJetCollection" ) ;
	// edm::InputTag vertices_ = iConfig.getParameter<edm::InputTag> ("vertexCollection" ) ;
	edm::InputTag jets_("ak4PFJets");
	edm::InputTag chsJets_("ak4PFJetsCHS");
	edm::InputTag vertices_("offlinePrimaryVertices");

	mJetsToken     = consumes<reco::PFJetCollection>  ( jets_     ) ;
	mChsJetsToken  = consumes<reco::PFJetCollection>  ( chsJets_  ) ;
	mVerticesToken = consumes<reco::VertexCollection> ( vertices_ ) ;

	// Register products
  /*
	 * // MET:
	 * // pfMet, pfMetT1, pfMetT1 w/ CHS Jets
	 * produces< double > ( "pfMetPt"         ).setBranchAlias ( "pfMetPt"         ) ;
	 * produces< double > ( "pfMetPhi"        ).setBranchAlias ( "pfMetPhi"        ) ;
	 * produces< double > ( "pfMetSumEt"      ).setBranchAlias ( "pfMetSumEt"      ) ;
	 * produces< double > ( "pfMetT1Pt"       ).setBranchAlias ( "pfMetT1Pt"       ) ;
	 * produces< double > ( "pfMetT1Phi"      ).setBranchAlias ( "pfMetT1Phi"      ) ;
	 * produces< double > ( "pfMetT1SumEt"    ).setBranchAlias ( "pfMetT1SumEt"    ) ;
	 * produces< double > ( "pfMetT1CHSPt"    ).setBranchAlias ( "pfMetT1CHSPt"    ) ;
	 * produces< double > ( "pfMetT1CHSPhi"   ).setBranchAlias ( "pfMetT1CHSPhi"   ) ;
	 * produces< double > ( "pfMetT1CHSSumEt" ).setBranchAlias ( "pfMetT1CHSSumEt" ) ;
	 * // Muons:
	 * // cleanPatMuons
	 * produces< doubleVector > ("muonPts"  ) .setBranchAlias("muonPts"  ) ;
	 * produces< doubleVector > ("muonEtas" ) .setBranchAlias("muonEtas" ) ;
	 * produces< doubleVector > ("muonPhis" ) .setBranchAlias("muonPhis" ) ;
	 * // Jets:
	 * // Jets, CHSJets
	 * produces< doubleVector > ( "JetPts"     ).setBranchAlias ( "JetPts"     ) ;
	 * produces< doubleVector > ( "JetEtas"    ).setBranchAlias ( "JetEtas"    ) ;
	 * produces< doubleVector > ( "JetPhis"    ).setBranchAlias ( "JetPhis"    ) ;
	 * produces< doubleVector > ( "ChsJetPts"  ).setBranchAlias ( "ChsJetPts"  ) ;
	 * produces< doubleVector > ( "ChsJetEtas" ).setBranchAlias ( "ChsJetEtas" ) ;
	 * produces< doubleVector > ( "ChsJetPhis" ).setBranchAlias ( "ChsJetPhis" ) ;
   */
	// Vertices:
	// Number of good offlinePrimaryVertices
	produces< int > ("nVtx").setBranchAlias("nVtx");
	produces< int > ("nGoodVtx").setBranchAlias("nGoodVtx");


}


MetFilter::~MetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  /*
	 * ////////////////////////////////////////////////
	 * // MET
	 * ////////////////////////////////////////////////
	 * std::auto_ptr< double > pfMetPt         ( new double);
	 * std::auto_ptr< double > pfMetPhi        ( new double);
	 * std::auto_ptr< double > pfMetSumEt      ( new double);
	 * std::auto_ptr< double > pfMetT1Pt       ( new double);
	 * std::auto_ptr< double > pfMetT1Phi      ( new double);
	 * std::auto_ptr< double > pfMetT1SumEt    ( new double);
	 * std::auto_ptr< double > pfMetT1CHSPt    ( new double);
	 * std::auto_ptr< double > pfMetT1CHSPhi   ( new double);
	 * std::auto_ptr< double > pfMetT1CHSSumEt ( new double);

	 * Handle<reco::PFMETCollection> pfMetHandle;
   * iEvent.getByLabel("pfMet", pfMetHandle);
	 * Handle<reco::PFMETCollection> pfMetT1Handle;
   * iEvent.getByLabel("pfMetT1", pfMetT1Handle);
	 * Handle<reco::PFMETCollection> pfMetT1CHSHandle;
   * iEvent.getByLabel("pfMetT1CHS", pfMetT1CHSHandle);

	 * *pfMetPt         = pfMetHandle.       product()->front().  pt();
	 * *pfMetPhi        = pfMetHandle.       product()->front().  phi();
	 * *pfMetSumEt      = pfMetHandle.       product()->front().  sumEt();
	 * *pfMetT1Pt       = pfMetT1Handle.     product()->front().  pt();
	 * *pfMetT1Phi      = pfMetT1Handle.     product()->front().  phi();
	 * *pfMetT1SumEt    = pfMetT1Handle.     product()->front().  sumEt();
	 * *pfMetT1CHSPt    = pfMetT1CHSHandle.  product()->front().  pt();
	 * *pfMetT1CHSPhi   = pfMetT1CHSHandle.  product()->front().  phi();
	 * *pfMetT1CHSSumEt = pfMetT1CHSHandle.  product()->front().  sumEt();

	 * iEvent.put (pfMetPt         , "pfMetPt"         ) ;
	 * iEvent.put (pfMetPhi        , "pfMetPhi"        ) ;
	 * iEvent.put (pfMetSumEt      , "pfMetSumEt"      ) ;
	 * iEvent.put (pfMetT1Pt       , "pfMetT1Pt"       ) ;
	 * iEvent.put (pfMetT1Phi      , "pfMetT1Phi"      ) ;
	 * iEvent.put (pfMetT1SumEt    , "pfMetT1SumEt"    ) ;
	 * iEvent.put (pfMetT1CHSPt    , "pfMetT1CHSPt"    ) ;
	 * iEvent.put (pfMetT1CHSPhi   , "pfMetT1CHSPhi"   ) ;
	 * iEvent.put (pfMetT1CHSSumEt , "pfMetT1CHSSumEt" ) ;

	 * ////////////////////////////////////////////////
	 * // Muons
	 * ////////////////////////////////////////////////
	 * std::auto_ptr< doubleVector > muonPts  ( new doubleVector);
	 * std::auto_ptr< doubleVector > muonEtas ( new doubleVector);
	 * std::auto_ptr< doubleVector > muonPhis ( new doubleVector);

	 * Handle<pat::MuonCollection> muonHandle;
   * iEvent.getByLabel("cleanPatMuons", muonHandle);

	 * muonPts->reserve(muonHandle->size());
	 * for (pat::MuonCollection::const_iterator it = muonHandle->begin(); it != muonHandle->end(); it++) {
	 * 	pat::Muon currentObj = *it;

	 * 	double pt  = currentObj.   pt();
	 * 	double eta = currentObj.  eta();
	 * 	double phi = currentObj.  phi();
	 * 	muonPts->push_back(pt);
	 * 	muonEtas->push_back(eta);
	 * 	muonPhis->push_back(phi);
	 * }

	 * iEvent.put (muonPts    , "muonPts"    );

	 * ////////////////////////////////////////////////
	 * // Jets
	 * ////////////////////////////////////////////////
	 * std::auto_ptr< doubleVector > jetPts             ( new doubleVector);
	 * std::auto_ptr< doubleVector > jetEtas            ( new doubleVector);
	 * std::auto_ptr< doubleVector > jetPhis            ( new doubleVector);
	 * std::auto_ptr< doubleVector > chsJetPts          ( new doubleVector);
	 * std::auto_ptr< doubleVector > chsJetEtas         ( new doubleVector);
	 * std::auto_ptr< doubleVector > chsJetPhis         ( new doubleVector);

	 * Handle<reco::PFJetCollection> jetCollection;
	 * Handle<reco::PFJetCollection> chsJetCollection;
	 * iEvent.getByToken(mJetsToken, jetCollection);
	 * iEvent.getByToken(mChsJetsToken, chsJetCollection);

	 * jetPts->reserve(jetCollection->size());
	 * jetEtas->reserve(jetCollection->size());
	 * jetPhis->reserve(jetCollection->size());
	 * for (reco::PFJetCollection::const_iterator it = jetCollection->begin(); it != jetCollection->end(); it++) {
	 * 	reco::PFJet currentObj = *it;
	 * 	double pt  = currentObj.   pt();
	 * 	double eta = currentObj.  eta();
	 * 	double phi = currentObj.  phi();
	 * 	jetPts->push_back(pt);
	 * 	jetEtas->push_back(eta);
	 * 	jetPhis->push_back(phi);
	 * }
	 * iEvent.put (jetPts,   "JetPts");
	 * iEvent.put (jetEtas,  "JetEtas");
	 * iEvent.put (jetPhis,  "JetPhis");

	 * chsJetPts->reserve(chsJetCollection->size());
	 * chsJetEtas->reserve(chsJetCollection->size());
	 * chsJetPhis->reserve(chsJetCollection->size());
	 * for (reco::PFJetCollection::const_iterator it = chsJetCollection->begin(); it != chsJetCollection->end(); it++) {
	 * 	reco::PFJet currentObj = *it;
	 * 	double pt  = currentObj.   pt();
	 * 	double eta = currentObj.  eta();
	 * 	double phi = currentObj.  phi();
	 * 	chsJetPts->push_back(pt);
	 * 	chsJetEtas->push_back(eta);
	 * 	chsJetPhis->push_back(phi);
	 * }
	 * iEvent.put (chsJetPts,   "ChsJetPts");
	 * iEvent.put (chsJetEtas,  "ChsJetEtas");
	 * iEvent.put (chsJetPhis,  "ChsJetPhis");
   */


	////////////////////////////////////////////////
	// Vertices
	////////////////////////////////////////////////
	std::auto_ptr< int > nVtx    ( new int);
	std::auto_ptr< int > nGoodVtx ( new int);
	*nVtx = 0;
	*nGoodVtx = 0;
	// int nVtx = 0;
	// int nAllVtx = 0;
	Handle<reco::VertexCollection> vertexCollection;
	iEvent.getByToken(mVerticesToken, vertexCollection);

	for (auto it = vertexCollection->begin(); it != vertexCollection->end(); it++) {
		auto currentVertex = *it;
		bool vertexIsGood = (!currentVertex.isFake()) && (currentVertex.ndof()>4) && (fabs(currentVertex.z()<=24.0)) && (currentVertex.position().Rho()<=2.0);
		// std::cout << "!currentVertex.isFake() ,  currentVertex.ndof()>4 ,  fabs(currentVertex.z()<=24.0) ,  currentVertex.position().Rho()<=2.0; " 
		// 	<< (!currentVertex.isFake()) << ", " << (currentVertex.ndof()>4) << ", " << (fabs(currentVertex.z()<=24.0)) << ", " << (currentVertex.position().Rho()<=2.0) << std::endl;
		// std::cout << "vertexIsGood " << vertexIsGood << std::endl;
		if (vertexIsGood) (*nGoodVtx)++;
		(*nVtx)++;
	}
	// std::cout << "nGoodVtx, nVtx:\t" << *nGoodVtx << ", " << *nVtx << "\n";

	iEvent.put (nVtx, "nVtx");
	iEvent.put (nGoodVtx, "nGoodVtx");

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MetFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MetFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MetFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MetFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MetFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MetFilter);
