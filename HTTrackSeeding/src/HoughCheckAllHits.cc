// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      HoughCheckAllHits
// 
/**\class HoughCheckAllHits HoughCheckAllHits.cc HoughTest/HoughTransChecks/src/HoughCheckAllHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Thu, May 16, 2013
// $Id: HoughCheckAllHits.cc
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/PerigeeConversions.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"

#include "TString.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheckAllHits : public edm::EDAnalyzer {
  public:
    explicit HoughCheckAllHits(const edm::ParameterSet&);
    ~HoughCheckAllHits();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  public:
    virtual void  beginJob ( );
    virtual void  analyze  (const edm::Event&, const edm::EventSetup&);
    virtual void  endJob   ( );
    virtual void  beginRun             ( edm::Run const&, edm::EventSetup const& );
    virtual void  endRun               ( edm::Run const&, edm::EventSetup const& );
    virtual void  beginLuminosityBlock ( edm::LuminosityBlock const&, edm::EventSetup const& );
    virtual void  endLuminosityBlock   ( edm::LuminosityBlock const&, edm::EventSetup const& );
  private:
    edm::InputTag               trackTags_;
    const vector<unsigned int>  algoSel_;
    const unsigned int          verbosity_;
};


HoughCheckAllHits::HoughCheckAllHits(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getParameter<edm::InputTag>("tracks")),
  algoSel_  (iConfig.getParameter<vector<unsigned int> >("algoSel")),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0))
{ }


HoughCheckAllHits::~HoughCheckAllHits()
{ } 


void
HoughCheckAllHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using  reco::TrackCollection;
  Handle<TrackCollection>  tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int ntrk = 0;
  int mtrk = 0;
  for ( TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack, ++ntrk ) {
    // Get fitted track parameters in perigee convention
    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    GlobalPoint refPoint(itTrack->vx(), itTrack->vy(), itTrack->vz());
    GlobalVector pTrack(itTrack->px(), itTrack->py(), itTrack->pz());
    FreeTrajectoryState fts(refPoint, pTrack, TrackCharge(itTrack->charge()), theMF.product());
    TSCPBuilderNoMaterial tscpBuilder;
    TrajectoryStateClosestToPoint tsAtClosestApproach = tscpBuilder(fts, GlobalPoint(0.,0.,0.));
    PerigeeTrajectoryParameters perigeePars = tsAtClosestApproach.perigeeParameters();
    if (verbosity_ > 1)
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", Nhits=" << itTrack->recHitsSize()
	   << ", (vx,vy,vz)=(" << itTrack->vx() << "," << itTrack->vy() << "," << itTrack->vz() << ")"
           << ", doca=" << perigeePars.transverseImpactParameter()
           << ", kappa=" << perigeePars.transverseCurvature()
	   << ", phi=" << perigeePars.phi()
           << ", z0=" << perigeePars.longitudinalImpactParameter()
	   << ", theta=" << perigeePars.theta()
	   << ", algo=" << itTrack->algoName(itTrack->algo()).c_str() << endl;
    if (algoSel_.size() != 0) {  // empty vector means 'keep all tracks'
      bool validAlgo = false;
      for ( vector<unsigned int>::const_iterator itAlgo = algoSel_.begin(); itAlgo != algoSel_.end(); itAlgo++ )
        if (itTrack->algo() == *itAlgo) {
          validAlgo = true;
          break;
        }
      if (!validAlgo)  continue;
    }
    mtrk++;
  }
  //LogInfo ("HTTrackSeeding") << "Tracks f unds after trackfitter" << mtrk;  
  std::cout << "Tracks founds after trackfitter " << mtrk << std::endl;  
}


void 
HoughCheckAllHits::beginJob()
{ } 


void 
HoughCheckAllHits::endJob() 
{ }


void 
HoughCheckAllHits::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{ }


void 
HoughCheckAllHits::endRun(edm::Run const&, edm::EventSetup const&)
{ }


void 
HoughCheckAllHits::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{ }


void 
HoughCheckAllHits::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{ }


void
HoughCheckAllHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<vector<unsigned int> >("algoSel", {});
  desc.addUntracked<unsigned int>("verbosity", 0);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HoughCheckAllHits);
