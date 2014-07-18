// -*- C++ -*-
//
// Package:    HTProdRecHits
// Class:      HTProdRecHits
// 
/**\class HTProdRecHits HTProdRecHits.cc MLoVetere/HTProdRecHits/src/HTProdRecHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maurizio Lo Vetere,559 R-009,+41227675905,
//         Created:  Fri Nov 30 21:19:49 CET 2012
// $Id: HTProdRecHits.cc,v 1.1 2012/12/06 22:29:23 mlv Exp $
//
//

#include <iostream>

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "TH1.h"
#include "TFile.h"


class HTProdRecHits : public edm::EDAnalyzer {
  public:
    explicit      HTProdRecHits    ( const edm::ParameterSet& );
                 ~HTProdRecHits    ( );
    static  void  fillDescriptions ( edm::ConfigurationDescriptions& descriptions );
  private:
    virtual void  analyze  ( const edm::Event&, const edm::EventSetup& );
    virtual void  beginRun ( edm::Run const&  , edm::EventSetup const& );
    virtual void  beginJob ( );
    virtual void  endJob   ( );
  private:
    //edm::InputTag pixelClusterTags_;
    //edm::InputTag stripClusterTags_;
    edm::InputTag pixelRecHitTags_;
    edm::InputTag stripMatchedRecHitTags_;
    edm::InputTag stripRphiRecHitTags_;
    edm::InputTag stripRphiRecHitUnmatchedTags_;
    edm::InputTag stripStereoRecHitTags_;
    edm::InputTag stripStereoRecHitUnmatchedTags_;
    edm::ESHandle<TransientTrackingRecHitBuilder> builder_;
    std::string   builderName_;
    std::string   fileName_;
    TFile * file_;
    TH1F  * hNumPixelClusters_;
    TH1F  * hNumStripClusters_;
};


HTProdRecHits::HTProdRecHits(const edm::ParameterSet& iConfig)
  : //pixelClusterTags_              (iConfig.getParameter<edm::InputTag>("pixelClusters"              )),
    //stripClusterTags_              (iConfig.getParameter<edm::InputTag>("stripClusters"              )),
    pixelRecHitTags_               (iConfig.getParameter<edm::InputTag>("pixelRecHits"               )),
    stripMatchedRecHitTags_        (iConfig.getParameter<edm::InputTag>("stripMatchedRecHits"        )),
    stripRphiRecHitTags_           (iConfig.getParameter<edm::InputTag>("stripRphiRecHits"           )),
    stripRphiRecHitUnmatchedTags_  (iConfig.getParameter<edm::InputTag>("stripRphiRecHitsUnmatched"  )),
    stripStereoRecHitTags_         (iConfig.getParameter<edm::InputTag>("stripStereoRecHits"         )),
    stripStereoRecHitUnmatchedTags_(iConfig.getParameter<edm::InputTag>("stripStereoRecHitsUnmatched")),
    builderName_                   (iConfig.getParameter<std::string>  ("TTRHBuilder"                )),
    fileName_                      (iConfig.getUntrackedParameter<std::string>("hfile"))
{
}


HTProdRecHits::~HTProdRecHits()
{
}


void
HTProdRecHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /*
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelClustersHandle;
  iEvent.getByLabel( pixelClusterTags_, pixelClustersHandle );
  hNumPixelClusters_->Fill( pixelClustersHandle.product()->size() );
  std::cout << "====> number of pixel clusters " << pixelClustersHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator iter = pixelClustersHandle->begin();
        iter != pixelClustersHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > stripClustersHandle;
  iEvent.getByLabel( stripClusterTags_, stripClustersHandle );
  hNumStripClusters_->Fill( stripClustersHandle.product()->size() );
  std::cout << "====> number of strip clusters " << stripClustersHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripCluster>::const_iterator iter = stripClustersHandle->begin();
        iter != stripClustersHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  */
  edm::Handle<edmNew::DetSetVector<SiPixelRecHit> > pixelRecHitsHandle;
  iEvent.getByLabel( pixelRecHitTags_, pixelRecHitsHandle );
  std::cout << "====> number of pixel reco hits " << pixelRecHitsHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiPixelRecHit>::const_iterator iter = pixelRecHitsHandle->begin();
        iter != pixelRecHitsHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripMatchedRecHit2D> > stripMatchedRecHitsHandle;
  iEvent.getByLabel( stripMatchedRecHitTags_, stripMatchedRecHitsHandle );
  std::cout << "====> number of strip matched reco hits " << stripMatchedRecHitsHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripMatchedRecHit2D>::const_iterator iter = stripMatchedRecHitsHandle->begin();
        iter != stripMatchedRecHitsHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripRphiRecHitsHandle;
  iEvent.getByLabel( stripRphiRecHitTags_, stripRphiRecHitsHandle );
  std::cout << "====> number of strip rphi reco hits " << stripRphiRecHitsHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripRecHit2D>::const_iterator iter = stripRphiRecHitsHandle->begin();
        iter != stripRphiRecHitsHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripRphiRecHitsUnmatchedHandle;
  iEvent.getByLabel( stripRphiRecHitUnmatchedTags_, stripRphiRecHitsUnmatchedHandle );
  std::cout << "====> number of strip unmatched rphi reco hits " << stripRphiRecHitsUnmatchedHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripRecHit2D>::const_iterator iter = stripRphiRecHitsUnmatchedHandle->begin();
        iter != stripRphiRecHitsUnmatchedHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripStereoRecHitsHandle;
  iEvent.getByLabel( stripStereoRecHitTags_, stripStereoRecHitsHandle );
  std::cout << "====> number of strip stereo reco hits " << stripStereoRecHitsHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripRecHit2D>::const_iterator iter = stripStereoRecHitsHandle->begin();
        iter != stripStereoRecHitsHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripStereoRecHitsUnmatchedHandle;
  iEvent.getByLabel( stripStereoRecHitUnmatchedTags_, stripStereoRecHitsUnmatchedHandle );
  std::cout << "====> number of strip unmatched stereo reco hits " << stripStereoRecHitsUnmatchedHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripRecHit2D>::const_iterator iter = stripStereoRecHitsUnmatchedHandle->begin();
        iter != stripStereoRecHitsUnmatchedHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  /*
  using namespace edm;
  using namespace std;
  int nhit = 0;
  for ( trackingRecHit_iterator i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++ ) {
    cout << "hit #" << nhit;
    TransientTrackingRecHit::RecHitPointer hit=builder->build(&**i );
    DetId hitId = hit->geographicalId();
    if ( hitId.det() == DetId::Tracker ) {
      if ( hitId.subdetId() == StripSubdetector::TIB )  
        cout << " - TIB " << TIBDetId(hitId).layer();
      else if ( hitId.subdetId() == StripSubdetector::TOB ) 
        cout << " - TOB " << TOBDetId(hitId).layer();
      else if ( hitId.subdetId() == StripSubdetector::TEC ) 
        cout << " - TEC " << TECDetId(hitId).wheel();
      else if ( hitId.subdetId() == StripSubdetector::TID ) 
        cout << " - TID " << TIDDetId(hitId).wheel();
      else if ( hitId.subdetId() == StripSubdetector::TID ) 
        cout << " - TID " << TIDDetId(hitId).wheel();
      else if ( hitId.subdetId() == (int) PixelSubdetector::PixelBarrel ) 
        cout << " - PixBar " << PXBDetId(hitId).layer();
      else if ( hitId.subdetId() == (int) PixelSubdetector::PixelEndcap )
        cout << " - PixFwd " << PXFDetId(hitId).disk();
      else 
        cout << " UNKNOWN TRACKER HIT TYPE ";
    }
    if ( hit->isValid() ) cout << " - globalPos =" << hit->globalPosition() << endl;
    else cout << " - invalid hit" << endl;
    nhit++;
  }
  cout << endl;
  */
}


void
HTProdRecHits::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
  std::cout << "Builder name: " << builderName_ << std::endl;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName_,theBuilder);
  builder_ = theBuilder.product();
}


void 
HTProdRecHits::beginJob()
{
  file_ = new TFile(fileName_.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  hNumPixelClusters_ = new TH1F("pixelSize","size",1000,0,100);
  hNumStripClusters_ = new TH1F("stripSize","size",1000,0,100);
  TH1::AddDirectory(oldAddDir);
}


void 
HTProdRecHits::endJob() 
{
  file_->Write();
  file_->Close();
}


void
HTProdRecHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //desc.add<edm::InputTag>("pixelClusters",edm::InputTag("siPixelClusters"));
  //desc.add<edm::InputTag>("stripClusters",edm::InputTag("siStripClusters"));
  desc.add<edm::InputTag>("pixelRecHits" ,edm::InputTag("siPixelRecHits" ));
  desc.add<edm::InputTag>("stripMatchedRecHits"        ,edm::InputTag("siStripMatchedRecHits:matchedRecHit"        ));
  desc.add<edm::InputTag>("stripRphiRecHits"           ,edm::InputTag("siStripMatchedRecHits:rphiRecHit"           ));
  desc.add<edm::InputTag>("stripRphiRecHitsUnmatched"  ,edm::InputTag("siStripMatchedRecHits:rphiRecHitUnmatched"  ));
  desc.add<edm::InputTag>("stripStereoRecHits"         ,edm::InputTag("siStripMatchedRecHits:stereoRecHit"         ));
  desc.add<edm::InputTag>("stripStereoRecHitsUnmatched",edm::InputTag("siStripMatchedRecHits:stereoRecHitUnmatched"));
  desc.add<std::string>  ("TTRHBuilder" ,std::string("WithTrackAngle"));
  desc.addUntracked<std::string>("hfile",std::string("outfile.hbook" ));
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HTProdRecHits);
