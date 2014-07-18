// -*- C++ -*-
//
// Package:    HTGetRecHits
// Class:      HTGetRecHits
// 
/**\class HTGetRecHits HTGetRecHits.cc MLoVetere/HTGetRecHits/src/HTGetRecHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maurizio Lo Vetere,559 R-009,+41227675905,
//         Created:  Fri Nov 30 21:19:49 CET 2012
// $Id$
//
//

#include <iostream>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "TH1.h"
#include "TFile.h"


class HTGetRecHits : public edm::EDAnalyzer {
  public:
    explicit      HTGetRecHits     ( const edm::ParameterSet& );
                 ~HTGetRecHits     ( );
    static  void  fillDescriptions ( edm::ConfigurationDescriptions& descriptions );
  private:
    virtual void  analyze  ( const edm::Event&, const edm::EventSetup& );
    virtual void  beginJob ( );
    virtual void  endJob   ( );
  private:
    edm::InputTag pixelClusterTags_;
    edm::InputTag stripClusterTags_;
    edm::InputTag pixelRecHitTags_;
    edm::InputTag stripMatchedRecHitTags_;
    edm::InputTag stripRphiRecHitTags_;
    edm::InputTag stripRphiRecHitUnmatchedTags_;
    edm::InputTag stripStereoRecHitTags_;
    edm::InputTag stripStereoRecHitUnmatchedTags_;
    std::string   fileName_;
    TFile * file_;
    TH1F  * hPixelSize_;
    TH1F  * hStripSize_;
};


HTGetRecHits::HTGetRecHits(const edm::ParameterSet& iConfig)
 :  pixelClusterTags_(iConfig.getUntrackedParameter<edm::InputTag>("pixelClusters")),
    stripClusterTags_(iConfig.getUntrackedParameter<edm::InputTag>("stripClusters")),
    pixelRecHitTags_ (iConfig.getUntrackedParameter<edm::InputTag>("pixelRecHits" )),
    stripMatchedRecHitTags_        (iConfig.getUntrackedParameter<edm::InputTag>("stripMatchedRecHits"        )),
    stripRphiRecHitTags_           (iConfig.getUntrackedParameter<edm::InputTag>("stripRphiRecHits"           )),
    stripRphiRecHitUnmatchedTags_  (iConfig.getUntrackedParameter<edm::InputTag>("stripRphiRecHitsUnmatched"  )),
    stripStereoRecHitTags_         (iConfig.getUntrackedParameter<edm::InputTag>("stripStereoRecHits"         )),
    stripStereoRecHitUnmatchedTags_(iConfig.getUntrackedParameter<edm::InputTag>("stripStereoRecHitsUnmatched")),
    fileName_        (iConfig.getUntrackedParameter<std::string>  ("hfile"        ))
{
}


HTGetRecHits::~HTGetRecHits()
{
}


void
HTGetRecHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelClustersHandle;
  iEvent.getByLabel( pixelClusterTags_, pixelClustersHandle );
  std::cout << "====> number of pixel clusters " << pixelClustersHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator iter = pixelClustersHandle->begin();
        iter != pixelClustersHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > stripClustersHandle;
  iEvent.getByLabel( stripClusterTags_, stripClustersHandle );
  std::cout << "====> number of strip clusters " << stripClustersHandle->size() << std::endl;
  for ( edmNew::DetSetVector<SiStripCluster>::const_iterator iter = stripClustersHandle->begin();
        iter != stripClustersHandle->end();
        iter++ ) {
    std::cout << "      " << iter->detId() << " [" << iter->size() << "]" << std::endl;
  }
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

  //using reco::TrackCollection;
  //Handle<TrackCollection> tracks;
  //iEvent.getByLabel(trackTags_,tracks);
  //for ( TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack ) {
  //  int charge = 0;
  //  charge = itTrack->charge();  
  //  histo_->Fill(itTrack->pt());
  //}
}


void 
HTGetRecHits::beginJob()
{
  file_ = new TFile(fileName_.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  hPixelSize_ = new TH1F("pixelSize","size",1000,0,100);
  hStripSize_ = new TH1F("stripSize","size",1000,0,100);
  TH1::AddDirectory(oldAddDir);
}


void 
HTGetRecHits::endJob() 
{
  file_->Write();
  file_->Close();
}


void
HTGetRecHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("pixelClusters",edm::InputTag("siPixelClusters"));
  desc.addUntracked<edm::InputTag>("stripClusters",edm::InputTag("siStripClusters"));
  desc.addUntracked<edm::InputTag>("pixelRecHits" ,edm::InputTag("siPixelRecHits" ));
  desc.addUntracked<edm::InputTag>("stripMatchedRecHits"        ,edm::InputTag("siStripMatchedRecHits:matchedRecHit"        ));
  desc.addUntracked<edm::InputTag>("stripRphiRecHits"           ,edm::InputTag("siStripMatchedRecHits:rphiRecHit"           ));
  desc.addUntracked<edm::InputTag>("stripRphiRecHitsUnmatched"  ,edm::InputTag("siStripMatchedRecHits:rphiRecHitUnmatched"  ));
  desc.addUntracked<edm::InputTag>("stripStereoRecHits"         ,edm::InputTag("siStripMatchedRecHits:stereoRecHit"         ));
  desc.addUntracked<edm::InputTag>("stripStereoRecHitsUnmatched",edm::InputTag("siStripMatchedRecHits:stereoRecHitUnmatched"));
  desc.addUntracked<std::string>  ("hfile"        ,std::string  ("outfile.hbook"  ));
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HTGetRecHits);
