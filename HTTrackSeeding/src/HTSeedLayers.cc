// -*- C++ -*-
//
// Package:    HTSeedLayers
// Class:      HTSeedLayers
// 
/**\class HTSeedLayers HTSeedLayers.cc MLoVetere/HTSeedLayers/src/HTSeedLayers.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maurizio Lo Vetere,559 R-009,+41227675905,
//         Created:  Fri Nov 30 21:19:49 CET 2012
// $Id: HTSeedLayers.cc,v 1.6 2013/04/09 07:39:09 mlv Exp $
//
//

#include <iostream>

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"

#include "TH1.h"
#include "TFile.h"

#include "RecoTracker/TkHitPairs/interface/LayerHitMapCache.h"


class HTSeedLayers : public edm::EDAnalyzer {
  public:
    explicit      HTSeedLayers     ( const edm::ParameterSet& );
                 ~HTSeedLayers     ( );
    static  void  fillDescriptions ( edm::ConfigurationDescriptions& descriptions );
  private:
    void          produceSeeds   ( ctfseeding::SeedingLayer::Hits hits );
    virtual void  analyze  ( const edm::Event&, const edm::EventSetup& );
    virtual void  beginRun ( edm::Run const&  , edm::EventSetup const& );
    virtual void  beginJob ( );
    virtual void  endJob   ( );
  private:
    std::string  theLayerListName_;
    std::string  fileName_;
    TFile     *  file_;
    TH1F      *  numPixelClusters_;
    TH1F      *  numStripClusters_;
    ctfseeding::SeedingLayerSets theLayerSets_;
};


HTSeedLayers::HTSeedLayers(const edm::ParameterSet& iConfig)
  : theLayerListName_ (iConfig.getParameter<std::string> ("seedingLayers")),
    fileName_         (iConfig.getUntrackedParameter<std::string>("hfile"))
{
}


HTSeedLayers::~HTSeedLayers()
{
}


void
HTSeedLayers::produceSeeds(ctfseeding::SeedingLayer::Hits hits)
{
  std::cout << "  SeedingLayers with " << hits.size() << " hits in the set" << std::endl; 
  if ( hits.size()==0 ) return;
  int nhit = 0;
  for ( ctfseeding::SeedingLayer::Hits::const_iterator hit = hits.begin(); 
        hit != hits.end();
        hit++ ) {
    std::cout << "    hit # " << nhit;
    DetId hitId = (*hit)->geographicalId();
    if ( hitId.det() == DetId::Tracker ) {
      if      ( hitId.subdetId() == StripSubdetector::TIB )  
         std::cout << " - TIB " << TIBDetId(hitId).layer();
      else if ( hitId.subdetId() == StripSubdetector::TOB ) 
         std::cout << " - TOB " << TOBDetId(hitId).layer();
      else if ( hitId.subdetId() == StripSubdetector::TID ) 
         std::cout << " - TID " << TIDDetId(hitId).wheel();
      else if ( hitId.subdetId() == StripSubdetector::TEC ) 
         std::cout << " - TEC " << TECDetId(hitId).wheel();
      else if ( hitId.subdetId() == PixelSubdetector::PixelBarrel ) 
         std::cout << " - PixBar " << PXBDetId(hitId).layer();
      else if ( hitId.subdetId() == PixelSubdetector::PixelEndcap )
         std::cout << " - PixFwd " << PXFDetId(hitId).disk();
      else 
         std::cout << " UNKNOWN TRACKER HIT TYPE ";
    }
    if ( (*hit)->isValid() ) std::cout << " - globalPos =" << (*hit)->globalPosition() << std::endl;
    else std::cout << " - invalid hit" << std::endl;
    nhit++;
  }
}


void
HTSeedLayers::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  for ( ctfseeding::SeedingLayerSets::const_iterator aSet = theLayerSets_.begin(); 
        aSet != theLayerSets_.end();
        aSet++ ) {
    ctfseeding::SeedingLayer::Hits oneLayerSetHits; 
    for ( ctfseeding::SeedingLayers::const_iterator aLayer = aSet->begin(); 
          aLayer != aSet->end();
          aLayer++ ) {
      ctfseeding::SeedingLayer::Hits hits = aLayer->hits(iEvent,iSetup);
      oneLayerSetHits.insert(oneLayerSetHits.end(),hits.begin(),hits.end());
    }
    produceSeeds(oneLayerSetHits);
  }
}


void
HTSeedLayers::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
  std::cout << "SeedingLayerSetsBuilder name: " << theLayerListName_ << std::endl;
  edm::ESHandle<SeedingLayerSetsBuilder> layerBuilder;
  setup.get<TrackerDigiGeometryRecord>().get(theLayerListName_.c_str(),layerBuilder);
  theLayerSets_ = layerBuilder->layers(setup); 
  int i=0;
  for ( ctfseeding::SeedingLayerSets::const_iterator aSet = theLayerSets_.begin(); 
        aSet != theLayerSets_.end();
        aSet++ ) {
    std::cout << "SeedingLayerSet number " << ++i << std::endl;
    for ( ctfseeding::SeedingLayers::const_iterator aLayer = aSet->begin(); 
           aLayer != aSet->end();
           aLayer++ ) {
      std::cout << "  " << aLayer->name() << std::endl; 
    }
  }
}


void 
HTSeedLayers::beginJob()
{
  file_ =  new TFile(fileName_.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  numPixelClusters_ = new TH1F("pixelSize","size",1000,0,100);
  numStripClusters_ = new TH1F("stripSize","size",1000,0,100);
  TH1::AddDirectory(oldAddDir);
}


void 
HTSeedLayers::endJob() 
{
  file_->Write();
  file_->Close();
}


void
HTSeedLayers::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>         ("seedingLayers",std::string("HoughTransformSeedLayersAllHitsOneSet"));
  desc.addUntracked<std::string>("hfile"        ,std::string("outfile.hbook"                        ));
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HTSeedLayers);
