// -*- C++ -*-
//
// Package:    HTGetClusters
// Class:      HTGetClusters
// 
/**\class HTGetClusters HTGetClusters.cc MLoVetere/HTGetClusters/src/HTGetClusters.cc

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

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "TH1.h"
#include "TFile.h"


class HTGetClusters : public edm::EDAnalyzer {
  public:
    explicit      HTGetClusters    ( const edm::ParameterSet& );
                 ~HTGetClusters    ( );
    static  void  fillDescriptions ( edm::ConfigurationDescriptions& descriptions );
  private:
    virtual void  analyze  ( const edm::Event&, const edm::EventSetup& );
    virtual void  beginJob ( );
    virtual void  endJob   ( );
  private:
    edm::InputTag pixelTags_;
    edm::InputTag stripTags_;
    std::string   hFileName_;
    TFile * file_;
    TH1F  * hPixelSize_;
    TH1F  * hStripSize_;
};


HTGetClusters::HTGetClusters(const edm::ParameterSet& iConfig)
  : pixelTags_(iConfig.getUntrackedParameter<edm::InputTag>("pixels")),
    stripTags_(iConfig.getUntrackedParameter<edm::InputTag>("strips")),
    hFileName_(iConfig.getUntrackedParameter<std::string>("hfile"))
{
}


HTGetClusters::~HTGetClusters()
{
}


void
HTGetClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixels;
  iEvent.getByLabel( pixelTags_, pixels );
  hPixelSize_->Fill( pixels.product()->size() );
  //std::cout << "====> number of pixel clusters " << pixels->size() << std::endl;
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > strips;
  iEvent.getByLabel( stripTags_, strips );
  hStripSize_->Fill( strips.product()->size() );
  //std::cout << "====> number of strip clusters " << strips->size() << std::endl;
}


void 
HTGetClusters::beginJob()
{
  file_ = new TFile(hFileName_.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  hPixelSize_ = new TH1F("pixelSize","size",1000,0,100);
  hStripSize_ = new TH1F("stripSize","size",1000,0,100);
  TH1::AddDirectory(oldAddDir);
}


void 
HTGetClusters::endJob() 
{
  file_->Write();
  file_->Close();
}


void
HTGetClusters::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 edm::ParameterSetDescription desc;
 desc.addUntracked<edm::InputTag>("pixels",edm::InputTag("siPixelClusters"));
 desc.addUntracked<edm::InputTag>("strips",edm::InputTag("siStripClusters"));
 desc.addUntracked<std::string>  ("hfile" ,std::string  ("outfile.hbook"  ));
 descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HTGetClusters);
