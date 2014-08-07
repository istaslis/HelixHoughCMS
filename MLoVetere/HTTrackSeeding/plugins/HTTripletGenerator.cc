#include "MLoVetere/HTTrackSeeding/plugins/HTTripletGenerator.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include <string>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"

typedef TransientTrackingRecHit::ConstRecHitPointer SeedingHit;

#include <map>

using namespace ctfseeding;


HTTripletGenerator::HTTripletGenerator( const edm::ParameterSet & cfg )
{
  edm::LogInfo("HTTrackSeeding|HTTripletGenerator") << "Constructing HTTripletGenerator";
  edm::ParameterSet nBinsPSet  = cfg.getParameter<edm::ParameterSet>("DivPSet");
  theNumBins.nCurv( nBinsPSet.getParameter<unsigned int>("nBinsCurv") );
  theNumBins.nEta ( nBinsPSet.getParameter<unsigned int>("nBinsEta" ) );
  theNumBins.nPhi ( nBinsPSet.getParameter<unsigned int>("nBinsPhi" ) );
  theNumBins.nTip ( nBinsPSet.getParameter<unsigned int>("nBinsTip" ) );
  theNumBins.nLip ( nBinsPSet.getParameter<unsigned int>("nBinsLip" ) );
  LogTrace("HTTrackSeeding") 
    << "theDivCurv="  << theNumBins.nCurv()
    << "; theDivEta=" << theNumBins.nEta ()
    << "; theDivPhi=" << theNumBins.nPhi ()
    << "; theDivTip=" << theNumBins.nTip ()
    << "; theDivLip=" << theNumBins.nLip ();
  edm::ParameterSet minResPSet = cfg.getParameter<edm::ParameterSet>("MinResPSet");
  theMinRes.dCurv( minResPSet.getParameter<double>("resCurv") );
  theMinRes.dEta ( minResPSet.getParameter<double>("resEta" ) );
  theMinRes.dPhi ( minResPSet.getParameter<double>("resPhi" ) );
  theMinRes.dTip ( minResPSet.getParameter<double>("resTip" ) );
  theMinRes.dLip ( minResPSet.getParameter<double>("resLip" ) );
  LogTrace("HTTrackSeeding") 
    << "theMinResCurv=" << theMinRes.dCurv() << "cm-1;" 
    << " theMinResEta=" << theMinRes.dEta () << ";"
    << " theMinResPhi=" << theMinRes.dPhi () << "rad;"
    << " theMinResTip=" << theMinRes.dTip () << "cm;"
    << " theMinResLip=" << theMinRes.dLip () << "cm";
  edm::ParameterSet maxResPSet = cfg.getParameter<edm::ParameterSet>("MaxResPSet");
  theMaxRes.dCurv( maxResPSet.getParameter<double>("resCurv") );
  theMaxRes.dEta ( maxResPSet.getParameter<double>("resEta" ) );
  theMaxRes.dPhi ( maxResPSet.getParameter<double>("resPhi" ) );
  theMaxRes.dTip ( maxResPSet.getParameter<double>("resTip" ) );
  theMaxRes.dLip ( maxResPSet.getParameter<double>("resLip" ) );
  LogTrace("HTTrackSeeding") 
    << "theMaxResCurv=" << theMaxRes.dCurv() << "cm-1;"  
    << " theMaxResEta=" << theMaxRes.dEta () << ";"
    << " theMaxResPhi=" << theMaxRes.dPhi () << "rad;"
    << " theMaxResTip=" << theMaxRes.dTip () << "cm;"
    << " theMaxResLip=" << theMaxRes.dLip () << "cm";
  edm::ParameterSet turnsPSet = cfg.getParameter<edm::ParameterSet>("HalfTurns");
  thePHTurns = turnsPSet.getParameter<unsigned int>("positive");
  theNHTurns = turnsPSet.getParameter<unsigned int>("negative");
  LogTrace("HTTrackSeeding")
    << "thePHTurns=" << thePHTurns
    << " theNHTurns=" << theNHTurns;
  theRequiredLayers = cfg.getParameter<unsigned int>("RequiredLayers");
  LogTrace("HTTrackSeeding")
    << "theRequiredLayers=" << theRequiredLayers;
  theLayerBuilderName = cfg.getParameter<std::string>("SeedingLayers");
  vertexSrc = cfg.getParameter<std::string>("VertexSrc");
   seedSrc = cfg.getParameter<std::string>("SeedSrc");
   seedfromHits = cfg.getParameter<bool>("SeedsFromHits");
}


void  HTTripletGenerator::init( const edm::EventSetup & es, const TrackingRegion & reg  )
{
  edm::ESHandle<SeedingLayerSetsBuilder > layerBuilder;
  es.get<TrackerDigiGeometryRecord>().get(theLayerBuilderName, layerBuilder);
  theLayerSets = layerBuilder->layers(es);
  edm::ESHandle<MagneticField > fieldESH;
  es.get<IdealMagneticFieldRecord>().get(fieldESH);
  theField = fieldESH->inTesla(GlobalPoint(0,0,0)).z() * 2.99792458e-3F;  // GeV/cm units
  LogDebug("HTTrackSeeding") << "SeedingLayerSets";
  int i=0;
  for ( ctfseeding::SeedingLayerSets::const_iterator aSet = theLayerSets.begin(); aSet != theLayerSets.end(); aSet++ ) {
    LogTrace("HTTrackSeeding") << "SeedingLayerSet number " << ++i;
    for ( ctfseeding::SeedingLayers::const_iterator aLayer = aSet->begin(); aLayer != aSet->end(); aLayer++ ) {
      LogTrace("HTTrackSeeding") << "  " << aLayer->name();
    std::cout <<" Layer name "<<aLayer->name()<<std::endl;
    }
  }
  theRefPoint      = reg.origin();
  float aCurvBound = theField/reg.ptMin();  // curv in [ -aCurvBound, aCurvBound ]
  float aZBound    = reg.originZBound();    // lip  in [ -aZBound   , aZBound    ]
  float aRBound    = reg.originRBound();    // tip  in [ -aRBound   , aRBound    ]
  TrackingRegion::Range aEtaRange;
  TrackingRegion::Range aPhiRange;
  if ( const RectangularEtaPhiTrackingRegion * etaPhiReg = dynamic_cast<const RectangularEtaPhiTrackingRegion * >(&reg) ) {
    aEtaRange = etaPhiReg->etaRange();
    RectangularEtaPhiTrackingRegion::Margin phiMargin = etaPhiReg->phiMargin();
    aPhiRange = TrackingRegion::Range(reg.direction().phi()-phiMargin.left(),reg.direction().phi()+phiMargin.right());
  } else {
    aEtaRange = TrackingRegion::Range( -2.5, 2.5);
    aPhiRange = TrackingRegion::Range(-M_PI,M_PI);
  }
  //LogTrace("HTTrackSeeding") 
  std::cout
  << "X=" << theRefPoint.x() << "cm; Y=" << theRefPoint.y() << "cm; Z=" << theRefPoint.z() << "cm; " << std::endl
    << "theCurvBound=" << aCurvBound << "cm^-1; " 
    << "theZBound="    << aZBound << "cm; " 
    << "theRBound="    << aRBound << "cm; "
    << "theEtaRange="  << aEtaRange << "; thePhiRange=" << aPhiRange << "rad";
  theRange = HelixParRange( -aCurvBound     , aCurvBound     , 
                             aEtaRange.min(), aEtaRange.max(),
                            -aZBound        , aZBound        ,
	  		     aPhiRange.min(), aPhiRange.max(),
                            -aRBound        , aRBound        ,
			    thePHTurns      , theNHTurns     ); 
  
  
  
  
  
 
  
  
}


#include "MLoVetere/HTTrackSeeding/interface/HelixHoughInterface.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"

PHENIXHough::SimpleHit3D MakeHit(GlobalPoint &pos, GlobalError& err, int index, int layer)
{
    return PHENIXHough::SimpleHit3D(pos.x(),sqrt(err.cxx())/2*3.14,
                                    pos.y(),sqrt(err.cyy())/2*3.14,
                                    pos.z(),sqrt(err.czz())/2*3.14, index, layer);

}

bool comparehits(PHENIXHough::SimpleHit3D h, PHENIXHough::SimpleHit3D g)
{
    return (h.x*h.x+h.y*h.y < g.x*g.x+g.y*g.y);
}


//void GetSeedHits( const edm::Event & ev, const edm::EventSetup & iSetup, std::string seedSrc, std::vector<PHENIXHough::SimpleTrack3D> &seedTracks, std::vector<TransientTrackingRecHitBuilder::RecHitPointer> &seedhits)
//{
//    std::string tripletSrc=seedSrc;//"hiPixel3PrimTracks";

//    edm::Handle<OrderedHitTriplets> triplets;
//    ev.getByLabel(tripletSrc,triplets);


//    int index = 0;
//    for(unsigned it=0; it<triplets->size(); ++it){
//        const OrderedHitTriplet & triplet = (*triplets)[it];
//        PHENIXHough::SimpleTrack3D seed;

//        edm::ESHandle<TransientTrackingRecHitBuilder> recHitBuilder;
//        iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",recHitBuilder);


//        std::vector<SeedingHitSet::ConstRecHitPointer> hitsintriplet;
//        hitsintriplet.push_back(triplet.inner());
//        hitsintriplet.push_back(triplet.middle());
//        hitsintriplet.push_back(triplet.outer());

//        for (int i=0;i<3;i++) {
//            const TrackingRecHit & hit = *(hitsintriplet[i]);
//            DetId detID = hit.geographicalId()();
//            if(!detID) continue;

//            TransientTrackingRecHitBuilder::RecHitPointer rechit = recHitBuilder->build(&hit);

//            GlobalPoint globalPos = rechit->globalPosition();
//            GlobalError globalErr = rechit->globalPositionError();

//            PHENIXHough::SimpleHit3D houghhit = MakeHit(globalPos, globalErr, index, i); // so far layer is the index in triplet...
//            seed.hits.push_back(houghhit);

//            seedhits.push_back(rechit);

//            index++;



//           // std::cout << "   hit "<<(unsigned)detID<<" x="<<globalPos.x()<<" y="<<globalPos.y()<<" z="<<globalPos.z()<<std::endl;
//        }


//            std::sort(seed.hits.begin(),seed.hits.end(),comparehits);
//            for (unsigned i=0;i<seed.hits.size();i++) seed.hits[i].layer = i; //change not only the position in array, but the logical layer!


//        seedTracks.push_back(seed);

//    }
//}

void GetSeedHits( const edm::Event & ev, const edm::EventSetup & iSetup, std::string seedSrc, std::vector<PHENIXHough::SimpleTrack3D> &seedTracks, std::vector<TransientTrackingRecHitBuilder::RecHitPointer> &seedhits)
{
    std::string tripletSrc=seedSrc;

    edm::Handle<TrajectorySeedCollection> triplets;
    ev.getByLabel(tripletSrc,triplets);

    TrajectorySeedCollection *seeds = const_cast<TrajectorySeedCollection *> (triplets.product());

    int index = 0;
    for(unsigned it=0; it<seeds->size(); ++it){
        const TrajectorySeed & theTrackerSeed = (*seeds)[it];
        PHENIXHough::SimpleTrack3D seed;

        edm::ESHandle<TransientTrackingRecHitBuilder> recHitBuilder;
        iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",recHitBuilder);


        TrajectorySeed::range theSeedRange=theTrackerSeed.recHits();
        TrajectorySeed::const_iterator theSeedRangeIteratorBegin = theSeedRange.first;
        TrajectorySeed::const_iterator theSeedRangeIteratorEnd   = theSeedRange.second;
        TrajectorySeed::const_iterator theSeedItr = theSeedRangeIteratorBegin;


        for ( ; theSeedItr != theSeedRangeIteratorEnd; ++theSeedItr ) {
            const TrackingRecHit & hit = *(theSeedItr);
            DetId detID = hit.geographicalId()();
            if(!detID) continue;

            TransientTrackingRecHitBuilder::RecHitPointer rechit = recHitBuilder->build(&hit);

            GlobalPoint globalPos = rechit->globalPosition();
            GlobalError globalErr = rechit->globalPositionError();

            PHENIXHough::SimpleHit3D houghhit = MakeHit(globalPos, globalErr, index, theSeedItr-theSeedRangeIteratorBegin); // so far layer is the index in triplet...
            seed.hits.push_back(houghhit);

            seedhits.push_back(rechit);

            index++;

        }


            std::sort(seed.hits.begin(),seed.hits.end(),comparehits);
            for (unsigned i=0;i<seed.hits.size();i++) seed.hits[i].layer = i; //change not only the position in array, but the logical layer!


        seedTracks.push_back(seed);

    }
}



void GetSeedHitsFromTracks( const edm::Event & ev, const edm::EventSetup & iSetup, std::string seedSrc, std::vector<PHENIXHough::SimpleTrack3D> &seedTracks, std::vector<TransientTrackingRecHitBuilder::RecHitPointer> &seedhits)
{
    std::string trkSrc=seedSrc;

    edm::Handle<std::vector<reco::Track> > tracks;
    ev.getByLabel(trkSrc,tracks);


    int index = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        const reco::Track & etrk = (*tracks)[it];
        PHENIXHough::SimpleTrack3D seed;

        edm::ESHandle<TransientTrackingRecHitBuilder> recHitBuilder;
        iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",recHitBuilder);

        for (trackingRecHit_iterator ith = etrk.recHitsBegin(); ith != etrk.recHitsEnd(); ++ith) {
            const TrackingRecHit * hit = ith->get();
            DetId detID = hit->geographicalId()();

            if(!detID) continue;

            TransientTrackingRecHitBuilder::RecHitPointer rechit = recHitBuilder->build(hit);

            GlobalPoint globalPos = rechit->globalPosition();
            GlobalError globalErr = rechit->globalPositionError();

            PHENIXHough::SimpleHit3D houghhit = MakeHit(globalPos, globalErr, index, ith-etrk.recHitsBegin()); // here layer is the index in triplet - enough for now
            seed.hits.push_back(houghhit);
            seedhits.push_back(rechit);
            index++;
        }


            std::sort(seed.hits.begin(),seed.hits.end(),comparehits);
            for (unsigned i=0;i<seed.hits.size();i++) seed.hits[i].layer = i; //change not only the position in array, but the logical layer!


        seedTracks.push_back(seed);

    }
}



void  HTTripletGenerator::hitTriplets ( const TrackingRegion & reg, OrderedHitTriplets & prs, const edm::Event & ev, const edm::EventSetup & es )
{
  
    std::vector<PHENIXHough::SimpleTrack3D> seedTracks;
    
    float xVtx = 0, xVtxErr = 0, yVtx = 0, yVtxErr =0, zVtx = 0, zVtxErr = 0;
	int nVtx = 0;


	edm::Handle<std::vector<reco::Vertex> > vertices;
	const std::vector<reco::Vertex> * recoVertices;
 	ev.getByLabel(vertexSrc,vertices);
	recoVertices = vertices.product();

    nVtx = recoVertices->size();;

    //going through the vertices and taking the last one
    for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
         xVtx =(*recoVertices)[i].position().x();
         yVtx =(*recoVertices)[i].position().y();
         zVtx =(*recoVertices)[i].position().z();
         xVtxErr = (*recoVertices)[i].xError();
         yVtxErr = (*recoVertices)[i].yError();
         zVtxErr = (*recoVertices)[i].zError();
     }


  init(es, reg);
  int i=0;
  for ( SeedingLayerSets::const_iterator aLayerSet=theLayerSets.begin(); aLayerSet!=theLayerSets.end(); aLayerSet++, i++ ) {
    //    if ( aLayerSet->size()<3 ) throw cms::Exception("HTTripletGenerator") << "You are trying to use a set with " << aLayerSet->size() << " layers instead of at least 3 ";	
    
    TrackingRegion::Hits hits;
    std::vector<unsigned int> hitsInLayers;
    int k=0;
    for ( SeedingLayers::const_iterator layer=aLayerSet->begin(); layer!=aLayerSet->end(); layer++ ) {
      TrackingRegion::Hits ht = reg.hits(ev,es, &(*layer));
      hits.insert(hits.end(),ht.begin(),ht.end());
    
      hitsInLayers.push_back(ht.size());
      std::cout <<"Layer num :"<<k<<"; Layer name: "<<layer->name()<<std::endl;
      k++;
    }
    if ( hits.size()<3 ) continue;
      
    std::vector<TransientTrackingRecHitBuilder::RecHitPointer> seedhits;
    if (seedfromHits)
        GetSeedHits(ev,es, seedSrc, seedTracks, seedhits);
    else
        GetSeedHitsFromTracks(ev,es, seedSrc, seedTracks, seedhits);


    HelixHough finder( theRefPoint, theRange, theNumBins, theMinRes, theMaxRes, theRequiredLayers );
    int min_hits =  3;
    int max_hits = 10;
    finder.nEv = (unsigned int)ev.id().event();
      
    finder.SetVertex(xVtx,xVtxErr,yVtx,yVtxErr,zVtx,zVtxErr);
    finder.SetSeeds(seedTracks, seedhits);
    finder.findHelices( hits, hitsInLayers, min_hits, max_hits, prs );
  }
  std::cout << "Triplet generator ending with " << prs.size() << " triplets" << std::endl;
}
