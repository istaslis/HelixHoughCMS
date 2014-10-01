#include "MLoVetere/HTTrackSeeding/interface/HelixHoughInterface.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixHoughEngine.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "MLoVetere/HTTrackSeeding/interface/sPHENIXTracker.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3DInterface.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3DInterface.h"

#include "TFile.h"
#include "TTree.h"


const float Bfield = 3.8;
// pT in GeV/c , B in T , kappa in cm^-1
static inline double pT_to_kappa(double pT) {
	return 0.003 * Bfield / pT;
}

static inline double kappa_to_pT(double kappa) {
	return 0.003 * Bfield / kappa;
}


static int GetLayerNumber(PHENIXHough::SimpleHit3D &hit) {
  //  float layerBounds[] = { 3., 6., 9., /*20., 30., 38., 45.,*/ 55., 65., 73., 82.,
  //		92., 100, 112. };

  float layerBounds[] = { 3., 6., 9., 20., 30., 38., 45., 55.,  65., 73., 82.,
			  92., 100, 112. };


	double radius = sqrt(hit.x * hit.x + hit.y * hit.y);
	//cout << radius<<endl;
	unsigned int layer;
	for (layer = 0; layer < sizeof(layerBounds) / sizeof(float); layer++)
		if (layerBounds[layer] > radius)
			break;

	return layer - 1;
}

void SetProfile(std::vector<unsigned int> &zoomprofile, unsigned int phi, unsigned int d, unsigned int kappa, unsigned int dzdl, unsigned int z0)
{
	zoomprofile[0] = phi;
	zoomprofile[1] = d;
	zoomprofile[2] = kappa;
	zoomprofile[3] = dzdl;
	zoomprofile[4] = z0;
}

void ConvertTracks(std::vector<PHENIXHough::SimpleTrack3D> &input, std::vector< ::SimpleHit3D> &hitsList, std::vector< ::SimpleTrack3D> &output)
{
  for (unsigned int i=0;i<input.size(); i++)
  {
    PHENIXHough::SimpleTrack3D &track = input[i];
    ::SimpleTrack3D outtrack;
    for (unsigned int j=0;j<track.hits.size();j++)
    {
        outtrack.hits.push_back(hitsList[track.hits[j].index]);
    }
      
    output.push_back(outtrack);
    
  }
  
}

void PrintMCTrack(PHENIXHough::SimpleTrack3D track)
{
    std::cout << " Track hits: "<<track.hits.size()<<" phi="<<track.phi<<" d="<<track.d<<" kappa="<<track.kappa<<" dzdl="<<track.dzdl<<" z0="<<track.z0<<" index="<<track.index<<std::endl;
}

void DoTheJob(std::vector< ::SimpleHit3D> &hitsList, 
	      std::vector<unsigned int> &HitsInLayers, 
	      std::vector<PHENIXHough::SimpleTrack3D> &seeds, 
	      std::vector< ::SimpleTrack3D> &tracks, 
	      float xVtx, float yVtx, float zVtx, float zVtxErr, 
	      unsigned int nEv, 
	      std::vector<std::vector<unsigned int> > &zoomprofile_2,
	      double phimin,
	      double phimax,
	      double SV,
	      double dzdlmin,
	      double dzdlmax,
	      double pTmin,
	      double debugOutput,
	      double HoughTransformSeeding,
	      double SeedRefitterOnly)
{
  std::vector<PHENIXHough::SimpleHit3D> hits;
  std::vector<PHENIXHough::SimpleHit3D> hits_seeded;

  bool OuterFiltering = !SeedRefitterOnly;
  bool houghseeding = HoughTransformSeeding;
  bool debug = debugOutput;  
  
  //float ro1, ro2, phi1, phi2, z1, z2;
  
  for (unsigned int i=1;i<HitsInLayers.size();i++)
    HitsInLayers[i]+=HitsInLayers[i-1];
  
  int  currentlayer = 0;
  for (unsigned int i=0;i<hitsList.size();i++)
  {
    
    float dx = fabs(hitsList[i].x().upper()-hitsList[i].x().lower())*3.4;
    float dy = fabs(hitsList[i].y().upper()-hitsList[i].y().lower())*3.4;
    float dz = fabs(hitsList[i].z().upper()-hitsList[i].z().lower())*3.4;

    if (dz<0.01) dz=0.01;
   
    float x = hitsList[i].x().center();
    float y = hitsList[i].y().center();
    float z = hitsList[i].z().center();


    PHENIXHough::SimpleHit3D hit(x-xVtx, dx, y-yVtx, dy, z, dz, i+40000, 0);
    
    //get layer info from original seeding sequence
    if (i==HitsInLayers[ currentlayer])
      currentlayer++;


    hit.layer = currentlayer;

    if (hit.layer<3)
      hits.push_back(hit);
    else 
      hits_seeded.push_back(hit);

  }
  
  if (debug) std::cout <<"hits size "<<hits.size()<<" hits_seeded size"<<hits_seeded.size()<<std::endl;

  
  //init tracker
  
  	int nlayers = 25;
	std::vector<float> radii { 4.1, 7.5, 9.92,/* 25.0, 34.0, 43.0, 52.0, */61.0, 69.6, 78.2, 86.8, 96.5, 108.0, 120,130,140,150,160,170,180,190,200,210,220,230 };

	std::vector<float> material;
    material.assign(nlayers, 0.01);

    float kappa_max = pT_to_kappa(pTmin);
	float rho_max = pow(kappa_max, 1.);
	if (debug) std::cout << "kappamax : " << kappa_max << std::endl;

	float zxmin =  zVtx-zVtxErr*3-SV;
	float zxmax = zVtx+zVtxErr*3+SV;


	//dout<< "Vertex at z :"<<matcher.eventData[ev].zVtx<<endl;

	//phi,d,kappa,dzdl,z0
	
	unsigned int nthreads = 1;
	unsigned int seedNumber = 3;
	
	if (debug) std::cout << "zxmin :"<<zxmin<<"; zxmax :"<<zxmax<<std::endl;
	
	
	//HelixRange top_range(0., 2. * M_PI, -4, 4, 0.0, rho_max, -1, 1, -20, 20);
	//HelixRange top_range(0., 2. * M_PI, -sv, sv, 0.0, rho_max, -1, 1, zxmin,zxmax);
	HelixRange top_range(phimin, phimax, -SV, SV, 0.0, rho_max, dzdlmin, dzdlmax, zxmin,zxmax);

	std::vector<unsigned int> onezoom(5, 0);
	std::vector<std::vector<unsigned int> > zoomprofile;

	int levels = 5;
	zoomprofile.assign(levels, onezoom);

	SetProfile(zoomprofile[0], 8, 1, 1, 1, 1);
	SetProfile(zoomprofile[1], 5, 5, 5, 5, 6);
	SetProfile(zoomprofile[2], 6, 4, 2, 6, 2);
	SetProfile(zoomprofile[3], 3, 3, 3, 3, 2);
	SetProfile(zoomprofile[4], 3, 3, 3, 3, 2);

	PHENIXHough::sPHENIXTracker tracker(zoomprofile, 2, top_range, material, radii, Bfield, true, nthreads);
//   sPHENIXTracker tracker(zoomprofile, 2, top_range, material, radii, Bfield);
	tracker.setNLayers(seedNumber);
	unsigned int max_hits_seed = 6;
	unsigned int min_hits_seed = seedNumber;
	tracker.setClusterStartBin(2);
	tracker.setRejectGhosts(false);
	tracker.setChi2Cut(2000);//1000000);//20
	tracker.setChi2RemovalCut(0.5);
	tracker.setPrintTimings(true);
	tracker.setVerbosity(1);
	tracker.setCutOnDca(false);
	tracker.requireLayers(seedNumber);
	tracker.setSmoothBack(false);
    tracker.setBinScale(1);
    tracker.setZBinScale(1);
	tracker.setRemoveHits(false);
	tracker.setSeparateByHelicity(true);
	tracker.setMaxHitsPairs(0);

	tracker.setkappaCut(kappa_max);//1000000);//kappa_max);

	tracker.hitsinzoomlevel.assign(levels, std::vector<int>());
	tracker.clustersinzoomlevel_phi.assign(levels, 0);
	tracker.clustersinzoomlevel_d.assign(levels, 0);
	tracker.clustersinzoomlevel_k.assign(levels, 0);
	tracker.clustersinzoomlevel_dzdl.assign(levels, 0);
	tracker.clustersinzoomlevel_z0.assign(levels, 0);


	// *********************************************************
    // setup for tracking based on seeds
	//    std::vector<std::vector<unsigned int> > zoomprofile_2;
    
    //         zoomprofile_2.assign(5, onezoom);
    // SetProfile(zoomprofile_2[0], 8, 1, 1, 1, 1);
    //        SetProfile(zoomprofile_2[1], 8, 1, 5, 5, 1);
    //        SetProfile(zoomprofile_2[2], 8, 1, 6, 6, 1);
    //        SetProfile(zoomprofile_2[3], 8, 1, 5, 5, 1);
    //        SetProfile(zoomprofile_2[4], 4, 1, 2, 6, 1);                
    /*       zoomprofile_2.assign(13, onezoom);
       SetProfile(zoomprofile_2[0], 8, 1, 1, 1, 1);
       SetProfile(zoomprofile_2[1], 5, 1, 5, 5, 1);
	SetProfile(zoomprofile_2[2], 1, 1, 1, 3, 1);
    	SetProfile(zoomprofile_2[3], 1, 1, 1, 2, 1);
    	SetProfile(zoomprofile_2[4], 3, 1, 5, 5, 1);
    	        SetProfile(zoomprofile_2[5], 1, 1, 1, 3, 1);
    	SetProfile(zoomprofile_2[6], 1, 1, 1, 2, 1);
    	SetProfile(zoomprofile_2[7], 3, 1, 2, 5, 1);
    	SetProfile(zoomprofile_2[8], 1, 1, 1, 3, 1);
    	SetProfile(zoomprofile_2[9], 1, 1, 1, 2, 1);
    	SetProfile(zoomprofile_2[10], 3, 1, 3, 5, 1);
    	SetProfile(zoomprofile_2[11], 3, 1, 3, 5, 1);
    	SetProfile(zoomprofile_2[12], 3, 1, 3, 5, 1);
    */

	//HelixRange top_range_2(0., 2. * M_PI, -sv, sv, 0.0, rho_max, -1, 1, zxmin,zxmax);
	HelixRange top_range_2(phimin, phimax, -SV, SV, 0.0, rho_max, dzdlmin, dzdlmax, zxmin,zxmax);
	sPHENIXTracker tracker_seeded(zoomprofile_2, 1, top_range_2, material, radii, Bfield);
	unsigned int max_hits = 50;
	unsigned int min_hits = 2;//4
    tracker_seeded.setClusterStartBin(2);
	tracker_seeded.setNLayers(30);
	tracker_seeded.setSeedLayer(seedNumber);
	tracker_seeded.setRejectGhosts(false); //loop from 1st hit in de-ghosting!!!
	tracker_seeded.setChi2Cut(2000.0);//20
	tracker_seeded.setPrintTimings(true);
	tracker_seeded.setVerbosity(1);
	tracker_seeded.setCutOnDca(false);
	tracker_seeded.setSmoothBack(false);
	tracker_seeded.setBinScale(0.7);
	tracker_seeded.setZBinScale(0.7);
	tracker_seeded.setSeparateByHelicity(false);
	tracker_seeded.setRemoveHits(false);
	tracker_seeded.setChi2RemovalCut(20.0);
	//		tracker_seeded.requireLayers(1);
    
    
    
	
		timeval t1, t2;
	double time1 = 0.;
	double time2 = 0.;
    //RECO!!!
    std::vector<PHENIXHough::SimpleTrack3D> tracks_seeds, tracks_seeded, filteredtracks;;
	
	tracks.clear();
	tracker.clear();

	gettimeofday(&t1, NULL);

    if (houghseeding)
        tracker.findHelices(hits, min_hits_seed, max_hits_seed, tracks_seeds);
    else //to move into tracker class!
    {
      if (debug) std::cout << "Original seeds :"<<seeds.size()<<std::endl;
        for (unsigned  i=0;i<seeds.size();i++)
        {
	     //shift on vertex
	     for (unsigned j=0;j<seeds[i].hits.size();j++)
            { seeds[i].hits[j].x-=xVtx;seeds[i].hits[j].y-=yVtx;
             //TEMPORARY SOLUTION FOR ENDCAPS?
             //hits[j].layer = GetLayerNumber(hits[j]);
             //std::cout << "?x:"<<seeds[i].hits[j].x<<"?y:"<<seeds[i].hits[j].y<<"?z:"<< seeds[i].hits[j].z<<std::endl;
         }

            tracker.findTracksBySegments(seeds[i].hits, tracks_seeds, top_range);//tracks_seeds or ADD to tracks_seeds?
        }
        if (debug) std::cout << "Converted seeds :"<<tracks_seeds.size()<<std::endl;


    }


    for (unsigned i=0;i<tracks_seeds.size();i++)
    {
        tracks_seeds[i].index = i;
    }

    gettimeofday(&t2, NULL);
	
    time1 = ((double) (t1.tv_sec) + (double) (t1.tv_usec) / 1000000.);
    time2 = ((double) (t2.tv_sec) + (double) (t2.tv_usec) / 1000000.);
    if (debug) std::cout << "nEv = "<<nEv<< " seed tracking time = " << (time2 - time1) << std::endl << std::endl;
    
    if (debug) std::cout << "Tracks! :"<<tracks_seeds.size()<<std::endl<<std::endl;

    if (hits_seeded.size()==0) OuterFiltering = false;
        
    if (OuterFiltering) {
        tracker_seeded.setSeedStates(tracker.getKalmanStates());

        gettimeofday(&t1, NULL);

	tracker_seeded.findSeededHelices(tracks_seeds, hits_seeded, min_hits, max_hits, tracks_seeded);

        gettimeofday(&t2, NULL);
        time1 = ((double) (t1.tv_sec) + (double) (t1.tv_usec) / 1000000.);
        time2 = ((double) (t2.tv_sec) + (double) (t2.tv_usec) / 1000000.);
        std::cout << "nEv = "<<nEv<< " tracking time stage 1 = " << (time2 - time1) << std::endl << std::endl;
        //     vector<float>& isolation = tracker_seeded.getIsolation();

        std::cout << "nEv = "<<nEv<<" ; "<< tracks_seeded.size() << " stage 1 tracks found" << std::endl << std::endl;


    }


    filteredtracks.clear();
    std::vector<float> filteredtracks_chi2; //chi2 of the seed refitting
    std::vector<float> filteredtracks_matched_chi2; //chi2 of the final track (with matching)
    std::vector<float> filteredtracks_nhit; //chi2 of the final track (with matching)

    if (debug) std::cout << "Lengths : tracks_seeds.size "<<tracks_seeds.size()<<" (tracker.getKalmanStates()).size "<<(tracker.getKalmanStates()).size()<<std::endl;

    for (unsigned int i=0;i<tracks_seeds.size();i++)
      if (OuterFiltering)
	{
	  if (tracker_seeded.seedWasUsed(i)) {
	filteredtracks.push_back(tracks_seeds[i]);
	filteredtracks_chi2.push_back( (tracker.getKalmanStates())[i].chi2);
   	filteredtracks_matched_chi2.push_back(tracker_seeded.getSeedMatchedChi2(i));
	filteredtracks_nhit.push_back(tracker_seeded.getSeedTrackNhit(i));
	  }   
	}
      else
	{
	  filteredtracks.push_back(tracks_seeds[i]);
	  filteredtracks_chi2.push_back( (tracker.getKalmanStates())[i].chi2);
	  filteredtracks_matched_chi2.push_back(-1);
	  filteredtracks_nhit.push_back(-1);
	} 


    if (debug) {

    float kappa, dzdl, d, phi, z0, seedchi2, matchedchi2, nhit;
    float x,y,z,dx,dy,dz, layer;


    TFile *file = new TFile("seedoutput.root","recreate");
    TTree* seedtree = new TTree("seedtree", "a tree of seeds");
    seedtree->Branch("kappa", &kappa,"kappa/F");
    seedtree->Branch("dzdl", &dzdl,"dzdl/F");
    seedtree->Branch("d", &d, "d/F");
    seedtree->Branch("phi", &phi, "phi/F");
    seedtree->Branch("z0", &z0, "z0/F");
    seedtree->Branch("seedchi2",&seedchi2, "seedchi2/F");
    seedtree->Branch("matchedchi2",&matchedchi2, "matchedchi2/F");
    seedtree->Branch("nhit",&nhit, "nhit/F");

    TTree* hittree = new TTree("hitree", "a tree of hits");
    hittree->Branch("x", &x,"x/F");
    hittree->Branch("y", &y,"y/F");
    hittree->Branch("z", &z,"z/F");
    hittree->Branch("dx",&dx,"dx/F");
    hittree->Branch("dy",&dy,"dy/F");
    hittree->Branch("dz",&dz,"dz/F");
    hittree->Branch("layer",&layer,"layer/F");
    

    for (unsigned int i=0;i<filteredtracks.size();i++)
    {
        kappa = filteredtracks[i].kappa;
        dzdl = filteredtracks[i].dzdl;
        d =  filteredtracks[i].d;
        phi= filteredtracks[i].phi;
        z0 =  filteredtracks[i].z0;
	seedchi2 = filteredtracks_chi2[i];
	matchedchi2 = filteredtracks_matched_chi2[i];
        nhit = filteredtracks_nhit[i];
        seedtree->Fill();
        
    }
        for (unsigned int i=0;i<hits_seeded.size();i++)
    {

         dx = hits_seeded[i].dx;
         dy = hits_seeded[i].dy;
         dz = hits_seeded[i].dz;

         x = hits_seeded[i].x;
         y = hits_seeded[i].y;
         z = hits_seeded[i].z;

	 layer = (float)hits_seeded[i].layer;

        hittree->Fill();


    }

    
    seedtree->Write();
    hittree->Write();
    file->Close();
    
    }
    
    
    hits.clear();
    hits_seeded.clear();

    
    if (debug) std::cout<<"Convert tracks"<<std::endl;
    
    //tracks_seeds    ->   tracks
    
    ConvertTracks(filteredtracks, hitsList, tracks);

    if (debug) std::cout<<"Converted! "<<tracks.size()<<" tracks"<<std::endl;
}









::HelixHough::HelixHough( GlobalPoint        origin         ,
                        HelixParRange      range          ,
                        HelixParNBins      nBins          ,
                        HelixParResolution minResolution  ,
                        HelixParResolution maxResolution  ,
                        unsigned int       requiredLayers ) 
  : _origin(origin), _range(range), _nBins(nBins), _minimumResolution(minResolution), _maximumResolution(maxResolution), _requiredLayers(requiredLayers),
    _decreasePerZoom(0.5), _voteTime(0), _voteTimeXY(0), _voteTimeZ(0)  
{
  _voteTime   = new SimpleTimer;
  _voteTimeXY = new SimpleTimer;
  _voteTimeZ  = new SimpleTimer;
  assert ( _voteTime && _voteTimeXY && _voteTimeZ ); 
}


::HelixHough::~HelixHough()
{ 
  delete _voteTime;
  delete _voteTimeXY;
  delete _voteTimeZ;
}

void  ::HelixHough::findHelices ( const TrackingRegion::Hits & hits      ,
				  std::vector<unsigned int> &HitsInLayers	,
				  unsigned int                 min_hits  ,
				  unsigned int                 max_hits  ,
				  OrderedHitTriplets         & tracks    ,
				  std::vector<std::vector<unsigned int> > &zoomlevels  ,
				  double phimin,
				  double phimax,
				  double SV,
				  double dzdlmin,
				  double dzdlmax,
				  double pTmin,
				  double debugOutput,
				  double HoughTransformSeeding,
				  double SeedRefitterOnly,
				  unsigned int                 maxtracks )
{

  voteTime  ().reset();
  voteTimeXY().reset();
  voteTimeZ ().reset();

  LogDebug("HTTrackSeeding") << "Hits in input to HelixHough findHelices";
  std::vector< ::SimpleHit3D>  hitsList;
  hitsList.reserve( hits.size() );
  unsigned int i = 0;
  for ( TrackingRegion::Hits::const_iterator hit =hits.begin(); hit!=hits.end(); hit++, i++ ) {


      ::SimpleHit3D ahit ((*hit),_origin,i);
     if ( !ahit.isValid() )  continue;
     LogTrace("HTTrackSeeding") << ahit;
     
     
     if (fabs(ahit.z().upper()-ahit.z().lower())>1) continue;//TODO: we don't need big errors. But do we have them???
     hitsList.push_back(ahit);
  }
  assert( hitsList.size()==hits.size() );
  if ( hitsList.size()<3 ) return;

  initEvent(hitsList, min_hits);
  
//   HelixHoughEngine engine(*this,_range,_nBins);
   std::vector< ::SimpleTrack3D> temp_tracks;
//   engine.findHelices(hitsList,min_hits,max_hits,temp_tracks,maxtracks);
//   
  
   DoTheJob(hitsList, HitsInLayers, seeds, temp_tracks, xVtx, yVtx, zVtx, zVtxErr, nEv, zoomlevels, phimin, phimax, SV, dzdlmin, dzdlmax, pTmin, debugOutput,HoughTransformSeeding, SeedRefitterOnly);
  
  for (unsigned int i=0;i<temp_tracks.size();i++)
  {
    if (temp_tracks[i].hits.size()==3)
    {
      unsigned int ind0=temp_tracks[i].hits[0].index();
      unsigned int ind1=temp_tracks[i].hits[1].index();
      unsigned int ind2=temp_tracks[i].hits[2].index();
      
      OrderedHitTriplet triplet (seedHits[ind0],seedHits[ind1],seedHits[ind2]);

      tracks.push_back(triplet);
    }
    
  }
  
  std::cout << "Resulting tracks :"<<tracks.size()<<std::endl;
  
//  finalize(hits,temp_tracks,tracks);

  LogDebug("HTTrackSeeding") << "Time spent in seeding";
  LogTrace("HTTrackSeeding")  << "   vote time = " << voteTime  ().lapse();
  LogTrace("HTTrackSeeding")  << "xy vote time = " << voteTimeXY().lapse();
  LogTrace("HTTrackSeeding")  << " z vote time = " << voteTimeZ ().lapse();
}


void ::HelixHough::findSeededHelices ( std::vector< ::SimpleTrack3D> & seeds     , 
                                     const TrackingRegion::Hits & hits      ,
                                     unsigned int                 min_hits  ,
                                     unsigned int                 max_hits  ,
                                     OrderedHitTriplets         & tracks    ,
                                     unsigned int                 maxtracks )
{
  voteTime  ().reset();
  voteTimeXY().reset();
  voteTimeZ ().reset();

  LogDebug("HTTrackSeeding") << "Hits in input to HelixHough findSeededHelices";
  std::vector< ::SimpleHit3D>  hitsList;
  hitsList.reserve( hits.size() );
  unsigned int i = 0;
  for ( TrackingRegion::Hits::const_iterator hit =hits.begin(); hit!=hits.end(); hit++, i++ ) {
      ::SimpleHit3D ahit ((*hit),_origin,i);
     if ( !ahit.isValid() )  continue;
     LogTrace("HTTrackSeeding") << ahit;
     hitsList.push_back(ahit);
  }
  assert( hitsList.size()==hits.size() );
  if ( hitsList.size()<3 ) return;

  initEvent(hitsList, min_hits);
  initSeeding();
  
  HelixHoughEngine engine(*this,_range,_nBins);
  std::vector< ::SimpleTrack3D> temp_tracks;
  engine.findSeededHelices(seeds,hitsList,min_hits,max_hits,temp_tracks,maxtracks);
  
  finalize(hits,temp_tracks,tracks);

  LogDebug("HTTrackSeeding") << "Time spent in seeding";
  LogTrace("HTTrackSeeding")  << "   vote time = " << voteTime  ().lapse();
  LogTrace("HTTrackSeeding")  << "xy vote time = " << voteTimeXY().lapse();
  LogTrace("HTTrackSeeding")  << " z vote time = " << voteTimeZ ().lapse();
}


bool  ::HelixHough::breakRecursion ( const std::vector< ::SimpleHit3D>   & hits   ,
                                   const HelixParRange              & range  )  const
{ 
  unsigned int  layers = numberOfLayers( hits );
  return layers<_requiredLayers;
}


/*
 *  Dovrei sistemare i parametri e gli indici di traccia, sempre che abbiano un senso
 */

void  ::HelixHough::findTracks ( const std::vector< ::SimpleHit3D>   & hits   ,
                               std::vector< ::SimpleTrack3D>       & tracks ,
                               const HelixParRange              & range  )
{
  //LogDebug("HTTrackSeeding") << "Track seeding";
  std::multimap< ::SimpleHit3D::TrkLayerKey,const  ::SimpleHit3D*> layers;
  for ( auto hit = hits.begin(); hit != hits.end(); hit++ )
    layers.insert( std::pair< ::SimpleHit3D::TrkLayerKey,const  ::SimpleHit3D*>( hit->layer(), &(*hit) ) );
  unsigned int  count = 0;
  for ( auto iter = layers.begin(); iter != layers.end(); iter = layers.equal_range(iter->first).second )
     count++ ;
  if ( count<_requiredLayers ) return;
  std::vector< ::SimpleTrack3D> newTracks;
  /*
  // qui la faccio piu' semplice al momento
  ::SimpleTrack3D track;
  track.curv = range.curv().center();
  track.eta  = range.eta ().center();
  track.lip  = range.lip ().center();
  track.phi  = range.phi ().center();
  track.tip  = range.tip ().center();
  for ( auto hit = hits.begin(); hit != hits.end(); hit++ )
    track.hits.push_back(*hit);
  newTracks.push_back( track );

  LogTrace("HTTrackSeeding") << std::fixed << std::setprecision(4) << std::setfill(' ')
                             << " [" << std::setw(7) << range.curv().lower() << "," << std::setw(7) << range.curv().upper() << "]"
                             << " [" << std::setw(7) << range.eta ().lower() << "," << std::setw(7) << range.eta ().upper() << "]"
                             << " [" << std::setw(8) << range.lip ().lower() << "," << std::setw(8) << range.lip ().upper() << "]"
                             << " [" << std::setw(7) << range.phi ().lower() << "," << std::setw(7) << range.phi ().upper() << "]"
                             << " [" << std::setw(7) << range.tip ().lower() << "," << std::setw(7) << range.tip ().upper() << "]"
			     << " with " << track.hits.size() << " hits";

  // fine modifica 
  */

  // generate ntuples keeping hits in the first _requiredLayers; higher index layers are generally outer
  // don't use two hits from the same layer in each ntuple
  unsigned int nlayers =0;
  for ( auto itk = layers.begin(); itk != layers.end() && nlayers<_requiredLayers ; nlayers++) {
    auto  itv = layers.equal_range(itk->first);
    if ( newTracks.empty() ) {
      ::SimpleTrack3D track;
      track.curv    = range.curv();
      track.eta     = range.eta ();
      track.lip     = range.lip ();
      track.phi     = range.phi ();
      track.tip     = range.tip ();
      track.nlayers = count;
      track.nhits   = hits.size ();
      track.shared  = 1;
      assert( !track.curv.isEmpty() );
      assert( !track.eta .isEmpty() );
      assert( !track.lip .isEmpty() );
      assert( !track.phi .isEmpty() );
      assert( !track.tip .isEmpty() );
      newTracks.push_back( track );
    }
    std::vector< ::SimpleTrack3D> list;
    list.reserve( newTracks.size()*std::distance(itv.first,itv.second) );
    for ( auto hit = itv.first; hit != itv.second; hit++ )
      for ( auto track = newTracks.begin(); track!= newTracks.end(); track++ ) {
	::SimpleTrack3D trk = *track;
        trk.hits.push_back(*(hit->second));
        list.push_back(trk);
      }
    std::swap(newTracks,list);
    itk = itv.second;
  }
  
  
//  tracks.insert(tracks.end(),newTracks.begin(),newTracks.end());
//instead of storing - output them!
//for (unsigned int i=0;i<newTracks.size();i++)
//  std::cout << "track :"<<newTracks[i].curv.center()<<std::endl;

}


void  ::HelixHough::finalize ( const TrackingRegion::Hits & hits, const std::vector< ::SimpleTrack3D> & input, OrderedHitTriplets & output )
{
  LogDebug("HTTrackSeeding") << "Track seeds found "  << input.size();
  /*
  int j = 0;
  for ( auto track = input.begin(); track!= input.end(); track++, j++ ) {
    LogTrace("HTTrackSeeding") << "Track number " << j << " with " << track->hits.size() << " hits";
    for ( auto hit = track->hits.begin(); hit != track->hits.end(); hit++ )
      LogTrace("HTTrackSeeding") << "  " << hit->index();
  }
  */
  
  
  std::map<int,::SimpleTrack3D> mapper;
  for ( auto track = input.begin(); track!= input.end(); track++ ) {
    int code = 0;
    for ( auto hit = track->hits.begin(); hit != track->hits.end(); hit++ )
      code = ( code<<10 ) | hit->index();
    ::SimpleTrack3D atrack=mapper[code];
    atrack.curv = smallestCovering(atrack.curv,track->curv); 
    atrack.eta  = smallestCovering(atrack.eta ,track->eta ); 
    atrack.lip  = smallestCovering(atrack.lip ,track->lip );
    atrack.phi  = smallestCovering(atrack.phi ,track->phi );
    atrack.tip  = smallestCovering(atrack.tip ,track->tip ); 
    assert( !atrack.curv.isEmpty() );
    assert( !atrack.eta .isEmpty() );
    assert( !atrack.lip .isEmpty() );
    assert( !atrack.phi .isEmpty() );
    assert( !atrack.tip .isEmpty() );
    atrack.nlayers = std::max(track->nlayers,atrack.nlayers);
    atrack.nhits   = std::max(track->nhits  ,atrack.nhits  );
    atrack.shared++;
    mapper[code] = atrack;
  }

  LogDebug("HTTrackSeeding") << "Track seeds after duplicate removal "  << mapper.size();
  int i = 0;
  for ( auto iter = mapper.begin(); iter != mapper.end(); iter++, i++ ) {
    int index1 = (iter->first)>>20 & 0x3ff;
    int index2 = (iter->first)>>10 & 0x3ff;
    int index3 = (iter->first)     & 0x3ff;
    output.push_back( OrderedHitTriplet(hits[index1],hits[index2],hits[index3]) );
    LogTrace("HTTrackSeeding") << std::fixed << std::setprecision(4) << std::setfill(' ')
                               << "Track number " << std::setw(4) << i
                               << " -> " << std::setw(3) << index1
                               << " "    << std::setw(3) << index2
                               << " "    << std::setw(4) << index3 << "  "
                               << " [" << std::setw(9) << iter->second.curv.lower() << "," << std::setw(9) << iter->second.curv.upper() << "]"
                               << " [" << std::setw(9) << iter->second.eta .lower() << "," << std::setw(9) << iter->second.eta .upper() << "]"
                               << " [" << std::setw(9) << iter->second.lip .lower() << "," << std::setw(9) << iter->second.lip .upper() << "]"
                               << " [" << std::setw(9) << iter->second.phi .lower() << "," << std::setw(9) << iter->second.phi .upper() << "]"
                               << " [" << std::setw(9) << iter->second.tip .lower() << "," << std::setw(9) << iter->second.tip .upper() << "]"
                               << "  " << std::setw(3) << iter->second.nlayers
                               << " "  << std::setw(3) << iter->second.nhits
                               << " "  << std::setw(5) << iter->second.shared;
  }
}
