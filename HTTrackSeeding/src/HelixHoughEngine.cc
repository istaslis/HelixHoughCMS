#include "MLoVetere/HTTrackSeeding/interface/HelixHoughEngine.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParBinId.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"
#include "MLoVetere/HTTrackSeeding/interface/RangeFinder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>

    
HelixHoughEngine::HelixHoughEngine ( HelixHough & context, HelixParRange & range, HelixParNBins & nbins )
  : HelixHoughProxy(context), HelixHoughEngineBase(range,nbins)
{ }


HelixHoughEngine::HelixHoughEngine ( HelixHoughEngine & other, HelixParRange & range, HelixParNBins & nbins, int level )
  : HelixHoughProxy(other), HelixHoughEngineBase(range,nbins), Level(level)
{ 
  //if (Level==1)
  //  std::cout << Level<<std::endl;
  
}


HelixHoughEngine::~HelixHoughEngine ( )
{ }


bool  HelixHoughEngine::decentResolution ( )  const 
{
  return ( nCurv()==1 || dCurv() <= minimumResolution().dCurv() ) && 
         ( nEta ()==1 || dEta () <= minimumResolution().dEta () ) &&
         ( nLip ()==1 || dLip () <= minimumResolution().dLip () ) &&
         ( nPhi ()==1 || dPhi () <= minimumResolution().dPhi () ) && 
         ( nTip ()==1 || dTip () <= minimumResolution().dTip () ); 
}
    

bool  HelixHoughEngine::insaneResolution ( )  const 
{
  return ( nCurv()==2 ) && 
         ( nEta ()==2 ) &&
         ( nLip ()==2 ) &&
         ( nPhi ()==2 ) && 
         ( nTip ()==2 );
}


void HelixHoughEngine::findHelices ( const std::vector<SimpleHit3D> & hits, unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks, bool forcezoom )
{
  
  if ( ((maxtracks!=0) && (tracks.size()>=maxtracks))) return;
  unsigned int tracks_at_start = tracks.size();
  voteTime().start();
  vote(hits);
  voteTime().stop();


 //std::cout << "NextLevel Binning! tracks :"<<tracks.size()<<" hits :"<<hits.size()<<" nCurv:"<<nCurv()<<" bins: "<<bins_vec.size()<<" Level :"<<Level<<std::endl;


  // compute next level binning
  HelixParNBins nextBins ( ( dCurv()>maximumResolution().dCurv() ) ? std::min( nCurv(), static_cast<unsigned int>(ceil( dCurv()/maximumResolution().dCurv())) ) : 1 ,
                           ( dEta ()>maximumResolution().dEta () ) ? std::min( nEta (), static_cast<unsigned int>(ceil( dEta ()/maximumResolution().dEta ())) ) : 1 ,
                           ( dLip ()>maximumResolution().dLip () ) ? std::min( nLip (), static_cast<unsigned int>(ceil( dLip ()/maximumResolution().dLip ())) ) : 1 ,
                           ( dPhi ()>maximumResolution().dPhi () ) ? std::min( nPhi (), static_cast<unsigned int>(ceil( dPhi ()/maximumResolution().dPhi ())) ) : 1 ,
                           ( dTip ()>maximumResolution().dTip () ) ? std::min( nTip (), static_cast<unsigned int>(ceil( dTip ()/maximumResolution().dTip ())) ) : 1 );
  
  
  std::unordered_map<HelixParBinId,int> distinct;
int distbins = 0;
  
  // scan over the bins in 5-D hough space
  for ( auto itk = bins_vec.begin(); itk != bins_vec.end();) {
    const auto  bin = itk->first;

    if (distinct[bin]>0) 
    {continue;
     std::cout << "continue!";
    }
      else distinct[bin]=1;
  
    
    auto  itv = bins_vec.equal_range(bin);
    HelixParRange  nextRange = range(bin);
    HelixHoughEngine * nextLevel = new HelixHoughEngine( *this, nextRange, nextBins, Level+1 );
    if (Level==0) {
        distbins++;
	std::cout << "bin # "<<distbins<<std::endl;
    }
    
    
    assert( nextLevel!=0 );
    std::vector<SimpleHit3D> nextHits;
    nextHits.reserve( std::distance(itv.first,itv.second) );
    for ( auto iter = itv.first; iter != itv.second; iter++ )
      nextHits.push_back( hits[ iter->second ] );    
    if ( nextHits.size()>=min_hits ) {
      unsigned int expected_hits = static_cast<unsigned int>( decreasePerZoom()*hits.size() );
      bool zooming_helps =  forcezoom || nextHits.size()<=expected_hits;
      if ( nextLevel->insaneResolution() || breakRecursion( nextHits, nextLevel->range() ) || Level>1    ) {
        nextLevel->findTracks( nextHits, tracks, nextLevel->range() );
      } else if ( nextLevel->decentResolution() && ( nextHits.size()<=max_hits || !zooming_helps ) ) {
        nextLevel->findTracks( nextHits, tracks, nextLevel->range() );
      } else {
        nextLevel->findHelices( nextHits, min_hits, max_hits, tracks, maxtracks, false );
      }
      if ( maxtracks!=0 ) {
        double curv_proportion = DCurv() / topRange().DCurv();
        double  phi_proportion = DPhi()  / topRange().DPhi() ;
        double  tip_proportion = DTip()  / topRange().DTip() ;
        unsigned int expected_tracks = (unsigned int) fabs( ((double)maxtracks)*curv_proportion*phi_proportion*tip_proportion ) +1;
        if ( (tracks.size()-tracks_at_start)>expected_tracks ) return;
      }
    }
    delete nextLevel;
    itk = itv.second;

   
  }
  
 // std::cout << "bins : "<<distbins<<" hits :"<<hits.size()<<" Level :"<<Level<<std::endl;


  
}


void HelixHoughEngine::findSeededHelices ( const std::vector<SimpleTrack3D> & seeds     ,
                                           const std::vector<SimpleHit3D>   & hits      , 
                                           unsigned int                       min_hits  ,
                                           unsigned int                       max_hits  ,
                                           std::vector<SimpleTrack3D>       & tracks    ,
                                           unsigned int                       maxtracks ,
                                           bool                               forcezoom )
{
  if ( (maxtracks!=0) && (tracks.size()>=maxtracks) ) return;
  unsigned int tracks_at_start = tracks.size();
  voteTime().start();
  vote(hits);
  voteTime().stop();
  std::unordered_multimap<HelixParBinId,const SimpleTrack3D*> seedsMap;
  map2(seeds,seedsMap);
  // compute next level binning
  HelixParNBins nextBins ( ( dCurv()>maximumResolution().dCurv() ) ? std::min( nCurv(), static_cast<unsigned int>(ceil( dCurv()/maximumResolution().dCurv())) ) : 1 ,
                           ( dEta ()>maximumResolution().dEta () ) ? std::min( nEta (), static_cast<unsigned int>(ceil( dEta ()/maximumResolution().dEta ())) ) : 1 ,
                           ( dLip ()>maximumResolution().dLip () ) ? std::min( nLip (), static_cast<unsigned int>(ceil( dLip ()/maximumResolution().dLip ())) ) : 1 ,
                           ( dPhi ()>maximumResolution().dPhi () ) ? std::min( nPhi (), static_cast<unsigned int>(ceil( dPhi ()/maximumResolution().dPhi ())) ) : 1 ,
                           ( dTip ()>maximumResolution().dTip () ) ? std::min( nTip (), static_cast<unsigned int>(ceil( dTip ()/maximumResolution().dTip ())) ) : 1 );
  // scan over the bins in 5-D hough space
  for ( auto itk = bins_vec.begin(); itk != bins_vec.end(); ) {
    const auto  bin = itk->first;
    auto  itv = bins_vec.equal_range(bin);
    HelixParRange  nextRange = range(bin);
    HelixHoughEngine * nextLevel = new HelixHoughEngine( *this, nextRange, nextBins, Level+1 );
    assert( nextLevel==0 );
    std::vector<SimpleTrack3D> nextSeeds;
    auto ret = seedsMap.equal_range(bin);
    nextSeeds.reserve( std::distance(ret.first,ret.second) );
    for ( auto iter = ret.first; iter != ret.second; iter++ )
      nextSeeds.push_back( *iter->second );
    std::vector<SimpleHit3D> nextHits;
    nextHits.reserve( std::distance(itv.first,itv.second) );
    for ( auto iter = itv.first; iter != itv.second; iter++ )
      nextHits.push_back( hits[ iter->second ] );    
    if ( (nextHits.size()>=min_hits) && (nextSeeds.size()!=0) ) {
      unsigned int expected_hits = static_cast<unsigned int>( decreasePerZoom()*hits.size() );
      bool zooming_helps =  forcezoom || nextHits.size()<=expected_hits;
      if ( nextLevel->insaneResolution() || breakRecursion( nextHits, nextLevel->range() ) ) {
        nextLevel->findSeededTracks ( nextSeeds, nextHits, tracks, nextLevel->range() );
      } else if ( nextLevel->decentResolution() && ( nextHits.size()<=max_hits || !zooming_helps ) ) {
        nextLevel->findSeededTracks ( nextSeeds, nextHits, tracks, nextLevel->range() );
      } else {
        nextLevel->findSeededHelices( nextSeeds, nextHits, min_hits, max_hits, tracks, maxtracks, false );
      }
      if ( maxtracks!=0 ) {
        double curv_proportion = DCurv() / topRange().DCurv();
        double  phi_proportion = DPhi()  / topRange().DPhi() ;
        double  tip_proportion = DTip()  / topRange().DTip() ;
        unsigned int expected_tracks = (unsigned int) fabs( ((double)maxtracks)*curv_proportion*phi_proportion*tip_proportion ) +1;
        if ( (tracks.size()-tracks_at_start)>expected_tracks ) return;
      }
    }
    delete nextLevel;
    itk = itv.second;
  }
}


void HelixHoughEngine::vote ( const std::vector<SimpleHit3D> & hits )
{
  voteTimeXY().start();
  voteTimeZ ().start();
  bins_vec    .clear();
  float curv_size = dCurv();
  float  lip_size = dLip ();
  float  tip_size = dTip ();
  for ( unsigned int i=0; i<hits.size(); i++ ) {
    Interval curv( range().minCurv(), range().minCurv()+curv_size );
    for ( unsigned int binCurv=0; binCurv<nCurv(); ++binCurv, curv.shift(curv_size) ) {
      Interval tip( range().minTip(), range().minTip()+tip_size );
      for ( unsigned int binTip=0; binTip<nTip(); ++binTip, tip.shift(tip_size) ) {
        RangeFinderExact finder(hits[i],curv,tip);
        if ( finder.dPhi().isEmpty() ) continue;
  	for ( unsigned int hturn=-inHalfTurns(); hturn<outHalfTurns(); hturn++ ) {
          std::pair<BinRange,BinRange>  phiNRange = phi2bin(finder.dPhi(hturn));
          /*
          LogTrace("HTTrackSeeding") << std::fixed << std::setprecision(4) << std::setfill(' ')
                                     << "Track no " << std::setw(3) << hits[i].index() << "  "
                                     << " " << std::setw(7) << finder.dPhi().lower() << " " << std::setw(7) << finder.dPhi().upper()
                                     << " [" << phiNRange.first .first << "," << phiNRange.first .second << "]"
                                     << " [" << phiNRange.second.first << "," << phiNRange.second.second << "]";
	  */
          Interval lip( range().minLip(), range().minLip()+lip_size );
	  for ( unsigned int binLip=0; binLip<nLip(); ++binLip, lip.shift(lip_size) ) {
            BinRange  etaNRange = eta2bin(finder.dEta(lip,hturn));
            /*
            LogTrace("HTTrackSeeding") << std::fixed << std::setprecision(4) << std::setfill(' ')
                                       << "Track no " << std::setw(3) << hits[i].index() << "  "
                                       << "*" << std::setw(9) << lip.lower() << " " << std::setw(9) << lip.upper()
                                       << " " << std::setw(7) << finder.dEta(lip,hturn).lower() << " " << std::setw(7) << finder.dEta(lip,hturn).upper()
                                       << " [" << etaNRange.first << "," << etaNRange.second << "]";
	    */
            for ( unsigned int binEta = etaNRange.first; binEta<etaNRange.second; binEta++ ) {
              for ( unsigned int binPhi = phiNRange.first.first; binPhi<phiNRange.first.second; binPhi++ )
                bins_vec.insert( std::pair<HelixParBinId,unsigned int>( HelixParBinId(binCurv,binEta,binLip,binPhi,binTip), i ) );
              for ( unsigned int binPhi = phiNRange.second.first; binPhi<phiNRange.second.second; binPhi++ )
                bins_vec.insert( std::pair<HelixParBinId,unsigned int>( HelixParBinId(binCurv,binEta,binLip,binPhi,binTip), i ) );
            }
          }
        }
      }
    }
  }
  voteTimeXY().stop();
  voteTimeZ ().stop();
  /*
  LogDebug("HTTrackSeeding") << "Voting matrix";
  for ( unsigned int binCurv =0; binCurv<nCurv(); binCurv++ )
    for ( unsigned int binEta =0; binEta<nEta(); binEta++ )
      for ( unsigned int binLip =0; binLip<nLip(); binLip++ )
        for ( unsigned int binPhi =0; binPhi<nPhi(); binPhi++ )
          for ( unsigned int binTip =0; binTip<nTip(); binTip++ ) {
	    HelixParBinId bin(binCurv,binEta,binLip,binPhi,binTip);
            int count = bins_vec.count(bin);
            if ( count>2 )  LogTrace("HTTrackSeeding") << std::setfill(' ')
                                                       << " " << std::setw(2) << binCurv
                                                       << " " << std::setw(2) << binEta
                                                       << " " << std::setw(2) << binLip 
                                                       << " " << std::setw(2) << binPhi
                                                       << " " << std::setw(2) << binTip
                                                       << "  -> " << std::setw(5) << count;
          }
  */
}
