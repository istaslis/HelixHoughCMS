#ifndef HTTrackSeeding_HelixHoughEngine_H
#define HTTrackSeeding_HelixHoughEngine_H

/*** \class  HelixHoughEngine
  *
  *  This class implement seed search using a Hough transform inspired algorithm
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/HelixHoughInterface.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixHoughProxy.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixHoughEngineBase.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParBinId.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParNBins.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParRange.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParResolution.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTimer.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"
#include "MLoVetere/HTTrackSeeding/interface/VotingArray5D.h"

#include <unordered_map>
#include <vector>
#include <utility>


class HelixHoughEngine : protected HelixHoughProxy, protected HelixHoughEngineBase
{
  public:
    HelixHoughEngine           ( HelixHough       & context ,
                                 HelixParRange    & range   ,
                                 HelixParNBins    & nbins   );
    HelixHoughEngine           ( HelixHoughEngine & other   ,
                                 HelixParRange    & range   ,
                                 HelixParNBins    & nbins   ,
				 int 		    level
 			      );
    virtual  ~HelixHoughEngine ( );
    void  findHelices          ( const std::vector<SimpleHit3D>   & hits            ,
                                 unsigned int                       min_hits        ,
                                 unsigned int                       max_hits        ,
                                 std::vector<SimpleTrack3D>       & tracks          ,
                                 unsigned int                       maxtracks =0    ,
                                 bool                               forcezoom =true );
    void  findSeededHelices    ( const std::vector<SimpleTrack3D> & seeds           ,
                                 const std::vector<SimpleHit3D>   & hits            ,
                                 unsigned int                       min_hits        ,
                                 unsigned int                       max_hits        ,
                                 std::vector<SimpleTrack3D>       & tracks          ,
                                 unsigned int                       maxtracks =0    ,
                                 bool                               forcezoom =true );
  
    int Level;
  private:
    bool  decentResolution     ( )  const;
    bool  insaneResolution     ( )  const;
    void  fillBins             ( unsigned int                       index,
                                 unsigned int                       curv,
                                 Interval                           eta,
                                 unsigned int                       lip,
                                 AngularInterval                    phi,
                                 unsigned int                       tip             );
    void  vote                 ( const std::vector<SimpleHit3D> &   hits            );
  private:
    std::unordered_multimap<HelixParBinId,unsigned int>  bins_vec;
};


#endif // HTTrackSeeding_HelixHoughEngine_H
