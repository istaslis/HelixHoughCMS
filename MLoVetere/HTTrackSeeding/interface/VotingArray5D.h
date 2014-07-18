#ifndef HTTrackSeeding_VotingArray5D_H
#define HTTrackSeeding_VotingArray5D_H

/*** \class  VotingArray5D
  *
  *  An implementation of the voting array. Perhaps not the best one.
  *
  *  \author Maurizio Lo Vetere
  */

#include <cassert>
#include <vector>

#include "MLoVetere/HTTrackSeeding/interface/HelixHoughEngineBase.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParBinId.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParNBins.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParRange.h"


template <typename T, T V>
class  VotingArray5D : public HelixHoughEngineBase, private std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > >
{
  public:
    VotingArray5D ( HelixParRange range, HelixParNBins nbins );
    T &   operator[] ( HelixParBinId bin );
    void  clear      ( );
  private:
    static std::vector<std::vector<std::vector<std::vector<T> > > >  zeroVector ( HelixParNBins nbins ) {
      std::vector<T>                                           vec1( nbins.nTip(), V    );
      std::vector<std::vector<T> >                             vec2( nbins.nPhi(), vec1 );
      std::vector<std::vector<std::vector<T> > >               vec3( nbins.nLip(), vec2 );
      std::vector<std::vector<std::vector<std::vector<T> > > > vec4( nbins.nEta(), vec3 );
      return vec4;
    }
};


template <typename T, T V> inline  VotingArray5D<T,V>::VotingArray5D ( HelixParRange range, HelixParNBins nbins )
  : HelixHoughEngineBase( range, nbins ), 
    std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > > ( nbins.nCurv(), zeroVector( nbins ) ) 
{
  assert( nbins.nCurv()>0 );
  assert( nbins.nEta ()>0 );
  assert( nbins.nLip ()>0 );
  assert( nbins.nPhi ()>0 );
  assert( nbins.nTip ()>0 );
}


template <typename T, T V> inline T &  VotingArray5D<T,V>::operator[] ( HelixParBinId bin )
{
  return (*this)[bin.nCurv()][bin.nEta()][bin.nLip()][bin.nPhi()][bin.nTip()];
}


template <typename T, T V> inline void  VotingArray5D<T,V>::clear ( )
{
  for ( unsigned int i=0; i<nCurv(); i++ )
    for ( unsigned int j=0; j<nEta(); j++ )
      for ( unsigned int k=0; k<nLip(); k++ )
        for ( unsigned int l=0; l<nPhi(); l++ )
          for ( unsigned int m=0; m<nTip(); m++ )
            (*this)[i][j][k][l][m] = 0;
}


#endif // HTTrackSeeding_VotingArray5D_H
