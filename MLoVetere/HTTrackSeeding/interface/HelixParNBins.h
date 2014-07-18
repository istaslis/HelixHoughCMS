#ifndef HTTrackSeeding_HelixParNBins_H
#define HTTrackSeeding_HelixParNBins_H

/*** \class  HelixParNBins
  *
  *  Almost useless. Introduced not to break convention: #include "...<classname>"
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/HelixParBinId.h"

#include <cassert>


class HelixParNBins : private HelixParBinId
{
  public:
    HelixParNBins ( )   { }
    HelixParNBins (unsigned int ncurv, unsigned int neta, unsigned int nlip, unsigned int nphi, unsigned int ntip )
      : HelixParBinId() { nCurv(ncurv); nEta(neta); nLip(nlip); nPhi(nphi); nTip(ntip); }
    unsigned int  nCurv ( )  const              { return HelixParBinId::nCurv()+1; }
    unsigned int  nEta  ( )  const              { return HelixParBinId::nEta ()+1; }
    unsigned int  nLip  ( )  const              { return HelixParBinId::nLip ()+1; }
    unsigned int  nPhi  ( )  const              { return HelixParBinId::nPhi ()+1; }
    unsigned int  nTip  ( )  const              { return HelixParBinId::nTip ()+1; }
    void          nCurv ( unsigned int value )  { assert(value!=0); HelixParBinId::nCurv(value-1); }
    void          nEta  ( unsigned int value )  { assert(value!=0); HelixParBinId::nEta (value-1); }
    void          nLip  ( unsigned int value )  { assert(value!=0); HelixParBinId::nLip (value-1); }
    void          nPhi  ( unsigned int value )  { assert(value!=0); HelixParBinId::nPhi (value-1); }
    void          nTip  ( unsigned int value )  { assert(value!=0); HelixParBinId::nTip (value-1); }
};


#endif // HTTrackSeeding_HelixParBinId_H
