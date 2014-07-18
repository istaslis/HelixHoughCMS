#ifndef HTTrackSeeding_HelixParRange_H
#define HTTrackSeeding_HelixParRange_H

/*** \class  HelixParRange
  *
  *  Helix segment parameters, including limits on positive and negative half turns 
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"

#include <algorithm>
#include <utility>


class HelixParRange
{
  public:
    HelixParRange ( ) { }
    HelixParRange ( Interval curv, Interval eta, Interval lip , AngularInterval phi, Interval tip,
                    unsigned int pHTurns =1, unsigned int nHTurns =0 )
      : _rCurv(curv), _rEta(eta), _rLip(lip), _rPhi(phi), _rTip(tip), _pHTurns(pHTurns), _nHTurns(nHTurns)
    { }
    HelixParRange ( float minCurv, float maxCurv,
                    float minEta , float maxEta ,
                    float minLip , float maxLip ,
                    float minPhi , float maxPhi ,
                    float minTip , float maxTip ,
                    unsigned int pHTurns =1     ,
                    unsigned int nHTurns =0     )
      : _rCurv(minCurv,maxCurv),
        _rEta (minEta ,maxEta ),
        _rLip (minLip ,maxLip ),
        _rPhi (minPhi ,maxPhi ),
        _rTip (minTip ,maxTip ),
        _pHTurns(pHTurns)      ,
        _nHTurns(nHTurns)
    { }
   ~HelixParRange ( ) { }
    float  DCurv   ( )  const  { return _rCurv.length(); }
    float  DEta    ( )  const  { return _rEta .length(); }
    float  DLip    ( )  const  { return _rLip .length(); }
    float  DPhi    ( )  const  { return _rPhi .length(); }
    float  DTip    ( )  const  { return _rTip .length(); }
    float  maxCurv ( )  const  { return _rCurv.upper();  }
    float  maxEta  ( )  const  { return _rEta .upper();  }
    float  maxLip  ( )  const  { return _rLip .upper();  }
    float  maxPhi  ( )  const  { return _rPhi .upper();  } 
    float  maxTip  ( )  const  { return _rTip .upper();  }
    float  minCurv ( )  const  { return _rCurv.lower();  }
    float  minEta  ( )  const  { return _rEta .lower();  }
    float  minLip  ( )  const  { return _rLip .lower();  }
    float  minPhi  ( )  const  { return _rPhi .lower();  }
    float  minTip  ( )  const  { return _rTip .lower();  }
    float                         rCurv   ( float           value )  const  { return _rCurv.normalize(value); }
    Interval                      rCurv   ( Interval        value )  const  { return _rCurv.normalize(value); }
    float                         rEta    ( float           value )  const  { return _rEta .normalize(value); }
    Interval                      rEta    ( Interval        value )  const  { return _rEta .normalize(value); }
    float                         rLip    ( float           value )  const  { return _rLip .normalize(value); }
    Interval                      rLip    ( Interval        value )  const  { return _rLip .normalize(value); }
    float                         rPhi    ( float           value )  const  { return _rPhi .normalize(value); }
    std::pair<Interval,Interval>  rPhi    ( AngularInterval value )  const  { return _rPhi .normalize(value); }
    float                         rTip    ( float           value )  const  { return _rTip .normalize(value); }
    Interval                      rTip    ( Interval        value )  const  { return _rTip .normalize(value); }
    float                         pHTurns ( )                        const  { return _pHTurns;                }
    float                         nHTurns ( )                        const  { return _nHTurns;                }
  public:
    Interval         curv ( )                  const  { return _rCurv; }
    Interval         eta  ( )                  const  { return _rEta;  }
    Interval         lip  ( )                  const  { return _rLip;  }
    AngularInterval  phi  ( )                  const  { return _rPhi;  }
    Interval         tip  ( )                  const  { return _rTip;  }
 // deprecated
    void          setCurv ( float min, float delta )  { _rCurv =        Interval( min, min+delta ); }
    void          setEta  ( float min, float delta )  { _rEta  =        Interval( min, min+delta ); }
    void          setLip  ( float min, float delta )  { _rLip  =        Interval( min, min+delta ); }
    void          setPhi  ( float min, float delta )  { _rPhi  = AngularInterval( min, min+delta ); }
    void          setTip  ( float min, float delta )  { _rTip  =        Interval( min, min+delta ); }
  private:
    Interval        _rCurv;
    Interval        _rEta;
    Interval        _rLip;
    AngularInterval _rPhi;
    Interval        _rTip;
    unsigned int    _pHTurns;
    unsigned int    _nHTurns;
};


#endif // HTTrackSeeding_HelixParRange_H
