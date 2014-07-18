#ifndef HTTrackSeeding_HelixParResolution_H
#define HTTrackSeeding_HelixParResolution_H

/*** \class  HelixParResolution
  *
  *  Requested resolution on Helix parameters
  *
  *  \author Maurizio Lo Vetere
  */


class HelixParResolution
{
  public:
    HelixParResolution ( )  { }
    HelixParResolution ( float dCurv, float dEta, float dLip, float dPhi, float dTip ) 
      : _dCurv(dCurv), _dEta(dEta), _dLip(dLip), _dPhi(dPhi), _dTip(dTip) { }
    float  dCurv ()  const        { return  _dCurv; }
    float  dEta  ()  const        { return  _dEta;  }
    float  dLip  ()  const        { return  _dLip;  }
    float  dPhi  ()  const        { return  _dPhi;  }
    float  dTip  ()  const        { return  _dTip;  }
    void   dCurv ( float value )  { _dCurv = value; }
    void   dEta  ( float value )  { _dEta  = value; }
    void   dLip  ( float value )  { _dLip  = value; }
    void   dPhi  ( float value )  { _dPhi  = value; }
    void   dTip  ( float value )  { _dTip  = value; }
  private:  
    float  _dCurv;
    float  _dEta;
    float  _dLip;
    float  _dPhi;
    float  _dTip;
};


#endif // HTTrackSeeding_HelixParResolution_H
