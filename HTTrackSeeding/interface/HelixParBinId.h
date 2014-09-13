#ifndef HTTrackSeeding_HelixParBinId_H
#define HTTrackSeeding_HelixParBinId_H

/*** \class  HelixParBinId
  *
  *  Identifier of the 5D parameter space bin
  *
  *  \author Maurizio Lo Vetere
  */

#include <cassert>
#include <cstddef>
#include <functional>


class HelixParBinId
{
  public:
    HelixParBinId ( ) : _index(0)  { }
    HelixParBinId ( unsigned int ncurv, unsigned int neta, unsigned int nlip, unsigned int nphi, unsigned int ntip )
      : _index(0) { nCurv(ncurv); nEta(neta); nLip(nlip); nPhi(nphi); nTip(ntip); }
    bool          operator==( const HelixParBinId & other )  const  { return  _index == other._index; }
    bool          operator!=( const HelixParBinId & other )  const  { return  _index != other._index; }
    unsigned int  nCurv     ( )                              const  { return get(binC_s,binC_m); }
    unsigned int  nEta      ( )                              const  { return get(binE_s,binE_m); }
    unsigned int  nLip      ( )                              const  { return get(binL_s,binL_m); }
    unsigned int  nPhi      ( )                              const  { return get(binP_s,binP_m); }
    unsigned int  nTip      ( )                              const  { return get(binT_s,binT_m); }
    void          nCurv     ( unsigned int value )                  { set(binC_s,binC_m,value);  };
    void          nEta      ( unsigned int value )                  { set(binE_s,binE_m,value);  };
    void          nLip      ( unsigned int value )                  { set(binL_s,binL_m,value);  };
    void          nPhi      ( unsigned int value )                  { set(binP_s,binP_m,value);  };
    void          nTip      ( unsigned int value )                  { set(binT_s,binT_m,value);  };
  private:
    enum { binC_s =   24, binE_s =   18, binL_s  =   12, binP_s  =   6, binT_s =    0 };
    enum { binC_m = 0x3f, binE_m = 0x3f, binL_m  = 0x3f, binP_m = 0x3f, binT_m = 0x3f };
    unsigned int  get       ( unsigned int shift, unsigned int mask )               const;
    void          set       ( unsigned int shift, unsigned int mask, unsigned int value );
  private:
    friend struct std::hash<HelixParBinId>;
    unsigned int _index;
};


namespace std { template<> struct hash<HelixParBinId> { std::size_t operator()( const HelixParBinId & c ) const { return hash<unsigned int>()(c._index); } }; }


inline unsigned int  HelixParBinId::get ( unsigned int shift, unsigned int mask )  const
{ 
  return (_index>>shift) & mask;
}


inline void  HelixParBinId::set ( unsigned int shift, unsigned int mask, unsigned int value )
{ 
  assert ( value<=mask );
  _index &= ~( mask << shift );
  _index |=  ( (value & mask) << shift );
}


#endif // HTTrackSeeding_HelixParBinId_H
