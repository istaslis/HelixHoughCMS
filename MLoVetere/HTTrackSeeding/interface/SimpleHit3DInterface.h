#ifndef HTTrackSeeding_SimpleHit3D_H
#define HTTrackSeeding_SimpleHit3D_H

/*** \class  SimpleHit3D
  *
  *  WARNING. This class has to go away quickly!
  *
  *  \author Maurizio Lo Vetere
  */

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <utility>
#include <vector>
#include <iostream>

class SimpleHit3D
{
  public:
    typedef int TrkLayerKey;
    SimpleHit3D ( const TrackingRegion::Hit & hit,  GlobalPoint origin, unsigned int index =0 );
    Interval           rho   ( )  const  { return Interval();     }
    AngularInterval    phi   ( )  const  { return AngularInterval();     }
    Interval           x     ( )  const  { return _x;       }
    Interval           y     ( )  const  { return _y;       }
    Interval           z     ( )  const  { return _z;       }
    unsigned int     index   ( )  const  { return _index;   }
    bool             isValid ( )  const  { return _isValid; }
    TrkLayerKey      layer   ( )  const  { return _layer;   }
    void             print   ( std::ostream&    s )  const;
  private:
    bool             _isValid;
    //Interval         _rho;
    //AngularInterval  _phi;
    Interval         _x;
    Interval         _y;
    Interval         _z;
    TrkLayerKey       _layer;
    unsigned int     _index;
};


inline std::ostream& operator<< ( std::ostream& s, const SimpleHit3D hit )
{ 
  hit.print(s);
  return s;
}


inline unsigned int numberOfLayers ( const std::vector<SimpleHit3D> & hits )
{ 
  std::set<SimpleHit3D::TrkLayerKey> layers;
  for ( std::vector<SimpleHit3D>::const_iterator hit = hits.begin();  hit != hits.end(); hit++ ) 
    layers.insert( hit->layer() );
  return layers.size();
}


inline SimpleHit3D::SimpleHit3D ( const TrackingRegion::Hit & hit,  GlobalPoint origin, unsigned int index ) 
  : _isValid(false), _index(index)  
{
  DetId id = hit->geographicalId();
  if ( id.det() != DetId::Tracker )  {
    edm::LogWarning("HTTrackSeeding") << "Cannot generate a valid SimpleHit3D - detId doesn't belong to the tracker";
    return;
  };
  switch ( id.subdetId() ) { 
    case PixelSubdetector::PixelBarrel:
      _layer = PXBDetId(id).layer();
      break;
    case PixelSubdetector::PixelEndcap:
      _layer = PXFDetId(id).disk();
      break;
    case StripSubdetector::TIB:
      _layer = TIBDetId(id).layer();
      break;
    case  StripSubdetector::TID:
      _layer = TIDDetId(id).wheel();
      break;
    case StripSubdetector::TOB:
      _layer = TOBDetId(id).layer();
      break;
    case StripSubdetector::TEC:
      _layer = TECDetId(id).wheel();
      break;
    default:
      edm::LogWarning("HTTrackSeeding") << "Cannot generate a valid SimpleHit3D - subdetId doesn't belong to the tracker";
      return;      
  }
  _layer |= (id.subdetId() << 4);
  _isValid = hit->isValid();
  if ( !_isValid ) {
    edm::LogWarning("HTTrackSeeding") << "Cannot generate a valid SimpleHit3D - invalid hit";
    return;
  }
    Vector3DBase<float,GlobalTag> p = hit->globalPosition()-GlobalPoint(0.,0.,0.);//-origin; - IMPORTANT! we don't need additional shift in HI since we have ONE vertex
//  double    rho = p.perp();
//  double    phi = p.phi ();                                     // could be p.barePhi();
   //std::cout << "hit Id "<<(int)id<<" r "<<(double)p.perp()<<std::endl;
  double    x   = p.x();
  double    y   = p.y();
  double    z   = p.z   ();
  GlobalError c = hit->globalPositionError();
//  double   drho = sqrt( c.rerr  ( p+GlobalPoint(0.,0.,0.) ) );  // I'm cheating because p is not in the global reference but it is ok to get errors.
//  double   dphi = sqrt( c.phierr ( p+GlobalPoint(0.,0.,0.) ) );  // I'm cheating because p is not in the global reference but it is ok to get errors.
  double   dx   = sqrt( c.cxx())/2*3.14;
  double   dy   = sqrt( c.cyy())/2*3.14;
  double   dz   = sqrt( c.czz())/2*3.14;
  // Workaround to cope with missing errors case
//  drho = std::max(drho,0.00001);
//  dphi = std::max(dphi,0.00001);
  dx = std::max(dx,0.00001);
  dy = std::max(dy,0.00001);
  dz   = std::max(dz  ,0.00001);
//  _rho = Interval( std::max(0.,rho-drho), rho+drho );
//  _phi = AngularInterval( phi-dphi, phi+dphi );
  _x   = Interval( x-dx, x+dx );
  _y   = Interval( y-dy, y+dy );
  _z   = Interval( z-dz, z+dz );
}


inline void  SimpleHit3D::print ( std::ostream& s )  const
{
  if ( _isValid ) {
    s << std::fixed << std::setprecision(4) << std::setfill(' ')
 //     <<   "r= [" << std::setw(7) << _rho.lower() << " cm, "  << std::setw(7) << _rho.upper() << " cm], "
 //     << "phi= [" << std::setw(7) << _phi.lower() << " rad, " << std::setw(7) << _phi.upper() << " rad], "
 //     <<   "z= [" << std::setw(9) <<   _z.lower() << " cm, "  << std::setw(9) <<   _z.upper() << " cm] ";
            ;
  } else
    s << "invalid hit";
}


#endif // HTTrackSeeding_SimpleHit3D_H
