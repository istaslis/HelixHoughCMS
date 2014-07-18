#ifndef __SIMPLEHIT3D__
#define __SIMPLEHIT3D__

#include <iostream>


namespace PHENIXHough {


class SimpleHit3D
{

 public:
  
   SimpleHit3D(float xx=0., float dxx=0., float yy=0., float dyy=0., float zz=0., float dzz=0., unsigned int ind=0, int lyr=-1) : x(xx), dx(dxx), y(yy), dy(dyy), z(zz), dz(dzz), index(ind), layer(lyr) {}
  ~SimpleHit3D(){};
    
  float x, dx;
  float y, dy;
  float z, dz;
  unsigned int index;
  int layer;
};

}
#endif // __SIMPLEHIT3D__
