#include "MLoVetere/HTTrackSeeding/interface/sPHENIXTracker.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sys/time.h>


using namespace std;
using namespace Eigen;
//using namespace SeamStress;


static inline double sign(double x)
{
  return ((double)(x > 0.)) - ((double)(x < 0.));
}

void sPHENIXTracker::projectToDisk(SimpleTrack3D& seed, float z, float& x, float& y)
{
  float phi = seed.phi;
  float d = seed.d;
  float k = seed.kappa;
  float z0 = seed.z0;
  float dzdl = seed.dzdl;

  float r = 1./k;
  float y0 = d*sin(phi);

  float chi = (z-z0)*k*sqrt(1-dzdl*dzdl)/dzdl;
  float chi0= asin(((d+r)*sin(phi)-y0)/r);

  
  x = (d+r)*cos(phi)-r*cos(chi-chi0);
  y = (d+r)*sin(phi)+r*sin(chi-chi0);

  cout << " chi :"<<chi<<" chi0 "<<chi0<<" x "<<x<<" y "<<y<<" z "<<z<<" R "<<sqrt(x*x+y*y)<<" R3d "<<sqrt(x*x+y*y+z*z)<<endl;
}

void sPHENIXTracker::projectToLayer(SimpleTrack3D& seed, unsigned int layer, float& x, float& y, float& z)
{
    vector<float> radii { 4.1, 7.5, 9.92, /*25.0, 34.0, 43.0, 52.0, */61.0, 69.6, 78.2, 86.8, 96.5, 108.0 };
    
  float phi = seed.phi;
  float d = seed.d;
  float k = seed.kappa;
  float z0 = seed.z0;
  float dzdl = seed.dzdl;
  
  float hitx = seed.hits.back().x;
  float hity = seed.hits.back().y;
  
  float rad_det = radii[layer];
  
  float cosphi = cos(phi);
  float sinphi = sin(phi);
  
  k = fabs(k);
  
  float kd = (d*k + 1.);
  float kcx = kd*cosphi;
  float kcy = kd*sinphi;
  float kd_inv = 1./kd;
  float R2 = rad_det*rad_det;
  float a = 0.5*(k*R2 + ( d*d*k + 2.*d ))*kd_inv;
  float tmp1 = a*kd_inv;
  float P2x = kcx*tmp1;
  float P2y = kcy*tmp1;
  
  float h = sqrt(R2 - a*a);
  
  float ux = -kcy*kd_inv;
  float uy = kcx*kd_inv;
  
  float x1 = P2x + ux*h;
  float y1 = P2y + uy*h;
  float x2 = P2x - ux*h;
  float y2 = P2y - uy*h;
  float diff1 = (x1-hitx)*(x1-hitx) + (y1-hity)*(y1-hity);
  float diff2 = (x2-hitx)*(x2-hitx) + (y2-hity)*(y2-hity);
  float signk = 0.;
  if(diff1 < diff2){signk = 1.;}
  else{signk = -1.;}
  x = P2x + signk*ux*h;
  y = P2y + signk*uy*h;
  
  double sign_dzdl = sign(dzdl);
 // double onedzdl2_inv = 1./(1. - dzdl*dzdl);
  double startx = d*cosphi;
  double starty = d*sinphi;
  double D = sqrt((startx-x)*(startx-x) + (starty-y)*(starty-y));
 // double D_inv = 1./D;
  double v = 0.5*k*D;
  z = 0.;
  if(v > 0.1)
  {
    if(v >= 0.999999){v=0.999999;}
    double s = 2.*asin(v)/k;
//    double s_inv = 1./s;
 //   double sqrtvv = sqrt(1-v*v);
    double dz = sqrt(s*s*dzdl*dzdl/(1. - dzdl*dzdl));
    z = z0 + sign_dzdl*dz;
  }
  else
  {
    double s = 0.;
    double temp1 = k*D*0.5;temp1*=temp1;
    double temp2 = D*0.5;
    s += 2.*temp2;
    temp2*=temp1;
    s += temp2/3.;
    temp2*=temp1;
    s += (3./20.)*temp2;
    temp2*=temp1;
    s += (5./56.)*temp2;
   // double s_inv = 1./s;
    double dz = sqrt(s*s*dzdl*dzdl/(1. - dzdl*dzdl));
    z = z0 + sign_dzdl*dz;
  }
}

void sPHENIXTracker::updatelayer_sorted(vector<SimpleHit3D>& hits)
{
  for (unsigned int i=0;i<n_layers;i++)
    layer_sorted[i].clear();
  for (unsigned int i=0;i<hits.size();i++)
    layer_sorted[hits[i].layer].push_back(hits[i]);
}

void sPHENIXTracker::findSeededTracksByProjection(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
    vector<float> radii { 4.1, 7.5, 9.92, 25.0, 34.0, 43.0, 52.0, 61.0, 69.6, 78.2, 86.8, 96.5, 108.0 };
    
  findtracksiter += 1;
  
  if(seeds.size() == 0){return;}
  //unsigned int first_new_layer = seed_layer;//seeds[0].hits.size();
  
  /*for(unsigned int l=0;l<n_layers;++l)
  {
    layer_sorted[l].clear();
  }
  for(unsigned int i=0;i<hits.size();++i)
  {
    layer_sorted[hits[i].layer].push_back(hits[i]);
  }
  */
  float layer_phi_error[20];
  float layer_z_error[20];
  float layer_r_error[20];
  for(unsigned int l=0;l<n_layers;++l)
  {
    // additional uncertainty from multiple scattering
    float max_k = range.max_k;
    float Bfield_inv = 1./detector_B_field;
    float p_inv=0.;
    prev_p_inv = 3.33333333333333314e+02*max_k*Bfield_inv*sqrt(1. - range.max_dzdl*range.max_dzdl);
    p_inv=prev_p_inv;
    float total_scatter_2 = 0.;
    for(unsigned int i=seed_layer+1;i<=(l);++i)
    {
      float this_scatter = 0;//detector_scatter[i-1]*(radii[i]-radii[i-1])/radii[i];
      total_scatter_2 += this_scatter*this_scatter;
    }
    float angle = p_inv*sqrt(total_scatter_2)*1.0;
    float dsize = 0.5*(range.max_d-range.min_d);
    float angle_from_d = dsize/radii[l];
    float msval = 0.;
    if(angle_from_d > angle){msval=0.;}
    else{msval = (angle - angle_from_d);}
    
    if(layer_sorted.size() == 0){return;}
    float dr_big = 0.;
    float dz_big = 0.;
    for(unsigned int i=0;i<layer_sorted[l].size();++i)
    {
      float dx = layer_sorted[l][i].dx;
      float dy = layer_sorted[l][i].dy;
      float dr = dx*dx + dy*dy; 
      if(dr > dr_big){dr_big = dr;}
      float dz = layer_sorted[l][i].dz;
      if(dz > dz_big){dz_big = dz;}
    }
    layer_phi_error[l] = 2.*sqrt(dr_big)/radii[l] + msval;
    layer_z_error[l] = 2.*dz_big;
    layer_r_error[l] = 2.*sqrt(dr_big);
    cout<<"hits: "<<layer_sorted[l].size()<<" layer_z_error["<<l<<"]="<<layer_z_error[l]<<";  layer_phi_error["<<l<<"]="<< layer_phi_error[l]<<"; layer_r_error["<<l<<"] = "<<layer_r_error[l]<<endl;
  }

  for(unsigned int i=0;i<n_layers;++i)
  {
    angle_list[i].vec.clear();
  }
  
  //  for(unsigned int i=0;i<hits.size();i++)
  // {
  //  AngleIndexPair temppair(atan2(hits[i].y, hits[i].x), i);
  //  angle_list[hits[i].layer].addPair( temppair );
  // }
  for (unsigned int l=0;l<n_layers;l++)
    for (unsigned int i=0;i<layer_sorted[l].size();i++)
  {
    AngleIndexPair temppair(atan2(layer_sorted[l][i].y, layer_sorted[l][i].x), i);
    angle_list[l].addPair( temppair );
  }

  float phi_tol = 0.05;//0.025;
  float z_tol = 10.0;
  float r_tot = 10.0;
  vector<AngleIndexPair*> hit_candidates;
  
  vector<SimpleHit3D> cur_hits;
  vector<SimpleTrack3D> one_seed;

  vector<float> layerz;
  for (unsigned int l=0;l<n_layers;l++)
    if (layer_sorted[l].size()==0) layerz.push_back(0); else
      layerz.push_back(layer_sorted[l][0].z);
  //from here layer)sorted is not needed anymore
  


  for(unsigned int s=0;s<seeds.size();++s)
  {
    one_seed.clear();
    one_seed.push_back(seeds[s]);
    cur_hits.clear();
    for(unsigned int l=0;l<n_layers;++l)
    {
      hit_candidates.clear(); //if (layer_sorted[l].size()==0) continue;// MUSTN'T BE!
      layer_sorted[l].clear();
      float x,y,z;
      z =  layerz[l];//in the case of EC take the next z
      //z>120 - EC
      bool EC = abs(z)>120; 
      if (!EC)
	projectToLayer(seeds[s], l, x,y,z);
      else{cout<<"EC  "<<l<<"  ";
	projectToDisk(seeds[s],z,x,y);
      }
      float r = sqrt(x*x+y*y);//consider using r^2
      //cout << "Projection r3 "<<sqrt(x*x+y*y+z*z)<<" z: "<<z<<endl;
      float phi_error = phi_tol + layer_phi_error[l];
      double intersect_phi = atan2(y, x);
      angle_list[l].getRangeList(intersect_phi, phi_error, hit_candidates);
      float z_error = z_tol + layer_z_error[l];
      float r_error = r_tot + layer_r_error[l];
      for(unsigned int h=0;h<hit_candidates.size();++h)
      {
        SimpleHit3D* hit = &(hits[ hit_candidates[h]->index ]);
	float hitr = sqrt(hit->x*hit->x+hit->y*hit->y);
        if( (!EC && fabs(z - hit->z) < z_error)
	    || (EC && fabs(r-hitr) < r_error) ) {cur_hits.push_back(*hit); layer_sorted[l].push_back(*hit);}
      }
    }
    cout <<"  Segments with "<<cur_hits.size()<<" hits"<<endl;
    //    updatelayer_sorted(cur_hits);
    findSeededTracksbySegments(one_seed, cur_hits, tracks, range);
  }
}
