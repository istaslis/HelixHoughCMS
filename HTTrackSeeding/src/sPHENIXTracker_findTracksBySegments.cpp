#include "MLoVetere/HTTrackSeeding/interface/vector_math_inline.h"
#include "MLoVetere/HTTrackSeeding/interface/sPHENIXTracker.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include "MLoVetere/HTTrackSeeding/interface/Eigen/LU"
#include "MLoVetere/HTTrackSeeding/interface/Eigen/Core"
#include "MLoVetere/HTTrackSeeding/interface/Eigen/Dense"
#include <sys/time.h>


using namespace std;
using namespace Eigen;

typedef unsigned int uint;

const float dummyerror = 100; //1,10,100 works the same
const unsigned int dummyindex = (unsigned) -1; //it shouldn't repeat in one segment
const bool EmptyLayerDummy = false; //add dummies if there are no hits in one layer
const bool AlwaysAddDummy = false; //dummies if there are no GOOD hit in layer
const bool PrintStuff = false; //cout-debugging, can be removed safely

struct AlignedStruct
{
  float x1_a[4] __attribute__((aligned(16)));
  float x2_a[4] __attribute__((aligned(16)));
  float x3_a[4] __attribute__((aligned(16)));
  float y1_a[4] __attribute__((aligned(16)));
  float y2_a[4] __attribute__((aligned(16)));
  float y3_a[4] __attribute__((aligned(16)));
  float z1_a[4] __attribute__((aligned(16)));
  float z2_a[4] __attribute__((aligned(16)));
  float z3_a[4] __attribute__((aligned(16)));

  float dx1_a[4] __attribute__((aligned(16)));
  float dx2_a[4] __attribute__((aligned(16)));
  float dx3_a[4] __attribute__((aligned(16)));
  float dy1_a[4] __attribute__((aligned(16)));
  float dy2_a[4] __attribute__((aligned(16)));
  float dy3_a[4] __attribute__((aligned(16)));
  float dz1_a[4] __attribute__((aligned(16)));
  float dz2_a[4] __attribute__((aligned(16)));
  float dz3_a[4] __attribute__((aligned(16)));
  
  float kappa_a[4] __attribute__((aligned(16)));
  float dkappa_a[4] __attribute__((aligned(16)));
  
  float ux_mid_a[4] __attribute__((aligned(16)));
  float uy_mid_a[4] __attribute__((aligned(16)));
  float ux_end_a[4] __attribute__((aligned(16)));
  float uy_end_a[4] __attribute__((aligned(16)));
  
  float dzdl_1_a[4] __attribute__((aligned(16)));
  float dzdl_2_a[4] __attribute__((aligned(16)));
  float ddzdl_1_a[4] __attribute__((aligned(16)));
  float ddzdl_2_a[4] __attribute__((aligned(16)));
  
  float cur_kappa_a[4] __attribute__((aligned(16)));
  float cur_dkappa_a[4] __attribute__((aligned(16)));
  float cur_ux_a[4] __attribute__((aligned(16)));
  float cur_uy_a[4] __attribute__((aligned(16)));
  float cur_chi2_a[4] __attribute__((aligned(16)));
  float chi2_a[4] __attribute__((aligned(16)));
  
  unsigned int hit1[4];
  unsigned int hit2[4];
  unsigned int hit3[4];

};



void sPHENIXTracker::calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a)
{
  static const __m128 two = {2., 2., 2., 2.};
  
  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);
  
  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);
  
  
  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);
  
  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);
  
  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);
  
  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;
  
  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);
  
  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);
  
  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);
  
  __m128 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec_sqrt_ps(dr3);
  
  __m128 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m128 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);
  
  __m128 ux12 = (x2 - x1)*D12_inv;
  __m128 uy12 = (y2 - y1)*D12_inv;
  __m128 ux23 = (x3 - x2)*D23_inv;
  __m128 uy23 = (y3 - y2)*D23_inv;
  __m128 ux13 = (x3 - x1)*D31_inv;
  __m128 uy13 = (y3 - y1)*D31_inv;
  
  __m128 cosalpha = ux12*ux13 + uy12*uy13;
  __m128 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m128 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m128 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);
  
  __m128 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m128 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha*sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);
  
  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one - sinalpha*sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);
  
  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2)*D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);
}

//works only with hit_counter == 4
void sPHENIXTracker::calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a, float sinang_cut, float cosang_diff_inv, float* cur_kappa_a, float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a, float* chi2_a)
{
  static const __m128 two = {2., 2., 2., 2.};
  
  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);
  
  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);
  
  
  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);
  
  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);
  
  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);
  
  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;
  
  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);
  
  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);
  
  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);
  
  __m128 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec_sqrt_ps(dr3);
  
  __m128 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m128 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);
  
  __m128 ux12 = (x2 - x1)*D12_inv;
  __m128 uy12 = (y2 - y1)*D12_inv;
  __m128 ux23 = (x3 - x2)*D23_inv;
  __m128 uy23 = (y3 - y2)*D23_inv;
  __m128 ux13 = (x3 - x1)*D31_inv;
  __m128 uy13 = (y3 - y1)*D31_inv;
  
  __m128 cosalpha = ux12*ux13 + uy12*uy13;
  __m128 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m128 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m128 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);
  
  __m128 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m128 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha*sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);
  
  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one - sinalpha*sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);
  
  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2)*D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);
  
  __m128 c_dk = _mm_load_ps(cur_dkappa_a);
  __m128 c_k = _mm_load_ps(cur_kappa_a);
  __m128 c_ux = _mm_load_ps(cur_ux_a);
  __m128 c_uy = _mm_load_ps(cur_uy_a);
  __m128 c_chi2 = _mm_load_ps(cur_chi2_a);
  __m128 sinang = _mm_load1_ps(&sinang_cut);
  __m128 cosdiff = _mm_load1_ps(&cosang_diff_inv);
  
  __m128 kdiff = c_k - k;
  __m128 n_dk = c_dk + dk + sinang*k;
  __m128 chi2_k = kdiff*kdiff/(n_dk*n_dk);
  __m128 cos_scatter = c_ux*ux_mid + c_uy*uy_mid;
  __m128 chi2_ang = (one-cos_scatter)*(one-cos_scatter)*cosdiff*cosdiff;
  tmp1 = dzdl_1*sinang;
  _vec_fabs_ps(tmp1);
  __m128 chi2_dzdl = (dzdl_1 - dzdl_2)/(ddzdl_1 + ddzdl_2 + tmp1);
  chi2_dzdl *= chi2_dzdl;
  chi2_dzdl *= one_o_2;
  
  __m128 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
  _mm_store_ps(chi2_a, n_chi2);
}

void calculateKappaTangentsStruct(AlignedStruct &a4)
{
  sPHENIXTracker::calculateKappaTangents(a4.x1_a, a4.y1_a, a4.z1_a, a4.x2_a, a4.y2_a, a4.z2_a, a4.x3_a, a4.y3_a, a4.z3_a, a4.dx1_a, a4.dy1_a, a4.dz1_a, a4.dx2_a, a4.dy2_a, a4.dz2_a, a4.dx3_a, a4.dy3_a, a4.dz3_a, a4.kappa_a, a4.dkappa_a, a4.ux_mid_a, a4.uy_mid_a, a4.ux_end_a, a4.uy_end_a, a4.dzdl_1_a, a4.dzdl_2_a, a4.ddzdl_1_a, a4.ddzdl_2_a);

}

void calculateKappaTangentsStructWithChi2(AlignedStruct &a4, float sinang_cut, float cosang_diff_inv)
{
  sPHENIXTracker::calculateKappaTangents(a4.x1_a, a4.y1_a, a4.z1_a, a4.x2_a, a4.y2_a, a4.z2_a, a4.x3_a, a4.y3_a, a4.z3_a, a4.dx1_a, a4.dy1_a, a4.dz1_a, a4.dx2_a, a4.dy2_a, a4.dz2_a, a4.dx3_a, a4.dy3_a, a4.dz3_a, a4.kappa_a, a4.dkappa_a, a4.ux_mid_a, a4.uy_mid_a, a4.ux_end_a, a4.uy_end_a, a4.dzdl_1_a, a4.dzdl_2_a, a4.ddzdl_1_a, a4.ddzdl_2_a, sinang_cut, cosang_diff_inv, a4.cur_kappa_a, a4.cur_dkappa_a, a4.cur_ux_a, a4.cur_uy_a, a4.cur_chi2_a, a4.chi2_a);

}

//shoud be changed to SimpleHit3D field
bool IsDummy(SimpleHit3D &hit)
{
  return hit.dx == dummyerror;
}


//make dummy hit as linear interpolation from other two hits.
//dummyhitindex = 1 (dummy-hitA-hitB)
//dummyhitindex = 2 (hitA-dummy-hitB)
//dummyhitindex = 3 (hitA-hitB-dummy)
void MakeDummyHit(int dummyhitindex, SimpleHit3D &hitA, SimpleHit3D &hitB, SimpleHit3D &dummy)
{
  if (dummyhitindex == 1)
  {
    dummy.x = 2 * hitA.x - hitB.x;
    dummy.y = 2 * hitA.y - hitB.y;
    dummy.z = 2 * hitA.z - hitB.z;
    dummy.dx = dummyerror;
    dummy.dy = dummyerror;
    dummy.dz = dummyerror;
    dummy.layer = hitA.layer - 1;
    dummy.index = dummyindex;
  }
  else if (dummyhitindex == 2)
  {
    dummy.x = (hitA.x + hitB.x) / 2.0;
    dummy.y = (hitA.y + hitB.y) / 2.0 + 0.01; //added - to check if straight line works, and it works...
    dummy.z = (hitA.z + hitB.z) / 2.0;
    dummy.dx = dummyerror;
    dummy.dy = dummyerror;
    dummy.dz = dummyerror;
    dummy.layer = hitA.layer + 1;
    dummy.index = dummyindex;  //(int)(hitA.index + hitC.index)/2;//dummyindex;
  }
  else if (dummyhitindex == 3)
  {
    dummy.x = 2 * hitB.x - hitA.x;
    dummy.y = 2 * hitB.y - hitA.y;
    dummy.z = 2 * hitB.z - hitA.z;
    dummy.dx = dummyerror;
    dummy.dy = dummyerror;
    dummy.dz = dummyerror;
    dummy.layer = hitB.layer + 1;
    dummy.index = dummyindex;
  }
 
 
}

void InitThreeHits(AlignedStruct &a4, unsigned int hit_counter, SimpleHit3D &hitA, SimpleHit3D &hitB, SimpleHit3D &hitC)
{
  
  a4.x1_a[hit_counter] = hitA.x;
  a4.y1_a[hit_counter] = hitA.y;
  a4.z1_a[hit_counter] = hitA.z;
  a4.dx1_a[hit_counter] = hitA.dx;
  a4.dy1_a[hit_counter] = hitA.dy;
  a4.dz1_a[hit_counter] = hitA.dz;

  a4.x2_a[hit_counter] = hitB.x;
  a4.y2_a[hit_counter] = hitB.y;
  a4.z2_a[hit_counter] = hitB.z;
  a4.dx2_a[hit_counter] = hitB.dx;
  a4.dy2_a[hit_counter] = hitB.dy;
  a4.dz2_a[hit_counter] = hitB.dz;

  a4.x3_a[hit_counter] = hitC.x;
  a4.y3_a[hit_counter] = hitC.y;
  a4.z3_a[hit_counter] = hitC.z;
  a4.dx3_a[hit_counter] = hitC.dx;
  a4.dy3_a[hit_counter] = hitC.dy;
  a4.dz3_a[hit_counter] = hitC.dz;

}


void sPHENIXTracker::findTracksBySegments(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
      if (PrintStuff) {
        cout << "-------findTracksBySegments---------"<<endl;
        cout << "Range: Phi_min "<<range.min_phi<<endl;
        cout << "Range: Phi_max "<<range.max_phi<<endl;
        cout << "Hits to fit :"<<hits.size()<<endl;
        cout << "-------f..................s---------"<<endl;
    }
  
  segments1.clear();
  segments2.clear();
  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  
  for(unsigned int l=0;l<n_layers;++l)
  {
    layer_sorted[l].clear();
  }
  for(unsigned int i=0;i<hits.size();++i)
  {
    layer_sorted[hits[i].layer].push_back(hits[i]);
  }
 uint layerwithouthit = -1;

  //check if all layers have hits. If one of them is missing - dummyhit will be there
  //two missing hits are forbidden.
  for (unsigned int l = 0; l < n_layers; ++l)
  {
    if (layer_sorted[l].size() == 0)
    {
      if (PrintStuff) cout << "LAYER WITHOUT HITS: " << l << endl;
      if (layerwithouthit != (uint)-1) //two missed layers are forbidden
        return;
      layerwithouthit = l;
      if (!EmptyLayerDummy) return;
    }
  }
  
  if (PrintStuff) cout <<"Start fitting!"<<endl;

  timeval t1,t2;
  double time1=0.;
  double time2=0.;
  
  gettimeofday(&t1, NULL);
  
  float cosang_cut = 0;//0.999;
  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1./cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut*cosang_cut);
  float easy_chi2_cut = 1000000;//10.0;
  
  float inv_layer[20];
  for(unsigned int l=3;l<n_layers;++l)
  {
    inv_layer[l] = 1./(((float)l) - 2.);
  }
  
  unsigned int hit_counter = 0;
  
  AlignedStruct a4;
  

  
  TrackSegment temp_segment;
  
  unsigned int sizei=layer_sorted[0].size();
  unsigned int sizej=layer_sorted[1].size();
  unsigned int sizek=layer_sorted[2].size();

  if (sizei==0) sizei = 1;
  if (sizej==0) sizej = 1;
  if (sizek==0) sizek = 1;
  
  // make segments out of first 3 layers
  for(unsigned int i=0;i<sizei;++i)
  {
    for(unsigned int j=0;j<sizej;++j)
    {
      for(unsigned int k=0;k<sizek;++k)
      {

        if (layerwithouthit <=2)
        {
	  SimpleHit3D dummy;
          uint a = 0 ,b = 1;

          if (layerwithouthit==0) {a=1;b=2;}
          if (layerwithouthit==1) b=2;
          if (layerwithouthit==2) b=1; //can be safely removed

          MakeDummyHit(layerwithouthit+1, layer_sorted[a][j],layer_sorted[b][k], dummy);
          layer_sorted[layerwithouthit].push_back(dummy);

          if (PrintStuff){
          cout << "three-hits : "<<layerwithouthit<<endl;
          cout << "{{"<<dummy.x<<","<<dummy.y<<","<<dummy.z<<"},";
          cout << "{"<<layer_sorted[a][j].x<<","<<layer_sorted[a][j].y<<","<<layer_sorted[a][j].z<<"},";
          cout << "{"<<layer_sorted[b][k].x<<","<<layer_sorted[b][k].y<<","<<layer_sorted[b][k].z<<"}}";
          }

          //add 3 hits to array in the order depending on dummy hit location
          if (layerwithouthit==0) {
          InitThreeHits(a4, hit_counter, dummy, layer_sorted[1][j],
                          layer_sorted[2][k]);
                      a4.hit1[hit_counter] = layer_sorted[0].size() - 1;
                      a4.hit2[hit_counter] = j;
                      a4.hit3[hit_counter] = k;
                      hit_counter++;
          }
          else if (layerwithouthit==1)
          {

            InitThreeHits(a4, hit_counter, layer_sorted[0][i], dummy,
                layer_sorted[2][k]);
            a4.hit1[hit_counter] = i;
            a4.hit2[hit_counter] = layer_sorted[1].size() - 1;
            a4.hit3[hit_counter] = k;
            hit_counter++;
          } else //layerwithouthit==2
          {
            InitThreeHits(a4, hit_counter, layer_sorted[0][i],
                           layer_sorted[1][j], dummy);
                       a4.hit1[hit_counter] = i;
                       a4.hit2[hit_counter] = j;
                       a4.hit3[hit_counter] = layer_sorted[2].size() - 1;
                       hit_counter++;

          }
        }
        else
        {
          InitThreeHits(a4, hit_counter, layer_sorted[0][i], layer_sorted[1][j],
              layer_sorted[2][k]);
          a4.hit1[hit_counter] = i;
          a4.hit2[hit_counter] = j;
          a4.hit3[hit_counter] = k;
          hit_counter++;

          //dummies if there are no GOOD hit in layer
          if (AlwaysAddDummy)
          {
            SimpleHit3D dummyhit;
            MakeDummyHit(1, layer_sorted[1][j], layer_sorted[2][k],dummyhit);
            layer_sorted[0].push_back(dummyhit);

            InitThreeHits(a4, hit_counter, dummyhit, layer_sorted[1][j],
                layer_sorted[2][k]);
            a4.hit1[hit_counter] = layer_sorted[0].size() - 1;
            a4.hit2[hit_counter] = j;
            a4.hit3[hit_counter] = k;
            hit_counter++;

            MakeDummyHit(2, layer_sorted[0][i], layer_sorted[2][k], dummyhit);
            layer_sorted[1].push_back(dummyhit);

            InitThreeHits(a4, hit_counter, layer_sorted[0][i], dummyhit,
                layer_sorted[2][k]);
            a4.hit1[hit_counter] = i;
            a4.hit2[hit_counter] = layer_sorted[1].size() - 1;
            a4.hit3[hit_counter] = k;
            hit_counter++;

            MakeDummyHit(3, layer_sorted[0][i], layer_sorted[1][j], dummyhit);
            layer_sorted[2].push_back(dummyhit);

            InitThreeHits(a4, hit_counter, layer_sorted[0][i],
                layer_sorted[1][j], dummyhit);
            a4.hit1[hit_counter] = i;
            a4.hit2[hit_counter] = j;
            a4.hit3[hit_counter] = layer_sorted[2].size() - 1;
            hit_counter++;
          }
        }

        
        if(hit_counter == 4)
        {
          calculateKappaTangentsStruct(a4);
          
          for(unsigned int h=0;h<hit_counter;++h)
          {
            temp_segment.chi2 = (a4.dzdl_1_a[h] - a4.dzdl_2_a[h])/(a4.ddzdl_1_a[h] + a4.ddzdl_2_a[h] + fabs(a4.dzdl_1_a[h]*sinang_cut));
            temp_segment.chi2 *= temp_segment.chi2;
            if(temp_segment.chi2 > 2000.0){ continue;}
            temp_segment.ux = a4.ux_end_a[h];
            temp_segment.uy = a4.uy_end_a[h];
            temp_segment.kappa = a4.kappa_a[h];
            temp_segment.dkappa = a4.dkappa_a[h];
            temp_segment.hits[0] = a4.hit1[h];
            temp_segment.hits[1] = a4.hit2[h];
            temp_segment.hits[2] = a4.hit3[h];
            next_seg->push_back(temp_segment);
          }
          
          hit_counter=0;
        }
      }
    }
  }
  if(hit_counter != 0)
  {
    if (PrintStuff) {cout << "hit counter wasn't multiple of 4!" << endl;
    cout <<a4.x1_a[0]
         <<" "<< a4.x2_a[0]
         <<" "<< a4.x3_a[0]
         <<" "<<   a4.y1_a[0]
         <<" "<<   a4.y2_a[0]
         <<" "<<   a4.y3_a[0]
         <<" "<<   a4.z1_a[0]
         <<" "<<   a4.z2_a[0]
         <<" "<<   a4.z3_a[0]<<endl;
}
    calculateKappaTangentsStruct(a4);
    
    for(unsigned int h=0;h<hit_counter;++h)
    {
      temp_segment.chi2 = (a4.dzdl_1_a[h] - a4.dzdl_2_a[h])/(a4.ddzdl_1_a[h] + a4.ddzdl_2_a[h] + fabs(a4.dzdl_1_a[h]*sinang_cut));
      temp_segment.chi2 *= temp_segment.chi2;
      if(temp_segment.chi2 > 2000.0){/*cout << "WTF? "<<temp_segment.chi2<<endl; */continue;}
      temp_segment.ux = a4.ux_end_a[h];
      temp_segment.uy = a4.uy_end_a[h];
      temp_segment.kappa = a4.kappa_a[h];
      temp_segment.dkappa = a4.dkappa_a[h];
      temp_segment.hits[0] = a4.hit1[h];
      temp_segment.hits[1] = a4.hit2[h];
      temp_segment.hits[2] = a4.hit3[h];
      next_seg->push_back(temp_segment);
    }
    
    hit_counter=0;
  }
  swap(cur_seg, next_seg);
  
    
  
 if (PrintStuff) {
     cout << "cur_seg->size() after 3 layers :"<<cur_seg->size()<<endl;
  for(unsigned int i=0,sizei=cur_seg->size();i<sizei;++i)
  {
    cout << "segments? ";
    for(unsigned int l=0;l<3;++l)
    {
      cout << layer_sorted[l][(*cur_seg)[i].hits[l]].dx << ":"<<(*cur_seg)[i].hits[l]<<" ";
    }
    cout << endl;
  }
  cout << endl;
 }
  
  
  
  
  // add hits to segments layer-by-layer, cutting out bad segments
  unsigned int whichseg[4];
  for(unsigned int l=3;l<n_layers;++l)
  {
    if(l == (n_layers-1)){easy_chi2_cut*=0.25;}
    next_seg->clear();
    for(unsigned int i=0,sizei=cur_seg->size();i<sizei;++i)
    {
      for(unsigned int j=0,sizej=layer_sorted[l].size();j<sizej;++j)
      {
       if (IsDummy(layer_sorted[l][j])) continue;

	  //take last 2 hits without dummy
		if (IsDummy(layer_sorted[l-2][(*cur_seg)[i].hits[l-2]]))
		{
		    if (PrintStuff) cout << "On layer "<<l-2<<"there is a dummy"<<endl;
			  InitThreeHits(a4, hit_counter, layer_sorted[l-3][(*cur_seg)[i].hits[l-3]],
			      layer_sorted[l-1][(*cur_seg)[i].hits[l-1]],
			      layer_sorted[l][j]);
		}
		else
		if (IsDummy(layer_sorted[l-1][(*cur_seg)[i].hits[l-1]]))
		  {
		    if (PrintStuff) cout << "On layer "<<l-1<<"there is a dummy"<<endl;
			  InitThreeHits(a4, hit_counter, layer_sorted[l-3][(*cur_seg)[i].hits[l-3]],
			      layer_sorted[l-2][(*cur_seg)[i].hits[l-2]],
			      layer_sorted[l][j]);
		  }
		else
		  {
		  if (PrintStuff) cout << "No dummies in prev layers, adding new good one"<<endl;
		    //No dummy in segment
		  InitThreeHits(a4, hit_counter, layer_sorted[l-2][(*cur_seg)[i].hits[l-2]],
			      layer_sorted[l-1][(*cur_seg)[i].hits[l-1]],
			      layer_sorted[l][j]);
		  }

      whichseg[hit_counter] = i;
      a4.hit1[hit_counter] = j;
      hit_counter += 1;
      
        if(hit_counter == 4)
        {
          calculateKappaTangentsStruct(a4);//, sinang_cut, cosang_diff_inv);
          
          for(unsigned int h=0;h<hit_counter;++h)
          {
	    if (PrintStuff) cout << "45chi2 : "<<(a4.chi2_a[h])*inv_layer[l]<<endl;
	    float kdiff = (*cur_seg)[whichseg[h]].kappa - a4.kappa_a[h];
            float dk = (*cur_seg)[whichseg[h]].dkappa + a4.dkappa_a[h];
            dk += sinang_cut*a4.kappa_a[h];
            float chi2_k = kdiff*kdiff/(dk*dk);
            float cos_scatter = (*cur_seg)[whichseg[h]].ux*a4.ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*a4.uy_mid_a[h];
            float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
            float chi2_dzdl = (a4.dzdl_1_a[h] - a4.dzdl_2_a[h])/(a4.ddzdl_1_a[h] + a4.ddzdl_2_a[h] + fabs(a4.dzdl_1_a[h]*sinang_cut));
            chi2_dzdl *= chi2_dzdl;
            if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl)/((float)l - 2.) < easy_chi2_cut )
            {
              temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
              temp_segment.ux = a4.ux_end_a[h];
              temp_segment.uy = a4.uy_end_a[h];
              temp_segment.kappa = a4.kappa_a[h];
              temp_segment.dkappa = a4.dkappa_a[h];
              for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
              temp_segment.hits[l] = a4.hit1[h];
              next_seg->push_back(temp_segment);
            }
          }
          hit_counter=0;
        }
      }
      
          
      //if there are no dummy hit in the 3-hit segment - add it to 4 and 5th

          if (layer_sorted[l].size()==0 || AlwaysAddDummy){
              SimpleHit3D dummyhit;

              //todo: change to TrackSegment.hasDummy()
              bool nodummybefore = true;
              for (unsigned int ind=0;ind<l;ind++) {
                if (IsDummy(layer_sorted[ind][(*cur_seg)[i].hits[ind]]))
                {
                  nodummybefore = false;
                }
              }

              if (nodummybefore)//(!(*cur_seg)[i].hasDummy())
              {

                if (PrintStuff) cout << "No dummies in prev layers, adding DUMMY "<<(*cur_seg)[i].hasDummy()<<endl;
                MakeDummyHit(3, layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]],
                    layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]], dummyhit);
                layer_sorted[l].push_back(dummyhit);

                InitThreeHits(a4, hit_counter,
                    layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]],
                    layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]],
                    layer_sorted[l].back());

                whichseg[hit_counter] = i;
                a4.hit1[hit_counter] = layer_sorted[l].size() - 1;

                hit_counter += 1;

              }
          }

          if(hit_counter == 4)
          {
            calculateKappaTangentsStruct(a4);

            for(unsigned int h=0;h<hit_counter;++h)
                 {
                   float kdiff = (*cur_seg)[whichseg[h]].kappa - a4.kappa_a[h];
                   float dk = (*cur_seg)[whichseg[h]].dkappa + a4.dkappa_a[h];
                   dk += sinang_cut*a4.kappa_a[h];
                   float chi2_k = kdiff*kdiff/(dk*dk);
                   float cos_scatter = (*cur_seg)[whichseg[h]].ux*a4.ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*a4.uy_mid_a[h];
                   float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
                   float chi2_dzdl = (a4.dzdl_1_a[h] - a4.dzdl_2_a[h])/(a4.ddzdl_1_a[h] + a4.ddzdl_2_a[h] + fabs(a4.dzdl_1_a[h]*sinang_cut));
                   chi2_dzdl *= chi2_dzdl;
                   if (PrintStuff) cout << "45chi2 old : "<<((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k)/((float)l - 2.) <<endl;
                   if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k)/((float)l - 2.) < easy_chi2_cut )
                   {
                     temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
                     temp_segment.ux = a4.ux_end_a[h];
                     temp_segment.uy = a4.uy_end_a[h];
                     temp_segment.kappa = a4.kappa_a[h];
                     temp_segment.dkappa = a4.dkappa_a[h];
                     for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
                     temp_segment.hits[l] = a4.hit1[h];
                     next_seg->push_back(temp_segment);
                   }
                 }
                      hit_counter=0;

          }


      
      
    }
    if(hit_counter != 0)
    {
      // ????????????????????????????????????????????????????????????????
//      calculateKappaTangentsStructWithChi2(a4, a4.kappa_a, a4.dkappa_a, a4.ux_mid_a, a4.uy_mid_a, a4.ux_end_a, a4.uy_end_a, a4.dzdl_1_a, a4.dzdl_2_a, a4.ddzdl_1_a, a4.ddzdl_2_a, sinang_cut, cosang_diff_inv, a4.cur_kappa_a, a4.cur_dkappa_a, a4.cur_ux_a, a4.cur_uy_a, a4.cur_chi2_a, a4.chi2_a);
//
//      for(unsigned int h=0;h<hit_counter;++h)
//      {
//	 if (PrintStuff) cout << "45chi2 : "<<(a4.chi2_a[h])*inv_layer[l]<<endl;
//        if( (a4.chi2_a[h])*inv_layer[l] < easy_chi2_cut )
//        {
//          temp_segment.chi2 = a4.chi2_a[h];
//          temp_segment.ux = a4.ux_end_a[h];
//          temp_segment.uy = a4.uy_end_a[h];
//          temp_segment.kappa = a4.kappa_a[h];
//          temp_segment.dkappa = a4.dkappa_a[h];
//          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
//          temp_segment.hits[l] = a4.hit1[h];
//          next_seg->push_back(temp_segment);
//        }
//      }
//      hit_counter=0;
      

     calculateKappaTangentsStruct(a4);

      for(unsigned int h=0;h<hit_counter;++h)
      {
        float kdiff = (*cur_seg)[whichseg[h]].kappa - a4.kappa_a[h];
        float dk = (*cur_seg)[whichseg[h]].dkappa + a4.dkappa_a[h];
        dk += sinang_cut*a4.kappa_a[h];
        float chi2_k = kdiff*kdiff/(dk*dk);
        float cos_scatter = (*cur_seg)[whichseg[h]].ux*a4.ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*a4.uy_mid_a[h];
        float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
        float chi2_dzdl = (a4.dzdl_1_a[h] - a4.dzdl_2_a[h])/(a4.ddzdl_1_a[h] + a4.ddzdl_2_a[h] + fabs(a4.dzdl_1_a[h]*sinang_cut));
        chi2_dzdl *= chi2_dzdl;
	if (PrintStuff) cout << "45chi2 old : "<<((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k)/((float)l - 2.) <<endl;
        if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k)/((float)l - 2.) < easy_chi2_cut )
        {
          temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
          temp_segment.ux = a4.ux_end_a[h];
          temp_segment.uy = a4.uy_end_a[h];
          temp_segment.kappa = a4.kappa_a[h];
          temp_segment.dkappa = a4.dkappa_a[h];
          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
          temp_segment.hits[l] = a4.hit1[h];
          next_seg->push_back(temp_segment);
        }
      }
      hit_counter=0;

    }
    swap(cur_seg, next_seg);
  }
  
  
  if (PrintStuff) {
	  cout << "cur_seg size() : " << cur_seg->size()<<endl;
	  for(unsigned int i=0,sizei=cur_seg->size();i<sizei;++i)
	  {
		  cout << "segments "<<n_layers<<" ? ";
		  for(unsigned int l=0;l<n_layers;++l)
		  {
			  cout << layer_sorted[l][(*cur_seg)[i].hits[l]].dx << ":"<<(*cur_seg)[i].hits[l]<<" ";
		  }
		  cout << endl;
	  }
	  cout << endl;
  }
  
  
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  CAtime += (time2 - time1);
  
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  SimpleTrack3D temp_track_3hits;
  temp_track_3hits.hits.assign(3, SimpleHit3D());
  for(unsigned int i=0,sizei=cur_seg->size();i<sizei;++i)
  {
    if (PrintStuff) cout << "temp_comp ";
    for(unsigned int l=0;l<n_layers;++l)
    {
      temp_comb[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].index;
      if (PrintStuff) cout << layer_sorted[l][(*cur_seg)[i].hits[l]].dx<< " ";
    }
    if (PrintStuff) cout <<endl;
    sort(temp_comb.begin(),temp_comb.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_comb);
    if(it != combos.end()){continue;}
    combos.insert(temp_comb);
    
    if (PrintStuff) cout << "temp_track to be fitted:"<<endl;
    for(unsigned int l=0;l<n_layers;++l)
    {
      temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
      if (PrintStuff) cout << temp_track.hits[l].dx<<" "<<temp_track.hits[l].dy<<" "<<temp_track.hits[l].dz<<endl;
      if (PrintStuff) cout << temp_track.hits[l].x<<" "<<temp_track.hits[l].y<<" "<<temp_track.hits[l].z<<endl;
    }
    if (PrintStuff) cout << endl;

    
    gettimeofday(&t1, NULL);
    
    if (PrintStuff) cout << "temp_track_3hits : ";
    uint realhits = 0;
    for (unsigned int j = 0; j < 4; ++j)
    {
      if (IsDummy(temp_track.hits[j])) continue;


      temp_track_3hits.hits[realhits] = temp_track.hits[j];

      if (PrintStuff) cout << temp_track_3hits.hits[realhits].dx << " ";


      realhits++;
      if (realhits==3) break;
    }
    if (PrintStuff) cout << endl;
    
    
    float chi2tot = fitTrack_3(temp_track_3hits);
    
        if (PrintStuff)
      cout << "CHI2 fit3: " << chi2tot << endl;


    if (PrintStuff)
      {
        cout << "fast fit state: "<<endl;
        cout << "       phi:"<<temp_track_3hits.phi<<endl;
        cout << "     kappa:"<<temp_track_3hits.kappa<<endl;
        cout << "        z0:"<<temp_track_3hits.z0<<endl;

      }

    
    HelixKalmanState state;
    state.C = Matrix<float,5,5>::Zero(5,5);
    state.C(0,0) = 1.1;
    state.C(1,1) = 1.1;
    state.C(2,2) = 0.1;
    state.C(3,3) = 1.1;
    state.C(4,4) = 1.1;
    state.phi = temp_track_3hits.phi;
    if(state.phi < 0.){state.phi += 2.*M_PI;}
    state.d = temp_track_3hits.d;
    state.kappa = temp_track_3hits.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track_3hits.z0;
    state.dzdl = temp_track_3hits.dzdl;
    
    if (PrintStuff) cout << "OK? Fitting!" << endl;
    
    bool goodtrack = true;
    
      uint firsthitindex = 0;
    if (IsDummy(temp_track.hits[0])) firsthitindex = 1;

    
    for(unsigned int h=firsthitindex;h<temp_track.hits.size();++h)
    {
      kalman->addHit(temp_track.hits[h], state);nfits+=1;
      
      if(h >= 4)
      {
        if(state.chi2 != state.chi2)
        {
          goodtrack=false;break;
        }
        if(state.chi2/(2.*((float)(h+1)) - 5.) > chi2_cut)
        {
          goodtrack=false;break;
        }
      }
    }
    
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    KALtime += (time2 - time1);
    
    //kappa cut!!! pt<1
    if (state.kappa>kappa_cut) goodtrack = false;

    if(state.chi2 != state.chi2)
    {
      goodtrack=false;
    }
    if(state.chi2/(2.*((float)(n_layers)) - 5.) > chi2_cut)
    {
      goodtrack=false;
    }

    if (PrintStuff) cout<<"goodtrack : "<<goodtrack << " CHI2 : "<<state.chi2<<" CHI2norm : "<<state.chi2/(2.*((float)(n_layers)) - 5.)<<" kappa_cut "<<kappa_cut<< endl;
      
    if(goodtrack==false){continue;}
    
    temp_track.phi = state.phi;
    if(temp_track.phi < 0.){temp_track.phi += 2.*M_PI;}
    if(temp_track.phi > 2.*M_PI){temp_track.phi -= 2.*M_PI;}
    temp_track.d = state.d;
    temp_track.kappa = state.kappa;
    temp_track.z0 = state.z0;
    temp_track.dzdl = state.dzdl;
    tracks.push_back(temp_track);
    track_states.push_back(state);
    if((remove_hits == true) && (state.chi2 < chi2_removal_cut))
    {
      for(unsigned int i=0;i<temp_track.hits.size();++i)
      {
        if (temp_track.hits[i].index!=dummyindex)
            (*hit_used)[temp_track.hits[i].index] = true;
      }
    }

  }
}

vector<SimpleHit3D> GetLast3NonDummy(SimpleTrack3D &track, uint startfrom)
{
  vector <SimpleHit3D> lastthree;

  uint hitcounter = 0;
  for (int i=startfrom;i>=0;i--) //track.hits.size()-1
    if (!IsDummy(track.hits[i]) && hitcounter<3)
    {
      lastthree.push_back(track.hits[i]);
      hitcounter++;
    }
  return lastthree;
}

void sPHENIXTracker::findSeededTracksbySegments(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
    bool debug = false;

  if(seeds.size() == 0){return;}
  if (debug) cout << " SEGMENTS "<<hits.size()<<" n_l "<<n_layers<<endl;
  if (PrintStuff) {
    cout << "-------findSeededTracksbySegments---------"<<endl;
    cout << "Range: Phi_min "<<range.min_phi<<endl;
    cout << "Range: Phi_max "<<range.max_phi<<endl;
    cout << "Hits to fit :"<<hits.size()<<endl;
    cout << "-------f..................s---------"<<endl;
}

  //   for(unsigned int l=0;l<n_layers;++l)
  //  {
  //    layer_sorted[l].clear();
  //  }
  //for(unsigned int i=0;i<hits.size();++i)
  //   {
   //   layer_sorted[hits[i].layer].push_back(hits[i]);
  //  }
  //  for(unsigned int l=0;l<n_layers;++l)
  //  {
  //    cout << "            layer_sorted["<<l<<"].size(): "<<layer_sorted[l].size()<<endl;
  //  }

    if (debug)
        for (unsigned int i=0;i<n_layers;i++)
            cout << " new layer_sorted["<<i<<"].size() = "<<layer_sorted[i].size()<<endl;



  timeval t1,t2;
  double time1=0.;
  double time2=0.;

  
  gettimeofday(&t1, NULL);
  
  unsigned int first_new_layer = seed_layer;//!!!seeds[0].hits.size();
  
  segments1.clear();
  segments2.clear();
  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;



  float cosang_cut = 0.99;
  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1./cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut*cosang_cut);
  float easy_chi2_cut = 10.;
  
  unsigned int hit_counter = 0;
  float x1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float x2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float x3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float dx1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dx2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dx3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float kappa_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dkappa_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float ux_mid_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float uy_mid_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ux_end_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float uy_end_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float dzdl_1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dzdl_2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ddzdl_1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ddzdl_2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  TrackSegment temp_segment;
  unsigned int whichhit[4];
  unsigned int whichseed[4];
  for(unsigned int seed=0,seedsize=seeds.size();seed<seedsize;++seed)
  {
     if (seed_used[seeds[seed].index]) continue;
     vector<SimpleHit3D> hits3 = GetLast3NonDummy(seeds[seed],first_new_layer - 1);

    x1_a[hit_counter] = hits3[2].x;
    y1_a[hit_counter] = hits3[2].y;
    z1_a[hit_counter] = hits3[2].z;
    x2_a[hit_counter] = hits3[1].x;
    y2_a[hit_counter] = hits3[1].y;
    z2_a[hit_counter] = hits3[1].z;
    x3_a[hit_counter] = hits3[0].x;
    y3_a[hit_counter] = hits3[0].y;
    z3_a[hit_counter] = hits3[0].z;
    
    dx1_a[hit_counter] = hits3[2].dx;
    dy1_a[hit_counter] = hits3[2].dy;
    dz1_a[hit_counter] = hits3[2].dz;
    dx2_a[hit_counter] = hits3[1].dx;
    dy2_a[hit_counter] = hits3[1].dy;
    dz2_a[hit_counter] = hits3[1].dz;
    dx3_a[hit_counter] = hits3[0].dx;
    dy3_a[hit_counter] = hits3[0].dy;
    dz3_a[hit_counter] = hits3[0].dz;
    
    whichseed[hit_counter] = seed;
    hit_counter += 1;
    if(hit_counter == 4)
    {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
      
      for(unsigned int h=0;h<hit_counter;++h)
      {
        temp_segment.chi2 = 0.;
        temp_segment.ux = ux_end_a[h];
        temp_segment.uy = uy_end_a[h];
        temp_segment.kappa = kappa_a[h];
        temp_segment.dkappa = dkappa_a[h];
        temp_segment.seed = whichseed[h];
        next_seg->push_back(temp_segment);
      }
      
      hit_counter = 0;
    }
  }
  if(hit_counter != 0)
  {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
    
    for(unsigned int h=0;h<hit_counter;++h)
    {
      temp_segment.chi2 = 0.;
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.seed = whichseed[h];
      next_seg->push_back(temp_segment);
    }
    
    hit_counter = 0;
  }
    swap(cur_seg, next_seg);
  unsigned int whichseg[4];
  for(unsigned int l=first_new_layer;l<n_layers;++l)
  {
    next_seg->clear();
    for(unsigned int j=0,sizej=cur_seg->size();j<sizej;++j)
    {
      for(unsigned int i=0,sizei=layer_sorted[l].size();i<sizei;++i)
      {
        if((l-2)<=(first_new_layer-1))
        {
          vector<SimpleHit3D> hits3 = GetLast3NonDummy(seeds[(*cur_seg)[j].seed],seeds[(*cur_seg)[j].seed].hits.size()-1);
          SimpleHit3D neededhit = hits3[first_new_layer - (l-2) -1];

          x1_a[hit_counter] = neededhit.x;
          y1_a[hit_counter] = neededhit.y;
          z1_a[hit_counter] = neededhit.z;
          dx1_a[hit_counter] = neededhit.dx;
          dy1_a[hit_counter] = neededhit.dy;
          dz1_a[hit_counter] = neededhit.dz;
        }
        else
        {
          SimpleHit3D neededhit = IsDummy(layer_sorted[l-2][(*cur_seg)[j].hits[l-2]]) ? layer_sorted[l-3][(*cur_seg)[j].hits[l-3]]:
              layer_sorted[l-2][(*cur_seg)[j].hits[l-2]];

          x1_a[hit_counter] = neededhit.x;
          y1_a[hit_counter] = neededhit.y;
          z1_a[hit_counter] = neededhit.z;
          dx1_a[hit_counter] = neededhit.dx;
          dy1_a[hit_counter] = neededhit.dy;
          dz1_a[hit_counter] = neededhit.dz;
        }
        if((l-1)<=(first_new_layer-1))
        {
          vector<SimpleHit3D> hits3 = GetLast3NonDummy(seeds[(*cur_seg)[j].seed],seeds[(*cur_seg)[j].seed].hits.size()-1);
          SimpleHit3D neededhit = hits3[first_new_layer - (l-1) -1];

          x2_a[hit_counter] = neededhit.x;
          y2_a[hit_counter] = neededhit.y;
          z2_a[hit_counter] = neededhit.z;
          dx2_a[hit_counter] = neededhit.dx;
          dy2_a[hit_counter] = neededhit.dy;
          dz2_a[hit_counter] = neededhit.dz;
        }
        else
        {
SimpleHit3D neededhit = IsDummy(layer_sorted[l-1][(*cur_seg)[j].hits[l-1]]) ? layer_sorted[l-2][(*cur_seg)[j].hits[l-2]]:
              layer_sorted[l-1][(*cur_seg)[j].hits[l-1]];

          x2_a[hit_counter] = neededhit.x;
          y2_a[hit_counter] = neededhit.y;
          z2_a[hit_counter] = neededhit.z;
          dx2_a[hit_counter] = neededhit.dx;
          dy2_a[hit_counter] = neededhit.dy;
          dz2_a[hit_counter] = neededhit.dz;
        }
        x3_a[hit_counter] = layer_sorted[l][i].x;
        y3_a[hit_counter] = layer_sorted[l][i].y;
        z3_a[hit_counter] = layer_sorted[l][i].z;
        dx3_a[hit_counter] = layer_sorted[l][i].dx;
        dy3_a[hit_counter] = layer_sorted[l][i].dy;
        dz3_a[hit_counter] = layer_sorted[l][i].dz;
        
        whichhit[hit_counter] = i;
        whichseg[hit_counter] = j;
        hit_counter += 1;
        if(hit_counter == 4)
        {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
          
          for(unsigned int h=0;h<hit_counter;++h)
          {
            float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
            float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
            dk += sinang_cut*kappa_a[h];
            float chi2_k = kdiff*kdiff/(dk*dk);
            float cos_scatter = (*cur_seg)[whichseg[h]].ux*ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*uy_mid_a[h];
            float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
            float chi2_dzdl = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
            chi2_dzdl *= chi2_dzdl;
            if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl)/((float)l - 2.) < easy_chi2_cut )
            {
              temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              temp_segment.dkappa = dkappa_a[h];
              for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
              temp_segment.hits[l] = whichhit[h];
              temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
              next_seg->push_back(temp_segment);
            }
          }
          hit_counter = 0;
        }
      }
    }
    if(hit_counter != 0)
    {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

      for(unsigned int h=0;h<hit_counter;++h)
      {
        float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
        float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
        dk += sinang_cut*kappa_a[h];
        float chi2_k = kdiff*kdiff/(dk*dk);
        float cos_scatter = (*cur_seg)[whichseg[h]].ux*ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*uy_mid_a[h];
        float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
        float chi2_dzdl = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
        chi2_dzdl *= chi2_dzdl;
        if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl)/((float)l - 2.) < easy_chi2_cut )
        {
          temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          temp_segment.dkappa = dkappa_a[h];
          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
          temp_segment.hits[l] = whichhit[h];
          temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
          next_seg->push_back(temp_segment);
        }
      }
      hit_counter = 0;
    }
    swap(cur_seg, next_seg);
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  CAtime += (time2 - time1);

  candidates.push_back(cur_seg->size());

  dzdl.push_back((range.max_dzdl + range.min_dzdl)/2);
  
  gettimeofday(&t1, NULL);
  
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  for(unsigned int i=0,sizei=cur_seg->size();i<sizei;++i)
  {
    if (seed_used[seeds[(*cur_seg)[i].seed].index]) continue;
    for(unsigned int l=0;l<first_new_layer;++l)
    {
      temp_track.hits[l] = seeds[(*cur_seg)[i].seed].hits[l];
    }
    for(unsigned int l=first_new_layer;l<n_layers;++l)
    {
      temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
    }

    HelixKalmanState state = seed_states[seeds[(*cur_seg)[i].seed].index];
    bool goodtrack = true;
    for(unsigned int l=first_new_layer;l<n_layers;++l)
    {
      kalman->addHit(temp_track.hits[l], state);nfits+=1;
      if(state.chi2 != state.chi2)
      {
        goodtrack=false;break;
      }
      if(state.chi2/(2.*((float)(l+1)) - 5.) > chi2_cut)
      {
        goodtrack=false;break;
      }
    }
    temp_track.phi = state.phi;
    if(temp_track.phi < 0.){temp_track.phi += 2.*M_PI;}
    if(temp_track.phi > 2.*M_PI){temp_track.phi -= 2.*M_PI;}
    temp_track.d = state.d;
    temp_track.kappa = state.kappa;
    temp_track.z0 = state.z0;
    temp_track.dzdl = state.dzdl;

    if(goodtrack==false){continue;}
    if((remove_hits == true) && (state.chi2 < chi2_removal_cut))
    {
      for(unsigned int l=first_new_layer;l<n_layers;++l)
      {
        (*hit_used)[temp_track.hits[l].index] = true;
        temp_track.hits[l].index = index_mapping[temp_track.hits[l].index];
      }
    }
    else
    {
      for(unsigned int l=first_new_layer;l<n_layers;++l)
      {
        temp_track.hits[l].index = index_mapping[temp_track.hits[l].index];
      }
    }
    //    tracks.push_back(temp_track);
    //track_states.push_back(state);
    seed_used[seeds[(*cur_seg)[i].seed].index] = true;
    seed_matched_chi2[seeds[(*cur_seg)[i].seed].index] = state.chi2;
    seed_tracknhit[seeds[(*cur_seg)[i].seed].index] = n_layers;
  }
  
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  KALtime += (time2 - time1);
}



