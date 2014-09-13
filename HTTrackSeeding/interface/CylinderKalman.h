#ifndef __CYLINDERKALMAN__
#define __CYLINDERKALMAN__

#include "HelixKalman.h"
#include <vector>

class HelixKalmanState;
//class SimpleHit3D;


class CylinderKalman : public HelixKalman
{
  public:
    CylinderKalman(std::vector<float>& detector_radii, std::vector<float>& detector_material, float B);
    virtual ~CylinderKalman();
    
  protected:
    void calculateProjections(PHENIXHough::SimpleHit3D& hit, HelixKalmanState& state, Eigen::Matrix<float,2,5>& H, Eigen::Matrix<float,2,1>& ha);
    void calculateMeasurements(PHENIXHough::SimpleHit3D& hit, Eigen::Matrix<float,2,1>& m, Eigen::Matrix<float,2,2>& G);
    bool calculateScatteringVariance(HelixKalmanState& state, float& var);
    void updateIntersection(HelixKalmanState& state, int layer);
    void subtractProjections(Eigen::Matrix<float,2,1>& m, Eigen::Matrix<float,2,1>& ha, Eigen::Matrix<float,2,1>& diff);
    
  private:
    void calculate_dxda(PHENIXHough::SimpleHit3D& hit, HelixKalmanState& state, Eigen::Matrix<float,3,5>& dxda, float& x, float& y, float& z);
    std::vector<float> det_rad;
    std::vector<float> det_scatter_variance;
    unsigned int nlayers;
    float signk_store;
};


#endif
