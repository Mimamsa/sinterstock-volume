/** @file surface_filter.cpp
 * @brief 只由實測建立局部料面，剔除雙向離群點。 / Measured-only robust local surfaces for two-sided outlier rejection.
 */
#include "cloud.h"
#include <Eigen/Dense>
#include <map>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace {
using Key=std::pair<int,int>;
struct Tile {std::vector<size_t> ids;Point representative{};};
struct Model {Eigen::Vector3d plane=Eigen::Vector3d::Zero();double threshold=0;bool valid=false;};
double median(std::vector<double> v) {auto m=v.begin()+v.size()/2;std::nth_element(v.begin(),m,v.end());return *m;}
}
SurfaceFilterReport filterSurfaceOutliers(std::vector<Point> &measured,std::vector<Point> &estimated,const SurfaceFilterOptions &o,std::vector<Point> &rejectedMeasured,std::vector<Point> &rejectedEstimated) {
 SurfaceFilterReport report;if(!o.enabled) return report;
 if(!std::isfinite(o.cell)||!std::isfinite(o.radius)||!std::isfinite(o.tolerance)||o.cell<.1||o.radius<o.cell||o.radius/o.cell>10||o.tolerance<=0) throw std::runtime_error("Invalid local surface filter settings");
 for(const auto *points:{&measured,&estimated}) for(const auto &p:*points)
    if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||std::abs(p.x)>1e6||std::abs(p.y)>1e6) throw std::runtime_error("Invalid surface filter coordinates");
 const auto key=[&](const Point &p){return Key{int(std::floor(p.x/o.cell)),int(std::floor(p.y/o.cell))};};
 std::map<Key,Tile> tiles;std::map<Key,Model> models;
 /** @brief 1：每格以實測高度中位點作等權支點，避免高密度點群及補點自行支持雜點。
  *         1: use one measured median-height point per cell; density and estimates never vote for themselves.
  */
 for(size_t i=0;i<measured.size();++i) tiles[key(measured[i])].ids.push_back(i);
 for(auto &entry:tiles) {auto &t=entry.second;auto mid=t.ids.begin()+t.ids.size()/2;std::nth_element(t.ids.begin(),mid,t.ids.end(),[&](auto a,auto b){return measured[a].z<measured[b].z;});t.representative=measured[*mid];models[entry.first];}
 for(const auto &p:estimated) models[key(p)];
 const int range=int(std::ceil(o.radius/o.cell));
 /** @brief 2：中心格不參與自身預測；鄰格不足、陡立或多層不一致處保留，不強行平滑邊緣。
  *         2: exclude the target cell; preserve insufficient, steep, or incoherent neighborhoods rather than smoothing edges blindly.
  */
 for(auto &entry:models) {
    const double x=(entry.first.first+.5)*o.cell,y=(entry.first.second+.5)*o.cell;std::vector<Eigen::Vector3d> samples;
    for(int dy=-range;dy<=range;++dy) for(int dx=-range;dx<=range;++dx) {
        if((dx==0&&dy==0)||std::hypot(dx*o.cell,dy*o.cell)>o.radius) continue;
        auto it=tiles.find({entry.first.first+dx,entry.first.second+dy});if(it==tiles.end()) continue;const auto &p=it->second.representative;samples.emplace_back(p.x-x,p.y-y,p.z);
    }
    if(samples.size()<6) continue;
    std::vector<double> zs;for(const auto &p:samples) zs.push_back(p.z());Eigen::Vector3d plane(0,0,median(zs));
    bool stable=true;
    for(int iteration=0;iteration<5;++iteration) {
        std::vector<double> residuals;for(const auto &p:samples) residuals.push_back(p.z()-plane.dot(Eigen::Vector3d(p.x(),p.y(),1)));
        const double bias=median(residuals);plane.z()+=bias;std::vector<double> deviations;for(auto r:residuals) deviations.push_back(std::abs(r-bias));const double cut=std::max(.15,2.5*1.4826*median(deviations));
        Eigen::Matrix3d a=Eigen::Matrix3d::Zero();Eigen::Vector3d b=Eigen::Vector3d::Zero();
        for(const auto &p:samples) {Eigen::Vector3d v(p.x(),p.y(),1);double r=std::abs(p.z()-plane.dot(v)),w=r<=cut?1:cut/r;a+=w*v*v.transpose();b+=w*v*p.z();}
        Eigen::LDLT<Eigen::Matrix3d> solver(a);if(!solver.isPositive()||solver.vectorD().minCoeff()<1e-6) {stable=false;break;}plane=solver.solve(b);
    }
    if(!stable||!plane.allFinite()||std::hypot(plane.x(),plane.y())>2) continue;
    std::vector<double> residuals;for(const auto &p:samples) residuals.push_back(p.z()-plane.dot(Eigen::Vector3d(p.x(),p.y(),1)));
    const double bias=median(residuals);plane.z()+=bias;std::vector<double> deviations;for(double r:residuals) deviations.push_back(std::abs(r-bias));const double sigma=1.4826*median(deviations);
    size_t coherent=0;for(auto r:residuals) if(std::abs(r-bias)<=std::max(.25,3*sigma)) ++coherent;
    if(sigma>.2||coherent*4<samples.size()*3) continue;
    entry.second.plane={plane.x(),plane.y(),plane.z()-plane.x()*x-plane.y()*y};entry.second.threshold=std::max(o.tolerance,4*sigma);entry.second.valid=true;
 }
 /** @brief 3：所有判斷使用同一份未刪除的實測支點，最後一次性刪除；不平均或移動任何點。
  *         3: all decisions use the same unmodified measured support, then erase once without averaging/moving points.
  */
 const auto filter=[&](std::vector<Point> &points,std::vector<Point> &rejected) {points.erase(std::remove_if(points.begin(),points.end(),[&](const Point &p){auto it=models.find(key(p));if(it==models.end()||!it->second.valid) return false;const auto &m=it->second;const bool remove=std::abs(p.z-m.plane.dot(Eigen::Vector3d(p.x,p.y,1)))>m.threshold;if(remove) rejected.push_back(p);return remove;}),points.end());};
 filter(measured,rejectedMeasured);filter(estimated,rejectedEstimated);report.removedMeasured=rejectedMeasured.size();report.removedEstimated=rejectedEstimated.size();return report;
}
