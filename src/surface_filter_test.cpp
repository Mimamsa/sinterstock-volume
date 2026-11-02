/** @file surface_filter_test.cpp
 * @brief 高低雜點、斜坡保留、補點不自我支持及關閉清理測試。 / Test high/low noise, slopes, non-self-supporting fill, and bypass.
 */
#include "cloud.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void require(bool x,const char *s){if(!x) throw std::runtime_error(s);}
int main() {try {
 SurfaceFilterOptions o;std::vector<Point> points;
 for(int y=0;y<15;++y) for(int x=0;x<15;++x) for(float d:{.1f,.2f,.3f}) {float xx=x*.5f+d,yy=y*.5f+.2f;points.push_back({xx,yy,1+.4f*xx+.2f*yy,42});}
 auto original=points;points.push_back({3.15f,3.15f,20,99});points.push_back({3.25f,3.25f,-10,99});
 std::vector<Point> generated{{3.2f,3.2f,30},{3.2f,3.2f,1+.4f*3.2f+.2f*3.2f}};std::vector<Point> rm,re;
 auto r=filterSurfaceOutliers(points,generated,o,rm,re);
 require(r.removedMeasured==2&&r.removedEstimated==1,"High/low measured or estimated outliers not removed");
 require(points.size()==original.size()&&generated.size()==1,"Slope points lost");for(size_t i=0;i<points.size();++i) require(points[i].x==original[i].x&&points[i].z==original[i].z&&points[i].intensity==42,"Surviving measurements changed");
 auto sparse=std::vector<Point>{{0,0,1},{.1f,.1f,10}};generated.clear();rm.clear();re.clear();r=filterSurfaceOutliers(sparse,generated,o,rm,re);require(sparse.size()==2,"Unsupported neighborhood erased genuine points");
 auto edge=original;for(auto &p:edge) p.z=p.x<3.5f?0:5;const auto edgeCount=edge.size();rm.clear();re.clear();
 r=filterSurfaceOutliers(edge,generated,o,rm,re);require(edge.size()==edgeCount,"Sharp real surface edge was smoothed away");
 auto line=std::vector<Point>{};for(int x=0;x<15;++x) line.push_back({x*.5f,.1f,1,42});line.push_back({3.1f,.1f,-10,99});rm.clear();re.clear();
 r=filterSurfaceOutliers(line,generated,o,rm,re);require(line.size()==16,"Degenerate support accepted as a 2D surface");
 o.enabled=false;auto unchanged=sparse;r=filterSurfaceOutliers(sparse,generated,o,rm,re);require(sparse.size()==unchanged.size()&&r.removedMeasured==0,"Disabled filter mutated points");
 std::cout<<"PASS: high/low residuals, measured-only support, slopes, coordinate preservation, sparse support and bypass\n";
 }catch(const std::exception &e){std::cerr<<e.what()<<"\n";return 1;}return 0;
}
