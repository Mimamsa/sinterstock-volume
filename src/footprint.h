/** @file footprint.h
 * @brief 凹多邊形料場邊界，與校準平面獨立。 / Concave stockyard footprint, independent of calibration planes.
 */
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

using Footprint = std::vector<std::array<double,2>>;

/** @brief 驗證有限且非退化的簡單邊界；空值代表未指定。 / Validate a finite, nondegenerate simple boundary; empty means unspecified. */
inline void validateFootprint(const Footprint &polygon) {
    if(polygon.empty()) return;
    if(polygon.size()<3 || polygon.size()>256) throw std::runtime_error("Footprint requires 3..256 vertices");
    double area=0;
    auto cross=[](const auto &a,const auto &b,const auto &c) {return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);};
    for(size_t i=0;i<polygon.size();++i) {
        const auto &a=polygon[i],&b=polygon[(i+1)%polygon.size()];
        if(!std::isfinite(a[0]) || !std::isfinite(a[1]) || std::abs(a[0])>1e6 || std::abs(a[1])>1e6 || std::hypot(a[0]-b[0],a[1]-b[1])<1e-8)
            throw std::runtime_error("Invalid footprint vertex/edge");
        area+=a[0]*b[1]-a[1]*b[0];
        for(size_t j=i+2;j<polygon.size();++j) {
            if((j+1)%polygon.size()==i) continue;
            const auto &c=polygon[j],&d=polygon[(j+1)%polygon.size()];
            if(std::max(std::min(a[0],b[0]),std::min(c[0],d[0]))>std::min(std::max(a[0],b[0]),std::max(c[0],d[0])) ||
               std::max(std::min(a[1],b[1]),std::min(c[1],d[1]))>std::min(std::max(a[1],b[1]),std::max(c[1],d[1]))) continue;
            if(cross(a,b,c)*cross(a,b,d)<=0 && cross(c,d,a)*cross(c,d,b)<=0) throw std::runtime_error("Self-intersecting footprint");
        }
    }
    if(!std::isfinite(area) || std::abs(area)<1e-8) throw std::runtime_error("Degenerate footprint");
}

/** @brief 射線奇偶測試及真實邊界內縮；不以凸包或包圍矩形取代凹角。
 *         Ray-parity containment with actual-edge setback; never replace recesses with a convex hull or rectangle.
 * @param polygon 已驗證的料場邊界，空值不限制。 / Validated footprint; empty imposes no restriction.
 * @param x 公尺制 X。 / X in meters.
 * @param y 公尺制 Y。 / Y in meters.
 * @param margin 從真實外緣內縮的距離。 / Setback from the actual perimeter.
 */
inline bool insideFootprint(const Footprint &polygon,double x,double y,double margin=0) {
    if(polygon.empty()) return true;
    bool inside=false;double distance2=1e300;
    for(size_t i=0,j=polygon.size()-1;i<polygon.size();j=i++) {
        const auto &a=polygon[j],&b=polygon[i];
        const double dx=b[0]-a[0],dy=b[1]-a[1];
        const double t=std::clamp(((x-a[0])*dx+(y-a[1])*dy)/(dx*dx+dy*dy),0.,1.);
        const double ex=x-a[0]-t*dx,ey=y-a[1]-t*dy;distance2=std::min(distance2,ex*ex+ey*ey);
        if((a[1]>y)!=(b[1]>y) && x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0]) inside=!inside;
    }
    return (inside || distance2<1e-16) && distance2+1e-16>=margin*margin;
}

/** @brief 線段不得穿越凹多邊形外部，避免跨接料場缺角。
 *         Reject segments crossing outside a concave polygon, preventing bridges across site recesses.
 */
inline bool footprintSegmentInside(const Footprint &polygon,double ax,double ay,double bx,double by) {
    if(polygon.empty()) return true;
    if(!insideFootprint(polygon,ax,ay) || !insideFootprint(polygon,bx,by) || !insideFootprint(polygon,(ax+bx)*.5,(ay+by)*.5)) return false;
    auto side=[](double x,double y,double xx,double yy,double px,double py) {return (xx-x)*(py-y)-(yy-y)*(px-x);};
    for(size_t i=0;i<polygon.size();++i) {
        const auto &c=polygon[i],&d=polygon[(i+1)%polygon.size()];
        const double u=side(ax,ay,bx,by,c[0],c[1]),v=side(ax,ay,bx,by,d[0],d[1]);
        const double s=side(c[0],c[1],d[0],d[1],ax,ay),t=side(c[0],c[1],d[0],d[1],bx,by);
        if(u*v<-1e-16 && s*t<-1e-16) return false;
    }
    return true;
}
