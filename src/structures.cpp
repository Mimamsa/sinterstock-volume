/** @file structures.cpp
 * @brief 固定建築外殼與下層量測表面過濾。 / Fixed building-shell and lower-measured-surface filtering.
 */
#include "cloud.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

size_t removeBuildingStructures(std::vector<Point> &points, const std::vector<std::array<double,4>> &planes, const StructureOptions &o) {
    if(!o.enabled) return 0;
    validateFootprint(o.footprint);
    if(planes.size()!=6 || !std::isfinite(o.wallMargin) || o.wallMargin<0 ||
       !std::isfinite(o.roofMargin) || o.roofMargin<0 || !std::isfinite(o.cellSize) || o.cellSize<.1 ||
       !std::isfinite(o.surfaceBand) || o.surfaceBand<=0) throw std::runtime_error("Invalid structure removal parameters/model");
    for(const auto &p:planes) {
        double norm=0;for(int i=0;i<4;++i) if(!std::isfinite(p[i])) throw std::runtime_error("Nonfinite structure plane");
        for(int i=0;i<3;++i) norm+=p[i]*p[i];
        if(std::abs(norm-1)>1e-3) throw std::runtime_error("Structure planes must have unit normals");
    }
    /** @brief 步驟 1：以四面牆的中心判定每面牆的室內側，屋頂則保留向下半空間。
     *         Step 1: use the real footprint when available, otherwise the four walls; retain downward roof half-spaces.
     * 平面順序沿用固定場域格式：兩面 X 牆、兩面 Y 牆、兩面屋頂。
     * Plane order follows the fixed site format: two X walls, two Y walls, then two roofs.
     */
    if(std::abs(planes[0][0])<.5 || std::abs(planes[1][0])<.5 ||
       std::abs(planes[2][1])<.5 || std::abs(planes[3][1])<.5 ||
       std::abs(planes[4][2])<.5 || std::abs(planes[5][2])<.5) throw std::runtime_error("Unexpected building plane orientation");
    const double cx=-(planes[0][3]/planes[0][0]+planes[1][3]/planes[1][0])/2;
    const double cy=-(planes[2][3]/planes[2][1]+planes[3][3]/planes[3][1])/2;
    if(!std::isfinite(cx) || !std::isfinite(cy)) throw std::runtime_error("Invalid wall orientation");
    auto inside=[&](const Point &p) {
        if(!insideFootprint(o.footprint,p.x,p.y,o.wallMargin)) return false;
        for(size_t i=o.footprint.empty()?0:4;i<6;++i) {
            const auto &a=planes[i]; const double d=a[0]*p.x+a[1]*p.y+a[2]*p.z+a[3];
            if(i<4) {const double side=a[0]*cx+a[1]*cy+a[3]; if((side>0?d:-d)<o.wallMargin) return false;}
            else if((a[2]>0?d:-d)>-o.roofMargin) return false;
        }
        return true;
    };
    /** @brief 步驟 2：用室內且低於屋頂排除帶的量測建立 XY 最低高度網格。
     *         Step 2: build XY minimum heights using interior measurements below the roof exclusion band.
     * 只用量測支點；以稀疏網格限制記憶體，不建立逐點鄰居清單。
     * Use measured supports only; a sparse grid avoids per-point neighbor lists.
     */
    using Key=std::pair<double,double>;
    auto key=[&](const Point &p) {return Key{std::floor(p.x/o.cellSize),std::floor(p.y/o.cellSize)};};
    std::map<Key,float> lower;
    for(const auto &p:points) if(inside(p)) {auto r=lower.emplace(key(p),p.z);if(!r.second) r.first->second=std::min(r.first->second,p.z);}
    /** @brief 步驟 3：原地保留下層表面帶，去掉同格更高的屋架、梁及垂直物件。
     *         Step 3: compact the lower surface band, removing higher frames, beams, and vertical objects in each cell.
     * 厚度容許斜坡在格內的高度差；過小會削掉陡坡，過大則保留更多鋼構。
     * Thickness allows within-cell slope variation; small values erode steep slopes, large values retain more structures.
     */
    const auto count=points.size();
    size_t write=0;
    for(size_t i=0;i<count;++i) {
        const auto &p=points[i];
        if(!inside(p) || p.z>lower.at(key(p))+o.surfaceBand) continue;
        points[write++]=p;
    }
    points.resize(write);
    return count-points.size();
}
