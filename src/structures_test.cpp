/** @file structures_test.cpp
 * @brief 合成料堆與建築外殼的分離驗證。 / Synthetic stockpile/building-shell separation checks.
 */
#include "cloud.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
int main() {
    try {
        // 10 x 20 m building with equal-pitch roof; retain sloping material and reject its overhead beam.
        // 10 x 20 公尺雙斜屋頂：保留斜料面，移除其上方鋼梁。
        std::vector<std::array<double,4>> planes{{1,0,0,0},{1,0,0,-10},{0,1,0,0},{0,1,0,-20}};
        const double n=std::sqrt(1.04);
        planes.push_back({-.2/n,0,1/n,0});planes.push_back({.2/n,0,1/n,-2/n});
        std::vector<Point> points;
        for(int x=10;x<90;++x) for(int y=10;y<190;++y) {
            const float px=x*.1f,py=y*.1f;
            points.push_back({px,py,-4+px*.1f,42}); // Sloping measured pile / 斜料面。
            points.push_back({px,py,-1,99}); // Lower steelwork with material below / 下方有料面的鋼構。
        }
        const auto expected=points.size()/2;
        points.push_back({0,5,-4,99});points.push_back({5,0,-4,99}); // Walls / 牆面。
        points.push_back({5,5,1,99});points.push_back({-1,5,-4,99}); // Roof and exterior / 屋頂及外部。
        const auto original=points.size();StructureOptions options;options.enabled=true;options.roofMargin=.5;
        const auto removed=removeBuildingStructures(points,planes,options);
        if(points.size()!=expected || removed!=original-expected) throw std::runtime_error("Unexpected retained surface count");
        for(const auto &p:points) if(p.intensity!=42 || std::abs(p.z-(-4+p.x*.1f))>1e-5) throw std::runtime_error("Material/intensity changed");
        // 兩端延伸料地及接縫必須保留，真正凹角則排除；不得改變校準平面。
        // Preserve both end wings and internal seams while excluding true recesses, without changing calibration planes.
        options.footprint={{0,0},{4,0},{4,-4},{10,-4},{10,26},{6,26},{6,24},{4,24},{4,20},{0,20}};
        std::vector<Point> wings{{8,-2,-4,1},{8,22,-4,2},{8,0,-4,3},{8,20,-4,4},{8,25,-4,5},
                                 {2,-2,-4,99},{2,22,-4,99},{5,25,-4,99},{11,10,-4,99},{8,22,1,99}};
        if(removeBuildingStructures(wings,planes,options)!=5 || wings.size()!=5) throw std::runtime_error("Wing footprint count failed");
        for(const auto &p:wings) if(p.intensity==99) throw std::runtime_error("Recess/roof survived");
        if(!insideFootprint(options.footprint,8,0,.4) || insideFootprint(options.footprint,2,-2)) throw std::runtime_error("Footprint seam or recess failed");
        if(footprintSegmentInside(options.footprint,2,19,5,23)) throw std::runtime_error("Segment crossed recess");
        bool badPolygon=false;try {validateFootprint({{0,0},{2,2},{0,2},{2,0}});} catch(const std::runtime_error &) {badPolygon=true;}
        if(!badPolygon) throw std::runtime_error("Self crossing footprint accepted");
        options.enabled=false; if(removeBuildingStructures(points,{},options)!=0 || points.size()!=expected) throw std::runtime_error("Disabled filter changed data");
        options.enabled=true;bool rejected=false;
        try {removeBuildingStructures(points,{},options);} catch(const std::runtime_error &) {rejected=true;}
        if(!rejected) throw std::runtime_error("Missing model accepted");
        std::cout<<"Structure separation passed\n";
    } catch(const std::exception &e) {std::cerr<<e.what()<<"\n";return 1;}
}
