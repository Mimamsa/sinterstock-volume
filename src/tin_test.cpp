/**
 * @file tin_test.cpp
 * @brief 大缺口、曲面、分層與遮擋驗證測試。 / Large-gap, curved-surface, layering, and holdout-validation tests.
 */
#include "cloud.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

/** @brief 驗證大型人工遮擋重建，不以新增點數作為唯一成功條件。
 *         Verify reconstruction of artificial large occlusions rather than relying on output counts alone.
 * @return 全部通過為 0。 / Zero if all checks pass.
 */
int main() {
    FillOptions o{true,0,26.2,0,26.2,0,12,.2,22,.15,1};
    o.method=FillMethod::LayeredTin; o.supportSpacing=.4; o.layerHeight=8; o.validate=true;
    std::vector<Point> measured;
    for (int y=0;y<131;++y) for (int x=0;x<131;++x) {
        if (x>=30 && x<100 && y>=30 && y<100) continue;
        const float xx=float((x+.5)*.2),yy=float((y+.5)*.2);
        measured.push_back({xx,yy,2+.1f*xx+.05f*yy,42});
    }
    FillReport report;
    const auto result=reconstructStockpileSurface(measured,o,report);
    if (result.size()<4500 || !report.patches || !report.validationPredicted) {
        std::cerr << "Large-hole coverage insufficient: " << result.size() << "\n"; return 1;
    }
    double squared=0;
    for (const auto &p : result) {
        const double error=p.z-(2+.1*p.x+.05*p.y); squared+=error*error;
        if (!std::isnan(p.intensity) || p.x<6 || p.x>20 || p.y<6 || p.y>20 || std::abs(error)>1e-4) return 1;
    }
    size_t holes=0;
    if (!fillStockpileHoles(measured,o,holes).empty()) {std::cerr << "Legacy baseline unexpectedly filled >4096 cells\n"; return 1;}
    if (std::sqrt(report.validationSquaredError/report.validationPredicted)>1e-4) return 1;
    std::cout << "Large 14 m square mask: " << result.size() << "/4900 points; RMSE " << std::sqrt(squared/result.size()) << " m\n";

    auto bounded=o; bounded.maxHole=2; bounded.validate=false;
    const auto limited=reconstructStockpileSurface(measured,bounded,report);
    for (const auto &p : limited) if (std::hypot(p.x-13.f,p.y-13.f)<3) return 1;
    if (!report.rejectedDistance) return 1;
    // 曲面不是用一張平面補過去；與已知真值比較誤差。 / Check a curved surface against known truth instead of filling one plane.
    std::vector<Point> curved;
    for (int y=0;y<41;++y) for (int x=0;x<41;++x) {
        if (x>=13 && x<=27 && y>=13 && y<=27) continue;
        const float xx=float((x+.5)*.2),yy=float((y+.5)*.2);
        curved.push_back({xx,yy,1+.02f*xx*xx+.01f*yy*yy,42});
    }
    auto curveOptions=o; curveOptions.xMax=curveOptions.yMax=8.2; curveOptions.maxHole=8;
    curveOptions.supportSpacing=.2; curveOptions.validate=false;
    const auto curve=reconstructStockpileSurface(curved,curveOptions,report);
    if (curve.size()<200) return 1;
    squared=0;
    for (const auto &p : curve) {const double error=p.z-(1+.02*p.x*p.x+.01*p.y*p.y); squared+=error*error;}
    const double curvedRMSE=std::sqrt(squared/curve.size());
    if (curvedRMSE>.2) {std::cerr << "Curved RMSE " << curvedRMSE << "\n"; return 1;}
    std::cout << "Curved mask: " << curve.size() << "/225 points; RMSE " << curvedRMSE << " m\n";

    std::vector<Point> layers;
    for (int y=0;y<9;++y) for (int x=0;x<9;++x) {
        if (x>=3 && x<=5 && y>=3 && y<=5) continue;
        const float xx=float((x+.5)*.2),yy=float((y+.5)*.2),z=1+.2f*xx+.3f*yy;
        layers.push_back({xx,yy,z,42}); layers.push_back({xx,yy,z+12,42});
    }
    auto layer=o; layer.xMax=layer.yMax=1.8; layer.zMax=20; layer.maxHole=1;
    layer.supportSpacing=.2; layer.layerHeight=4; layer.validate=false;
    const auto multi=reconstructStockpileSurface(layers,layer,report);
    if (multi.size()!=18) {std::cerr << "Layered output " << multi.size() << "\n"; return 1;}
    for (const auto &p : multi) {
        const double base=1+.2*p.x+.3*p.y;
        if (std::min(std::abs(p.z-base),std::abs(p.z-base-12))>1e-4) return 1;
    }
    // 兩個高度不同但坡向皆水平的平台，不應因長三角形坡度不大而被接成斜坡。
    // Two flat platforms at different heights must not be joined into a ramp merely because a long triangle has a moderate slope.
    std::vector<Point> platforms;
    for (int y=0;y<21;++y) for (int x=0;x<61;++x) {
        if (x>20 && x<40) continue;
        platforms.push_back({float((x+.5)*.2),float((y+.5)*.2),x<=20 ? 1.f : 5.f});
    }
    auto flat=o; flat.xMax=12.2; flat.yMax=4.2; flat.zMax=12; flat.layerHeight=20;
    flat.maxHole=10; flat.maxSlope=3; flat.maxError=.2; flat.validate=false;
    const auto separated=reconstructStockpileSurface(platforms,flat,report);
    for (const auto &p : separated) if (p.x>5 && p.x<7) return 1;
    if (!report.rejectedConsistency) return 1;
    // 凹角不是待補洞區域；內部料堆缺口仍可補，三角網格不可跨接邊界凹角。
    // Site recesses are not holes: retain interior filling while rejecting triangles crossing the concave perimeter.
    auto concave=o;concave.validate=false;
    concave.footprint={{0,0},{26.2,0},{26.2,26.2},{18,26.2},{18,12},{8,12},{8,26.2},{0,26.2}};
    auto notched=measured;
    notched.erase(std::remove_if(notched.begin(),notched.end(),[&](const Point &p) {return !insideFootprint(concave.footprint,p.x,p.y);}),notched.end());
    const auto notchedFill=reconstructStockpileSurface(notched,concave,report);
    if(notchedFill.empty() || !report.rejectedFootprint) {std::cerr<<"Concave footprint test produced no interior fill/rejections\n";return 1;}
    for(const auto &p:notchedFill) if(!insideFootprint(concave.footprint,p.x,p.y) || std::abs(p.z-(2+.1*p.x+.05*p.y))>1e-4) return 1;
    std::vector<Point> collinear{{0,0,1},{.2f,0,1},{.4f,0,1}};
    if (!reconstructStockpileSurface(collinear,layer,report).empty()) return 1;
    auto invalid=o; invalid.layerHeight=.0001;
    bool rejected=false;
    try {reconstructStockpileSurface(measured,invalid,report);} catch (const std::runtime_error &) {rejected=true;}
    if (!rejected) return 1;
    std::cout << "PASS: distance limits, multiple height surfaces, no extrapolation, degenerate inputs, holdout error\n";
    return 0;
}
