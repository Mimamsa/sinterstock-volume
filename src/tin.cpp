/**
 * @file tin.cpp
 * @brief 以分層、受限 Delaunay TIN 補償較大的遮擋區域。 / Fill larger occlusions using layered, bounded Delaunay TINs.
 */
#include "cloud.h"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {
/** @brief 有效量測格；XYZ 總和保留支點的實際位置。 / Measured cell; XYZ sums retain the actual support location. */
struct MeasuredCell {
    double x=0,y=0,z=0;
    float low=std::numeric_limits<float>::infinity(),high=-std::numeric_limits<float>::infinity();
    unsigned count=0;
};
/** @brief 接受的三角形及其線性高度梯度。 / Accepted triangle and its linear height gradient. */
struct Triangle { Point a,b,c; double determinant,gx,gy; };
/** @brief 由鄰近量測擬合的局部坡向。 / Local slope fitted from neighboring measurements. */
struct LocalSlope { double gx=0,gy=0; bool valid=false; };

/**
 * @brief 自適應 5x5 至 9x9 鄰域截尾擬合坡向；拒絕缺乏支撐、多層或過陡的支點。
 *        Fit trimmed local slopes in adaptive 5x5..9x9 support bins, rejecting unsupported, mixed, or overly steep sites.
 * @param sites 本高度帶量測支點。 / Measured supports in this height band.
 * @param o 同實際三角網格的限制。 / Limits shared with the actual TIN.
 * @return 每個支點的局部坡向與有效旗標。 / Local slope and validity for each site.
 */
std::vector<LocalSlope> localSlopes(const std::vector<Point> &sites,const FillOptions &o) {
    const double step=std::max(o.supportSpacing,std::sqrt((o.xMax-o.xMin)*(o.yMax-o.yMin)/60000.));
    const int columns=int(std::ceil((o.xMax-o.xMin)/step))+4;
    auto coords=[&](const Point &p) {return std::pair<int,int>{int((double(p.x)-o.xMin)/step),int((double(p.y)-o.yMin)/step)};};
    auto key=[&](int x,int y) {return qint64(y+2)*columns+x+2;};
    std::unordered_map<qint64,size_t> bins;
    for (size_t i=0;i<sites.size();++i) {auto [x,y]=coords(sites[i]); bins[key(x,y)]=i;}
    std::vector<LocalSlope> slopes(sites.size());
    for (size_t i=0;i<sites.size();++i) {
        auto [x,y]=coords(sites[i]); std::vector<size_t> neighbors;
        // 缺口邊緣支點不足時擴大至 9x9；仍只使用實際量測，不使用補點。
        // Expand up to 9x9 bins near sparsely sampled hole edges; never use synthetic points as evidence.
        for(int radius=2;radius<=4;++radius) {
            neighbors.clear();
            for(int dy=-radius;dy<=radius;++dy) for(int dx=-radius;dx<=radius;++dx) {
                if(x+dx<0 || y+dy<0 || x+dx>=columns-4) continue;
                auto it=bins.find(key(x+dx,y+dy));if(it!=bins.end()) neighbors.push_back(it->second);
            }
            if(neighbors.size()>=8) break;
        }
        if(neighbors.size()<5) continue;
        const size_t minimum=std::max(size_t(5),(neighbors.size()*7+9)/10);
        double gx=0,gy=0,meanX=0,meanY=0,meanZ=0;
        bool fitted=false;
        // 三次截尾擬合：少量殘留鋼構或掃描跳點不應讓整段缺口邊界的坡向失效。
        // Three trimmed fits prevent a few structure remnants or scan jumps from invalidating an entire hole edge.
        for(int iteration=0;iteration<3;++iteration) {
            meanX=meanY=meanZ=0;
            for(auto id:neighbors) {meanX+=sites[id].x;meanY+=sites[id].y;meanZ+=sites[id].z;}
            meanX/=neighbors.size();meanY/=neighbors.size();meanZ/=neighbors.size();
            double xx=0,xy=0,yy=0,xz=0,yz=0;
            for(auto id:neighbors) {
                const auto &p=sites[id];const double dx=p.x-meanX,dy=p.y-meanY,dz=p.z-meanZ;
                xx+=dx*dx;xy+=dx*dy;yy+=dy*dy;xz+=dx*dz;yz+=dy*dz;
            }
            const double determinant=xx*yy-xy*xy;
            if(determinant<1e-10) {fitted=false;break;}
            gx=(xz*yy-yz*xy)/determinant;gy=(yz*xx-xz*xy)/determinant;fitted=true;
            if(iteration==2) break;
            std::vector<std::pair<double,size_t>> residuals;
            for(auto id:neighbors) {const auto &p=sites[id];residuals.push_back({std::abs(p.z-(meanZ+gx*(p.x-meanX)+gy*(p.y-meanY))),id});}
            std::sort(residuals.begin(),residuals.end());
            const double threshold=std::max(.1,3*residuals[residuals.size()/2].first);
            size_t keep=residuals.size();while(keep>minimum && residuals[keep-1].first>threshold) --keep;
            if(keep==neighbors.size()) break;
            neighbors.clear();for(size_t n=0;n<keep;++n) neighbors.push_back(residuals[n].second);
        }
        if(!fitted || std::hypot(gx,gy)>o.maxSlope) continue;
        double squared=0,maxError=0;
        for(auto id:neighbors) {
            const auto &p=sites[id];const double error=std::abs(p.z-(meanZ+gx*(p.x-meanX)+gy*(p.y-meanY)));
            squared+=error*error;maxError=std::max(maxError,error);
        }
        const auto &center=sites[i];
        if(std::sqrt(squared/neighbors.size())>std::min(.25,o.maxError*.5) || maxError>o.maxError ||
           std::abs(center.z-(meanZ+gx*(center.x-meanX)+gy*(center.y-meanY)))>o.maxError) continue;
        slopes[i]={gx,gy,true};
    }
    return slopes;
}

/**
 * @brief 從實際支點建立 Delaunay 三角網格並排除過長、退化、陡峭三角形。
 *        Build a Delaunay TIN from measured supports and reject long, degenerate, or steep triangles.
 * @param supports 各支點只代表真實量測，不含生成點。 / Supports representing measurements only, never generated points.
 * @param o 跨接距離與坡度限制。 / Bridging distance and slope limits.
 * @param report 累加幾何拒絕統計。 / Accumulated rejection statistics.
 * @return 可用於三角形內部重心內插的集合。 / Triangles usable for interior barycentric interpolation.
 */
std::vector<Triangle> triangulate(const std::vector<Point> &supports,const FillOptions &o,FillReport &report) {
    if (supports.size()<3) return {};
    const auto &origin=supports.front();
    const auto farthest=std::max_element(supports.begin(),supports.end(),[&](Point a,Point b) {
        return std::hypot(a.x-origin.x,a.y-origin.y)<std::hypot(b.x-origin.x,b.y-origin.y);
    });
    bool noncollinear=false;
    for (const auto &p : supports) if (std::abs((double(farthest->x)-origin.x)*(double(p.y)-origin.y)-
        (double(farthest->y)-origin.y)*(double(p.x)-origin.x))>1e-8) {noncollinear=true; break;}
    if (!noncollinear) return {};
    const auto slopes=localSlopes(supports,o);
    std::vector<double> xy; xy.reserve(2*supports.size());
    for (auto p : supports) {xy.push_back(double(p.x)-origin.x); xy.push_back(double(p.y)-origin.y);}
    orgQhull::Qhull hull;
    std::ostringstream errors; hull.setErrorStream(&errors); hull.disableOutputStream();
    try {hull.runQhull("measured supports",2,int(supports.size()),xy.data(),"d Qbb Qt Qc Qz");}
    catch (const orgQhull::QhullError &e) {throw std::runtime_error(std::string("TIN triangulation failed: ")+e.what());}
    std::vector<Triangle> triangles; triangles.reserve(2*supports.size());
    for (auto facet=hull.beginFacet();facet!=hull.endFacet();facet=facet.next()) {
        if (facet.isUpperDelaunay()) continue;
        const auto vertices=facet.vertices().toStdVector();
        if (vertices.size()!=3) continue;
        std::array<int,3> id{}; bool valid=true;
        for (int i=0;i<3;++i) {id[i]=vertices[i].point().id(); if (id[i]<0 || size_t(id[i])>=supports.size()) valid=false;}
        if (!valid) continue;
        const auto a=supports[id[0]],b=supports[id[1]],c=supports[id[2]];
        // 凹多邊形不能當凸包處理；三邊都必須留在實際料場內。
        // A concave site is not its convex hull; all three triangle edges must remain inside.
        if(!footprintSegmentInside(o.footprint,a.x,a.y,b.x,b.y) || !footprintSegmentInside(o.footprint,a.x,a.y,c.x,c.y) ||
           !footprintSegmentInside(o.footprint,b.x,b.y,c.x,c.y)) {++report.rejectedFootprint;continue;}

        const double ux=double(b.x)-a.x,uy=double(b.y)-a.y,vx=double(c.x)-a.x,vy=double(c.y)-a.y;
        const double determinant=ux*vy-uy*vx;
        if (std::abs(determinant)<1e-10) continue;
        const double ab=std::hypot(ux,uy),ac=std::hypot(vx,vy),bc=std::hypot(double(b.x)-c.x,double(b.y)-c.y);
        const double circumradius=ab*ac*bc/(2*std::abs(determinant));
        // 外接圓限制同時拒絕橫跨大缺口的細長三角形。 / Circumradius also rejects skinny triangles spanning large voids.
        if (std::max({ab,ac,bc})>o.maxHole || circumradius>o.maxHole*.6) {++report.rejectedDistance; continue;}
        const double uz=double(b.z)-a.z,vz=double(c.z)-a.z;
        const double gx=(uz*vy-vz*uy)/determinant,gy=(ux*vz-vx*uz)/determinant;
        if (std::hypot(gx,gy)>o.maxSlope) {++report.rejectedSlope; continue;}
        /** @brief 不只檢查三角形本身，也要求至少兩個端點的量測坡向支持此接面。
         *         Require at least two measured endpoint slopes to support the triangle before bridging.
         * 端點法向與面法向差不得超過 45 度；端點外推至其他頂點的偏差亦受限。
         * Endpoint and face normals may differ by at most 45 degrees; cross-vertex continuation errors are also bounded.
         */
        bool consistent=true;
        int supportedVertices=0;
        const Point verticesPoints[]{a,b,c};
        const double tolerance=o.maxError+.05*std::max({ab,ac,bc});
        for (int i=0;i<3;++i) {
            const auto &s=slopes[id[i]];
            if (!s.valid) continue;
            ++supportedVertices;
            const double cosine=(1+gx*s.gx+gy*s.gy)/(std::sqrt(1+gx*gx+gy*gy)*std::sqrt(1+s.gx*s.gx+s.gy*s.gy));
            if (cosine<.7071067811865476) {consistent=false; break;}
            for (int j=0;j<3;++j) {
                const double z=verticesPoints[i].z+s.gx*(double(verticesPoints[j].x)-verticesPoints[i].x)+s.gy*(double(verticesPoints[j].y)-verticesPoints[i].y);
                if (std::abs(z-verticesPoints[j].z)>tolerance) consistent=false;
            }
        }
        if (!consistent || supportedVertices<2) {++report.rejectedConsistency; continue;}
        triangles.push_back({a,b,c,determinant,gx,gy});
    }
    return triangles;
}

/** @brief 只在三角形內部估計高度，禁止外插。 / Estimate height only inside a triangle, with no extrapolation.
 * @param t 已通過限制的三角形。 / Accepted triangle.
 * @param x,y 查詢世界 XY。 / Query world XY.
 * @param z 回傳估計高度。 / Receives estimated height.
 * @return 查詢點是否位於三角形內。 / Whether the query is inside the triangle.
 */
bool predict(const Triangle &t,double x,double y,double &z) {
    const double dx=x-t.a.x,dy=y-t.a.y;
    const double wb=(dx*(double(t.c.y)-t.a.y)-dy*(double(t.c.x)-t.a.x))/t.determinant;
    const double wc=((double(t.b.x)-t.a.x)*dy-(double(t.b.y)-t.a.y)*dx)/t.determinant;
    if (wb< -1e-8 || wc< -1e-8 || wb+wc>1+1e-8) return false;
    z=t.a.z+t.gx*dx+t.gy*dy; return std::isfinite(z);
}

/**
 * @brief 移除一個已知支點圓形區域，再以剩餘量測預測被移除的真值。
 *        Remove a known circular support patch and predict its withheld measured values.
 * @param sites 該高度帶量測支點。 / Measured supports in the band.
 * @param o 與實際補點相同的三角形限制。 / Same triangle limits used for filling.
 * @param report 累計覆蓋率與高度誤差；不是未知區域的準確度保證。
 *               Accumulated coverage and height errors, not an accuracy guarantee for unseen holes.
 */
void validatePatch(const std::vector<Point> &sites,const FillOptions &o,FillReport &report) {
    if (sites.size()<100 || report.validationPatches>=3) return;
    const double radius=std::min(5.,o.maxHole*.25);
    std::vector<Point> withheld,training;
    Point selectedCenter{};
    // 嘗試不同空間位置，避免選到稀疏孤點；整塊遮掉而非隨機留一點。
    // Try several locations to avoid sparse singletons; withhold contiguous patches rather than random single points.
    for (size_t attempt=1;attempt<=7;++attempt) {
        const auto center=sites[sites.size()*attempt/8];
        selectedCenter=center;
        withheld.clear(); training.clear(); training.reserve(sites.size());
        for (const auto &p : sites) {
            if (std::hypot(double(p.x)-center.x,double(p.y)-center.y)<radius) withheld.push_back(p);
            else training.push_back(p);
        }
        if (withheld.size()>=20 && withheld.size()<=1024 && training.size()>=30) break;
    }
    if (withheld.size()<20 || withheld.size()>1024 || training.size()<30) return;
    FillReport ignored;
    const auto triangles=triangulate(training,o,ignored);
    const double previousSquared=report.validationSquaredError;
    const size_t previousPredicted=report.validationPredicted;
    double patchMax=0;
    std::vector<unsigned char> used(withheld.size(),0);
    for (const auto &t : triangles) {
        const double xmin=std::min({t.a.x,t.b.x,t.c.x}),xmax=std::max({t.a.x,t.b.x,t.c.x});
        const double ymin=std::min({t.a.y,t.b.y,t.c.y}),ymax=std::max({t.a.y,t.b.y,t.c.y});
        for (size_t i=0;i<withheld.size();++i) {
            const auto &p=withheld[i]; if (used[i] || p.x<xmin || p.x>xmax || p.y<ymin || p.y>ymax) continue;
            double z;
            if (!predict(t,p.x,p.y,z)) continue;
            used[i]=1; ++report.validationPredicted;
            const double error=std::abs(z-p.z);
            report.validationSquaredError+=error*error; report.validationAbsoluteError+=error;
            report.validationMaxError=std::max(report.validationMaxError,error);
            patchMax=std::max(patchMax,error);
        }
    }
    ++report.validationPatches; report.validationTargets+=withheld.size();
    const size_t predicted=report.validationPredicted-previousPredicted;
    report.validation.push_back({selectedCenter,radius,withheld.size(),predicted,
        predicted ? std::sqrt((report.validationSquaredError-previousSquared)/predicted) : 0,patchMax});
}
}

std::vector<Point> reconstructStockpileSurface(const std::vector<Point> &points,const FillOptions &o,FillReport &report) {
    report={};
    validateFootprint(o.footprint);
    if (!o.enabled || points.empty()) return {};
    for (double v : {o.xMin,o.xMax,o.yMin,o.yMax,o.zMin,o.zMax,o.spacing,o.maxHole,o.maxError,o.maxSlope,o.supportSpacing,o.layerHeight,o.footprintMargin})
        if (!std::isfinite(v)) throw std::runtime_error("Non-finite TIN parameter");
    if (o.xMin>=o.xMax || o.yMin>=o.yMax || o.zMin>=o.zMax || o.spacing<.1 || o.maxHole<o.spacing ||
        o.footprintMargin<0 || o.maxError<=0 || o.maxSlope<=0 || o.supportSpacing<o.spacing || o.layerHeight<.5 || o.layerHeight<=2*o.maxError ||
        std::abs(o.zMin)>1e15 || std::abs(o.zMax)>1e15)
        throw std::runtime_error("Invalid TIN bounds, spacing, distance, slope, or height band");
    const double nxValue=std::ceil((o.xMax-o.xMin)/o.spacing),nyValue=std::ceil((o.yMax-o.yMin)/o.spacing);
    if (nxValue<2 || nyValue<2 || nxValue*nyValue>2000000 || points.size()>std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("TIN output grid exceeds 2,000,000 cells; increase spacing or reduce bounds");
    const int nx=int(nxValue),ny=int(nyValue);
    const double supportStep=std::max(o.supportSpacing,std::sqrt((o.xMax-o.xMin)*(o.yMax-o.yMin)/60000.));
    const int sx=int(std::ceil((o.xMax-o.xMin)/supportStep)),sy=int(std::ceil((o.yMax-o.yMin)/supportStep));
    const double bandStep=o.layerHeight*.5;
    const qint64 firstBand=qint64(std::floor(o.zMin/bandStep)),lastBand=qint64(std::ceil(o.zMax/bandStep));
    if (lastBand-firstBand>128) throw std::runtime_error("TIN requires over 128 height bands; increase layer height or limit Z");
    /** @brief 步驟 A：建立 Z 排序的索引，各量測只需進入約兩個重疊高度帶。
     *         Step A: sort indices by Z so each measurement enters approximately two overlapping height bands.
     */
    std::vector<uint32_t> order; order.reserve(points.size());
    for (size_t i=0;i<points.size();++i) {
        const auto &p=points[i];
        if (std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) &&
            p.x>=o.xMin && p.x<o.xMax && p.y>=o.yMin && p.y<o.yMax && p.z>=o.zMin && p.z<=o.zMax) order.push_back(uint32_t(i));
    }
    std::sort(order.begin(),order.end(),[&](uint32_t a,uint32_t b) {return points[a].z<points[b].z;});
    std::vector<qint64> validationBands;
    if (o.validate && !order.empty()) {
        // 按量測高度分布選 P20/P50/P80，避免只驗證最下方稀疏區。
        // Select bands near measured-height P20/P50/P80 instead of testing only the lowest sparse bands.
        for (double fraction : {.2,.5,.8}) validationBands.push_back(qint64(std::llround(points[order[size_t((order.size()-1)*fraction)]].z/bandStep)));
        std::sort(validationBands.begin(),validationBands.end());
        validationBands.erase(std::unique(validationBands.begin(),validationBands.end()),validationBands.end());
    }
    std::vector<MeasuredCell> cells(size_t(nx)*ny);
    std::vector<float> lastGenerated(cells.size(),std::numeric_limits<float>::quiet_NaN());
    std::vector<Point> output;
    for (qint64 band=firstBand;band<=lastBand;++band) {
        const double center=band*bandStep,low=center-o.layerHeight*.5,high=center+o.layerHeight*.5;
        const auto begin=std::lower_bound(order.begin(),order.end(),low,[&](uint32_t i,double z) {return points[i].z<z;});
        const auto end=std::lower_bound(begin,order.end(),high,[&](uint32_t i,double z) {return points[i].z<z;});
        if (end-begin<3) continue;
        std::fill(cells.begin(),cells.end(),MeasuredCell{});
        for (auto it=begin;it!=end;++it) {
            const auto &p=points[*it];
            int x=int(std::floor((double(p.x)-o.xMin)/o.spacing)),y=int(std::floor((double(p.y)-o.yMin)/o.spacing));
            if (x<0 || x>=nx || y<0 || y>=ny) continue;
            auto &c=cells[size_t(y)*nx+x]; c.x+=p.x; c.y+=p.y; c.z+=p.z; ++c.count;
            c.low=std::min(c.low,p.z); c.high=std::max(c.high,p.z);
        }
        /** @brief 步驟 B：在較粗 XY 網格保留最密集的可信量測格，限制三角網格記憶體。
         *         Step B: select the densest trustworthy measured cell per coarse XY bin to bound TIN memory.
         * 同格跨越多層高度的點不取平均當支點，降低屋頂／牆面與堆體混接。
         * Do not average large within-cell vertical spans into supports, reducing roof/wall/pile mixing.
         */
        std::vector<int> selected(size_t(sx)*sy,-1);
        for (int id=0;id<nx*ny;++id) {
            const auto &c=cells[id]; if (!c.count || double(c.high)-c.low>2*o.maxError) continue;
            const double x=c.x/c.count,y=c.y/c.count;
            int cx=int((x-o.xMin)/supportStep),cy=int((y-o.yMin)/supportStep);
            if (cx<0 || cx>=sx || cy<0 || cy>=sy) continue;
            auto &chosen=selected[size_t(cy)*sx+cx];
            if (chosen<0 || c.count>cells[chosen].count) chosen=id;
        }
        std::vector<Point> sites; sites.reserve(selected.size());
        for (int id : selected) if (id>=0) {
            const auto &c=cells[id]; sites.push_back({float(c.x/c.count),float(c.y/c.count),float(c.z/c.count)});
        }
        if (sites.size()<3) continue;
        ++report.bands; report.supports+=sites.size();
        if (o.validate && std::binary_search(validationBands.begin(),validationBands.end(),band)) validatePatch(sites,o,report);
        const auto triangles=triangulate(sites,o,report); report.triangles+=triangles.size();
        /** @brief 步驟 C：只在受限三角形內部、且該高度帶沒有量測的網格中心新增點。
         *         Step C: generate points only inside bounded triangles and cells lacking measurements in this band.
         * 每個帶只輸出中央半帶的高度，避免重疊帶重複補點；不要求空區完全封閉。
         * Emit only the central half-band in Z to avoid duplicates; empty regions need not be completely enclosed.
         */
        for (const auto &t : triangles) {
            int x0=std::max(0,int(std::ceil((std::min({t.a.x,t.b.x,t.c.x})-o.xMin)/o.spacing-.5)));
            int x1=std::min(nx-1,int(std::floor((std::max({t.a.x,t.b.x,t.c.x})-o.xMin)/o.spacing-.5)));
            int y0=std::max(0,int(std::ceil((std::min({t.a.y,t.b.y,t.c.y})-o.yMin)/o.spacing-.5)));
            int y1=std::min(ny-1,int(std::floor((std::max({t.a.y,t.b.y,t.c.y})-o.yMin)/o.spacing-.5)));
            const size_t previous=output.size();
            for (int y=y0;y<=y1;++y) for (int x=x0;x<=x1;++x) {
                const size_t id=size_t(y)*nx+x; if (cells[id].count) continue;
                const double xx=o.xMin+(x+.5)*o.spacing,yy=o.yMin+(y+.5)*o.spacing; double z;
                if (!insideFootprint(o.footprint,xx,yy,o.footprintMargin) || xx>=o.xMax || yy>=o.yMax || !predict(t,xx,yy,z) || z<o.zMin || z>o.zMax ||
                    z<center-bandStep*.5 || z>=center+bandStep*.5) continue;
                if (std::isfinite(lastGenerated[id]) && std::abs(z-lastGenerated[id])<std::max(.05,o.spacing)) continue;
                lastGenerated[id]=float(z); output.push_back({float(xx),float(yy),float(z)});
                if (output.size()>=2000000) {report.outputLimited=true; return output;}
            }
            if (output.size()>previous) ++report.patches;
        }
    }
    return output;
}
