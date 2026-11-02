/** @file volume.cpp
 * @brief 量測支撐的近水平地面 RANSAC 與固定 10 cm 體積積分。 / Measured near-horizontal ground RANSAC and fixed 10 cm volume integration.
 */
#include "cloud.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <stdexcept>

namespace {
struct Cell { double minimum=std::numeric_limits<double>::infinity();uint32_t count=0;float height=0;double netVolume=0,netArea=0;int label=0;bool valid=false,estimated=false; };

/** @brief 以等面積粗格低分位量測找地面，避免密集料堆點壓過稀疏露出地面。
 *         Fit ground from equal-area low-quantile measured samples, avoiding density bias toward piles.
 */
void estimateGround(const std::vector<Point> &points,const Footprint &footprint,double margin,const VolumeOptions &o,VolumeReport &r) {
    std::map<std::pair<int,int>,std::vector<Point>> tiles;
    for(const auto &p:points) if(std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)&&insideFootprint(footprint,p.x,p.y,margin)) {
        if(std::abs(p.x)>1e6 || std::abs(p.y)>1e6) throw std::runtime_error("Ground coordinates exceed supported range");
        tiles[{int(std::floor(p.x)),int(std::floor(p.y))}].push_back(p);
    }
    std::vector<Eigen::Vector3d> samples;
    for(auto &tile:tiles) {
        auto &z=tile.second;if(z.size()<3) continue;
        const size_t k=(z.size()-1)/10;std::nth_element(z.begin(),z.begin()+k,z.end(),[](const Point &a,const Point &b) {return a.z<b.z;});
        samples.emplace_back(z[k].x,z[k].y,z[k].z);
    }
    r.groundSamples=samples.size();
    if(samples.size()<30) throw std::runtime_error("Insufficient measured floor support; provide a manual ground plane");
    std::vector<double> heights;for(const auto &p:samples) heights.push_back(p.z());
    const size_t lowCount=(heights.size()*45)/100;std::nth_element(heights.begin(),heights.begin()+lowCount,heights.end());
    const double cutoff=heights[lowCount];
    std::vector<size_t> candidates;for(size_t i=0;i<samples.size();++i) if(samples[i].z()<=cutoff) candidates.push_back(i);
    std::mt19937 random(73129);std::uniform_int_distribution<size_t> pick(0,candidates.size()-1);
    size_t best=0;double bestScore=-1e30;std::array<double,3> plane{};
    auto votes=[&](const auto &a) {
        size_t n=0,below=0;
        for(const auto &p:samples) {const double residual=p.z()-a[0]*p.x()-a[1]*p.y()-a[2];
            if(std::abs(residual)<=o.groundTolerance) ++n;else if(residual< -o.groundTolerance) ++below;}
        // 地面不應浮在大量實際量測上方；懲罰地下量測，避免誤選平坦料面或斜坡。
        // Penalize measurements below a candidate floor, avoiding flat pile tops or ramps above exposed ground.
        return std::pair<size_t,double>{n,double(n)-1.5*below};
    };
    // 含水平候選，避免三個樣本幾乎共線；只接受與 Z 朝上相符、坡度小於 5 度的地面。
    // Include horizontal candidates and reject near-collinear triples; ground slope must be within 5 degrees of XY.
    for(int iteration=0;iteration<600;++iteration) {
        const auto &a=samples[candidates[pick(random)]];
        std::array<double,3> trial{0,0,a.z()};
        if(iteration%4) {
            const auto &b=samples[candidates[pick(random)]],&c=samples[candidates[pick(random)]];
            const double ux=b.x()-a.x(),uy=b.y()-a.y(),vx=c.x()-a.x(),vy=c.y()-a.y(),det=ux*vy-uy*vx;
            if(std::abs(det)<1.) continue;
            trial[0]=((b.z()-a.z())*vy-(c.z()-a.z())*uy)/det;
            trial[1]=(ux*(c.z()-a.z())-vx*(b.z()-a.z()))/det;
            trial[2]=a.z()-trial[0]*a.x()-trial[1]*a.y();
        }
        if(std::hypot(trial[0],trial[1])>std::tan(5.*3.141592653589793/180.)) continue;
        const auto score=votes(trial);
        if(score.second>bestScore) {best=score.first;bestScore=score.second;plane=trial;}
    }
    if(best<std::max(size_t(30),samples.size()/20)) throw std::runtime_error("No supported near-horizontal ground; provide a manual plane");
    for(int iteration=0;iteration<3;++iteration) {
        Eigen::Matrix3d normal=Eigen::Matrix3d::Zero();Eigen::Vector3d rhs=Eigen::Vector3d::Zero();
        for(const auto &p:samples) if(std::abs(p.z()-plane[0]*p.x()-plane[1]*p.y()-plane[2])<=o.groundTolerance) {
            const Eigen::Vector3d v(p.x(),p.y(),1);normal+=v*v.transpose();rhs+=v*p.z();
        }
        if(normal.fullPivLu().rank()<3) throw std::runtime_error("Ground support is collinear; provide a manual plane");
        const auto solution=normal.ldlt().solve(rhs).eval();
        if(!solution.allFinite() || std::hypot(solution.x(),solution.y())>std::tan(5.*3.141592653589793/180.))
            throw std::runtime_error("Unstable fitted ground; provide a manual plane");
        plane={solution.x(),solution.y(),solution.z()};
    }
    double squared=0,x0=1e30,x1=-1e30,y0=1e30,y1=-1e30,allX0=1e30,allX1=-1e30,allY0=1e30,allY1=-1e30;
    for(const auto &p:samples) {
        allX0=std::min(allX0,p.x());allX1=std::max(allX1,p.x());allY0=std::min(allY0,p.y());allY1=std::max(allY1,p.y());
        const double residual=p.z()-plane[0]*p.x()-plane[1]*p.y()-plane[2];
        if(std::abs(residual)>o.groundTolerance) continue;
        ++r.groundInliers;squared+=residual*residual;x0=std::min(x0,p.x());x1=std::max(x1,p.x());y0=std::min(y0,p.y());y1=std::max(y1,p.y());
    }
    r.groundSpanX=x1-x0;r.groundSpanY=y1-y0;
    if(r.groundInliers<std::max(size_t(30),samples.size()/20) || r.groundSpanX<std::max(2.,(allX1-allX0)*.3) || r.groundSpanY<std::max(2.,(allY1-allY0)*.3))
        throw std::runtime_error("Ground support lacks spatial coverage; provide a manual ground plane");
    r.ground=plane;r.groundRMSE=std::sqrt(squared/r.groundInliers);
}
}

VolumeReport calculateStockpileVolumes(const std::vector<Point> &measured,const std::vector<Point> &estimated,
                                      const Footprint &footprint,double margin,const VolumeOptions &o) {
    VolumeReport r;r.enabled=o.enabled;if(!o.enabled) return r;
    try {
        validateFootprint(footprint);
        if(footprint.empty()) throw std::runtime_error("An explicit site footprint is required; rectangular bounds are not a site boundary");
        if(!std::isfinite(margin)||margin<0 || !std::isfinite(o.groundTolerance)||o.groundTolerance<=0 ||
           !std::isfinite(o.pileThreshold)||o.pileThreshold<0 || !std::isfinite(o.minimumArea)||o.minimumArea<0)
            throw std::runtime_error("Invalid volume settings");
        /** @brief 步驟 1：僅從清理後量測估計地面；手動模式則使用指定平面。
         *         Step 1: estimate ground using cleaned measurements only, or use the specified manual plane.
         */
        if(o.automaticGround) estimateGround(measured,footprint,margin,o,r);
        else {for(double a:o.ground) if(!std::isfinite(a)) throw std::runtime_error("Nonfinite ground plane");r.ground=o.ground;r.manualGround=true;}
        r.groundValid=true;
        double x0=1e30,y0=1e30,x1=-1e30,y1=-1e30;
        for(const auto &p:footprint) {x0=std::min(x0,p[0]);x1=std::max(x1,p[0]);y0=std::min(y0,p[1]);y1=std::max(y1,p[1]);}
        if(x0>x1 || y0>y1 || std::max({std::abs(x0),std::abs(x1),std::abs(y0),std::abs(y1)})>1e6) throw std::runtime_error("Invalid volume bounds");
        constexpr double step=VolumeReport::spacing,area=step*step;
        x0=std::floor(x0/step)*step;y0=std::floor(y0/step)*step;
        const double nxValue=std::floor((x1-x0)/step)+1,nyValue=std::floor((y1-y0)/step)+1;
        if(nxValue*nyValue>4000000) throw std::runtime_error("10 cm volume raster exceeds 4,000,000 cells; restrict the footprint (grid is never enlarged)");
        const int nx=int(nxValue),ny=int(nyValue);std::vector<Cell> cells(size_t(nx)*ny);
        auto center=[&](size_t id) {return std::array<double,2>{x0+(id%nx+.5)*step,y0+(id/nx+.5)*step};};
        for(size_t id=0;id<cells.size();++id) {const auto p=center(id);if(insideFootprint(footprint,p[0],p[1],margin)) {
            // 整個 10 cm 方格須在凹多邊形內，避免以中心判斷而把邊界外面積算入。
            // Require the entire square inside the concave polygon, never counting exterior area via a center-only test.
            const double l=p[0]-step*.5,b=p[1]-step*.5;
            const std::array<std::array<double,2>,4> corners{{{l,b},{l+step,b},{l+step,b+step},{l,b+step}}};
            bool inside=true;for(size_t i=0;i<4;++i) {const auto &a=corners[i],&z=corners[(i+1)%4];
                if(!insideFootprint(footprint,a[0],a[1],margin)||!footprintSegmentInside(footprint,a[0],a[1],z[0],z[1])) {inside=false;break;}}
            if(inside) {cells[id].valid=true;++r.footprintCells;}
        }}
        auto add=[&](const std::vector<Point> &points,bool synthetic) {
            for(const auto &p:points) {
                if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||!insideFootprint(footprint,p.x,p.y,margin)) continue;
                const int x=int(std::floor((double(p.x)-x0)/step)),y=int(std::floor((double(p.y)-y0)/step));
                if(x<0||y<0||x>=nx||y>=ny) continue;
                auto &c=cells[size_t(y)*nx+x];if(!c.valid || (synthetic && c.count && !c.estimated)) continue;
                const double height=p.z-r.groundZ(p.x,p.y);
                if(!std::isfinite(height)||std::abs(height)>std::numeric_limits<float>::max()) throw std::runtime_error("Ground-relative height overflow");
                c.minimum=std::min(c.minimum,height);++c.count;c.estimated=synthetic;
            }
        };
        /** @brief 步驟 2：固定 10 cm XY 柵格，量測優先，僅在無量測格使用補點的最低離地高度。
         *         Step 2: fixed 10 cm XY cells use minimum measured heights first, estimates only where measurements are absent.
         * 地面高度在每個點的 XY 位置計算，斜地面不以全場單一 Z 常數取代。
         * Evaluate ground at each point's XY; sloped ground is never replaced by one site-wide Z constant.
         */
        add(measured,false);add(estimated,true);
        for(auto &c:cells) if(c.valid) {
            if(c.count) {++r.coveredCells;c.height=float(c.minimum);if(c.height<=o.pileThreshold) ++r.groundCells;if(c.height< -o.groundTolerance) ++r.belowGroundCells;}
            else ++r.unknownCells;
        }
        /** @brief 步驟 2b：以確認結構的三維聯集與每個料柱求交，露出部分與梁下空隙不扣。
         *         Step 2b: intersect the confirmed structural union with each material prism, preserving exposed parts and beam clearance.
         * 無料面資料的格子不先補出假體積再扣除；部分覆蓋格依交集面積扣，圓柱／工字形不以外接矩形扣。
         * Unknown cells never invent volume; partial cells use intersection area, and circular/I sections are not bounding boxes.
         */
        const ColumnVolumeMask columnMask(o.columns);
        std::map<int,double> columnTotals;
        for(size_t id=0;id<cells.size();++id) {auto &c=cells[id];if(!c.valid||!c.count||c.height<=o.pileThreshold) continue;
            const auto p=center(id);const auto excluded=columnMask.intersect(p[0]-step*.5,p[1]-step*.5,step,c.height,r.ground);
            for(const auto &part:excluded.contributions) columnTotals[part.first]+=part.second;
            r.grossRasterVolume+=area*c.height;r.columnExcludedVolume+=excluded.volume;r.columnExcludedArea+=excluded.fullHeightArea;
            c.netVolume=std::max(0.,area*c.height-excluded.volume);c.netArea=std::max(0.,area-excluded.fullHeightArea);
            if(c.netVolume<1e-10) c.netVolume=0;
            if(c.netArea<1e-10) c.netArea=0;
        }
        for(const auto &part:columnTotals) r.columnVolumes.push_back(part);
        /** @brief 步驟 3：高於地面門檻的格子以 8 鄰接分成料堆；小區不計入體積。
         *         Step 3: segment cells above the ground threshold into 8-connected piles; exclude small components.
         */
        std::vector<size_t> queue;int label=0;
        for(size_t id=0;id<cells.size();++id) {
            auto &initial=cells[id];if(!initial.valid||!initial.count||initial.height<=o.pileThreshold||initial.netVolume<=1e-12||initial.label) continue;
            ++label;initial.label=label;queue.clear();queue.push_back(id);
            PileVolume pile;
            for(size_t q=0;q<queue.size();++q) {
                const size_t cur=queue[q];const int x=int(cur%nx),y=int(cur/nx);const auto &c=cells[cur];const auto p=center(cur);
                pile.volume+=c.netVolume;pile.area+=c.netArea;pile.maxHeight=std::max(pile.maxHeight,double(c.height));pile.centerX+=p[0];pile.centerY+=p[1];
                if(c.estimated) {++pile.estimatedCells;pile.estimatedVolume+=c.netVolume;} else {++pile.measuredCells;pile.measuredVolume+=c.netVolume;}
                for(int dy=-1;dy<=1;++dy) for(int dx=-1;dx<=1;++dx) {
                    if(x+dx<0||x+dx>=nx||y+dy<0||y+dy>=ny) continue;
                    const size_t next=size_t(y+dy)*nx+x+dx;auto &n=cells[next];
                    if(n.valid&&n.count&&n.height>o.pileThreshold&&n.netVolume>1e-12&&!n.label) {n.label=label;queue.push_back(next);}
                }
            }

            if(pile.area+1e-9<o.minimumArea) {r.discardedSmallCells+=queue.size();r.discardedSmallVolume+=pile.volume;for(auto cur:queue) cells[cur].label=-1;continue;}
            pile.id=r.piles.size()+1;pile.centerX/=queue.size();pile.centerY/=queue.size();pile.meanHeight=pile.volume/std::max(pile.area,1e-12);
            // 臨時 label 改成連續的報告 ID；已拜訪的格子保持非零。
            // Replace provisional labels with consecutive report IDs, preserving visited state.
            for(auto cur:queue) cells[cur].label=int(pile.id);
            r.totalVolume+=pile.volume;r.totalArea+=pile.area;r.piles.push_back(pile);
        }
        /** @brief 步驟 4：未知格獨立統計；只有被單一料堆完整圍住的未知區才歸入該堆缺測面積。
         *         Step 4: track unknowns separately; attribute an enclosed unknown region only when bordered by one pile.
         * 未知格不補零、不憑空加入體積；連到場外的缺口僅列入全場缺測統計。
         * Unknowns neither become zero-height data nor contribute invented volume; exterior-connected gaps remain site-wide unknowns.
         */
        std::vector<unsigned char> visited(cells.size(),0);
        const int dx4[]{-1,1,0,0},dy4[]{0,0,-1,1};
        for(size_t id=0;id<cells.size();++id) if(cells[id].valid&&!cells[id].count&&!visited[id]) {
            queue.clear();queue.push_back(id);visited[id]=1;int border=0;bool open=false,mixed=false;
            for(size_t q=0;q<queue.size();++q) {
                const size_t cur=queue[q];const int x=int(cur%nx),y=int(cur/nx);
                for(int k=0;k<4;++k) {
                    const int xx=x+dx4[k],yy=y+dy4[k];
                    if(xx<0||xx>=nx||yy<0||yy>=ny) {open=true;continue;}
                    const size_t next=size_t(yy)*nx+xx;const auto &c=cells[next];
                    if(!c.valid) {open=true;continue;}
                    if(!c.count) {if(!visited[next]) {visited[next]=1;queue.push_back(next);}continue;}
                    if(c.label<=0) mixed=true;else if(!border) border=c.label;else if(border!=c.label) mixed=true;
                }
            }
            if(!open&&!mixed&&border>0) r.piles[size_t(border-1)].enclosedUnknownCells+=queue.size();
        }
        /** @brief 步驟 5：從柵格標籤取得邊界線；顯示上限不影響面積、周長或體積計算。
         *         Step 5: extract raster-label boundary edges; display limits never change area, perimeter, or volume.
         */
        constexpr size_t maximumEdges=5000;
        size_t edgeIndex=0;
        auto eachEdge=[&](auto action) {
            for(size_t id=0;id<cells.size();++id) if(cells[id].label>0) {
                const int x=int(id%nx),y=int(id/nx);const double xx=x0+x*step,yy=y0+y*step;
                for(int k=0;k<4;++k) {
                    const int xn=x+dx4[k],yn=y+dy4[k];
                    if(xn>=0&&xn<nx&&yn>=0&&yn<ny&&cells[size_t(yn)*nx+xn].label==cells[id].label) continue;
                    action(cells[id].label,k,xx,yy);
                }
            }
        };
        eachEdge([&](int pile,int,double,double) {++edgeIndex;r.piles[size_t(pile-1)].boundaryLength+=step;});
        const size_t stride=std::max(size_t(1),(edgeIndex+maximumEdges-1)/maximumEdges);r.boundaryDisplayLimited=edgeIndex>maximumEdges;edgeIndex=0;
        eachEdge([&](int pile,int k,double x,double y) {
            if(edgeIndex++%stride) return;
            const double x1=x+(k==1?step:0),y1=y+(k==3?step:0);
            r.boundaries.push_back({float(x1),float(y1),float(x1+(k>=2?step:0)),float(y1+(k<2?step:0)),size_t(pile)});
        });
    } catch(const std::exception &e) {r.error=QString::fromUtf8(e.what());r.totalVolume=0;r.totalArea=0;r.piles.clear();r.boundaries.clear();}
    return r;
}
