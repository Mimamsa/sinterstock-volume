/** @file volume_test.cpp
 * @brief 解析體積、斜地面、獨立料堆、缺測與地面失敗驗證。 / Analytic volumes, sloped ground, separate piles, unknowns, and ground rejection.
 */
#include "cloud.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void require(bool value,const char *message) {if(!value) throw std::runtime_error(message);}
int main() {
    try {
        VolumeOptions o;o.enabled=true;o.automaticGround=false;o.ground={.02,-.01,-5};o.minimumArea=.01;
        Footprint footprint{{0,0},{4,0},{4,2},{0,2}};
        std::vector<Point> measured,estimated;
        for(int y=0;y<20;++y) for(int x=0;x<40;++x) {
            const float xx=(x+.5f)*.1f,yy=(y+.5f)*.1f;const double ground=.02*xx-.01*yy-5;
            const double h=x<10?2.:(x>=30?3.:0.);
            measured.push_back({xx,yy,float(ground+h),42});
            if(x<10) estimated.push_back({xx,yy,float(ground+99)}); // Must not double-count / 不得重複計算。
        }
        const auto r=calculateStockpileVolumes(measured,estimated,footprint,0,o);
        require(r.error.isEmpty()&&r.piles.size()==2,"Two piles not segmented");
        require(std::abs(r.totalVolume-10)<1e-5&&std::abs(r.totalArea-4)<1e-8,"Analytic prism volume failed");
        require(std::abs(r.piles[0].volume-4)<1e-5&&std::abs(r.piles[1].volume-6)<1e-5,"Per-pile volume failed");
        require(r.piles[0].estimatedCells==0,"Estimates replaced measurements");
        // 同格取最低離地高度，兩點高度 1、3 m，體積應為 0.01*1。
        // Use the minimum: heights 1 and 3 in one cell yield 0.01*1.
        auto flat=o;flat.ground={0,0,0};flat.minimumArea=0;
        const auto minimum=calculateStockpileVolumes({{.025f,.025f,1},{.075f,.075f,3}}, {},{{0,0},{.1,0},{.1,.1},{0,.1}},0,flat);
        require(minimum.error.isEmpty()&&std::abs(minimum.totalVolume-.01)<1e-8,"Cell minimum failed");
        const auto estimatedMin=calculateStockpileVolumes({},{{.025f,.025f,3},{.075f,.075f,1}},{{0,0},{.1,0},{.1,.1},{0,.1}},0,flat);
        require(estimatedMin.error.isEmpty()&&std::abs(estimatedMin.totalVolume-.01)<1e-8&&estimatedMin.piles[0].estimatedCells==1,"Synthetic cell minimum failed");
        const auto priority=calculateStockpileVolumes({{.025f,.025f,2}},{{.075f,.075f,1}},{{0,0},{.1,0},{.1,.1},{0,.1}},0,flat);
        require(priority.error.isEmpty()&&std::abs(priority.totalVolume-.02)<1e-8,"Measurement priority changed");
        // 用補點支撐缺格，仍分開統計來源；未知格不當作地面。
        // Estimates support missing cells, retaining provenance; unknowns never become ground.
        auto hole=measured;const Point removed=hole[5*40+5];hole.erase(hole.begin()+5*40+5);
        const auto unknown=calculateStockpileVolumes(hole,{},footprint,0,o);
        require(std::abs(unknown.totalVolume-9.98)<1e-5&&unknown.piles[0].enclosedUnknownCells==1,"Unknown-cell accounting failed");
        const auto filled=calculateStockpileVolumes(hole,{removed},footprint,0,o);
        require(std::abs(filled.totalVolume-10)<1e-5&&filled.piles[0].estimatedCells==1,"Fill provenance failed");
        // 自動地平面必須使用量測，且可還原緩斜地面；補點中的低離群值不得參與擬合。
        // Recover a gently sloped measured floor; synthetic low outliers never affect fitting.
        std::vector<Point> terrain;
        for(int y=0;y<100;++y) for(int x=0;x<100;++x) {
            const float xx=(x+.5f)*.1f,yy=(y+.5f)*.1f;
            const double base=.02*xx-.01*yy-5;
            terrain.push_back({xx,yy,float(base+((x>=40&&x<60&&y>=40&&y<60)?2.:0.)),42});
        }
        auto automatic=o;automatic.automaticGround=true;
        const auto autoReport=calculateStockpileVolumes(terrain,{{2,2,-100}},{{0,0},{10,0},{10,10},{0,10}},0,automatic);
        require(autoReport.error.isEmpty()&&autoReport.groundInliers>=30,"Automatic ground failed");
        require(std::abs(autoReport.ground[0]-.02)<1e-5&&std::abs(autoReport.ground[1]+.01)<1e-5&&std::abs(autoReport.ground[2]+5)<1e-4,"Sloped ground fit inaccurate");
        require(std::abs(autoReport.totalVolume-8)<1e-4,"Ground-relative automatic volume failed");
        const auto concave=calculateStockpileVolumes(measured,estimated,{{0,0},{4,0},{4,1},{3,1},{3,2},{0,2}},0,o);
        require(concave.error.isEmpty()&&std::abs(concave.totalVolume-7)<1e-5,"Concave footprint was filled as a rectangle");
        // 非格線凹角：格中心在場內仍不夠，整格跨到缺角外就不能算。
        // A non-grid-aligned notch must reject straddling cells, even when their centers are inside.
        std::vector<Point> square;
        for(int y=0;y<10;++y) for(int x=0;x<10;++x) square.push_back({(x+.5f)*.1f,(y+.5f)*.1f,1,42});
        Footprint notch{{0,0},{1,0},{1,1},{.6,1},{.6,.65},{.4,.65},{.4,1},{0,1}};
        const auto clipped=calculateStockpileVolumes(square,{{.5f,.8f,1000}},notch,0,flat);
        require(clipped.error.isEmpty()&&clipped.footprintCells==92&&std::abs(clipped.totalVolume-.92)<1e-6,"Boundary-straddling cells counted exterior area");
        require(!calculateStockpileVolumes(square,{}, {},0,flat).error.isEmpty(),"Missing footprint silently used bounding rectangle");
        // 真實倒凸輪廓：兩端料地保留，兩側缺角與場外高點不增加體積。
        // Actual concave footprint: retain end wings, reject recesses and exterior high observations.
        Footprint yard{{0,0},{19.9,0},{19.9,-13},{54.02,-13},{54.02,263},{31.45,263},{31.45,251},{19.9,251},{19.9,240},{0,240}};
        const auto site=calculateStockpileVolumes({{.05f,.05f,2},{25.05f,-12.95f,2},{25.05f,245.05f,2},{40.05f,262.05f,2},
                                                 {10.05f,-5.05f,1000},{10.05f,245.05f,1000},{25.05f,260.05f,1000}}, {},yard,0,flat);
        require(site.error.isEmpty()&&site.coveredCells==4&&std::abs(site.totalVolume-.08)<1e-6,"Inverted concave yard boundary leaked or lost wings");
        require(!calculateStockpileVolumes({{.05f,.05f,1}}, {},{{0,0},{1000,0},{1000,1000},{0,1000}},0,flat).error.isEmpty(),"Oversized fixed raster accepted");
        std::vector<Point> collinear;
        for(int i=0;i<50;++i) for(int k=0;k<3;++k) collinear.push_back({float(i+.05+k*.1),float(i+.05+k*.1),-5});
        require(!calculateStockpileVolumes(collinear,{}, {},0,automatic).error.isEmpty(),"Collinear ground support accepted");
        auto disabled=o;disabled.enabled=false;require(!calculateStockpileVolumes({}, {},{},0,disabled).enabled,"Disabled volume ran");
        require(!calculateStockpileVolumes({{0,0,0}}, {},{},0,automatic).error.isEmpty(),"Unsupported ground accepted");
        std::cout<<"PASS: exact 0.1 m cells, minima, independent piles, sloped ground, provenance, unknowns and rejection\n";
    } catch(const std::exception &e) {std::cerr<<e.what()<<"\n";return 1;}
}
