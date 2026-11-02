/** @file mesh_volume_test.cpp
 * @brief 解析料面積分、凹邊界、柱頂裁切與獨立小區測試。 / Analytic surface integration, concave boundaries, column tops and small components.
 */
#include "cloud.h"
#include <QTemporaryDir>
#include <QFile>
#include <cmath>
#include <iostream>
void require(bool ok,const char *message){if(!ok) throw std::runtime_error(message);}
int main(){try {
 VolumeOptions o;o.enabled=true;o.automaticGround=false;o.pileThreshold=0;o.minimumArea=0;o.mesh.supportSpacing=.1;o.mesh.maxEdge=2;
 Footprint site{{0,0},{1,0},{1,1},{0,1}};std::vector<Point> points{{0,0,2},{1,0,2},{1,1,2},{0,1,2}};
 auto calc=[&](const std::vector<Point> &p){auto r=calculateMeshVolume(p,{},site,0,{0,0,0},o);if(!r.error.isEmpty()) throw std::runtime_error(r.error.toStdString());return r;};
 auto r=calc(points);require(std::abs(r.totalVolume-2)<1e-8&&std::abs(r.coveredArea-1)<1e-8,"Flat mesh volume failed");
 auto slope=points;for(auto &p:slope) p.z=1+p.x+2*p.y;r=calc(slope);require(std::abs(r.totalVolume-2.5)<1e-8,"Linear surface integral failed");
 auto slopedGround=calculateMeshVolume(slope,{},site,0,{.2,.1,-1},o);require(slopedGround.error.isEmpty()&&std::abs(slopedGround.totalVolume-3.35)<1e-8,"Sloped ground mismatch");
 Footprint notch{{0,0},{1,0},{1,1},{.6,1},{.6,.6},{.4,.6},{.4,1},{0,1}};
 r=calculateMeshVolume(points,{},notch,0,{0,0,0},o);require(r.error.isEmpty()&&std::abs(r.totalVolume-1.84)<1e-7,"Concave clipping filled a notch");
 o.pileThreshold=.2;auto crossing=points;for(auto &p:crossing) p.z=p.x;r=calc(crossing);require(std::abs(r.totalVolume-.48)<1e-8,"Threshold clipping failed");o.pileThreshold=0;
 ColumnSemantic col;col.id=1;col.state=ColumnState::Confirmed;col.top=10;col.footprint=columnSection("rectangle",.25,.75,.25,.75);o.columns={col};r=calc(points);require(std::abs(r.totalVolume-1.5)<1e-8,"Column XY intersection failed");
 o.columns[0].top=1;r=calc(points);require(std::abs(r.totalVolume-1.75)<1e-8,"Buried short column failed");
 o.columns[0].footprint=site;o.columns[0].top=2;auto crossingTop=points;for(auto &p:crossingTop) p.z=1+2*p.x;r=calc(crossingTop);require(std::abs(r.totalVolume-.25)<1e-7,"Column top crossing the surface failed");
 o.columns[0].top=0;auto high=points;for(auto &p:high) p.z=3;auto crossingFloor=calculateMeshVolume(high,{},site,0,{2,0,-1},o);require(crossingFloor.error.isEmpty()&&std::abs(crossingFloor.totalVolume-2.75)<1e-7,"Column top crossing ground failed");
 o.columns[0].top=-1;r=calc(points);require(std::abs(r.totalVolume-2)<1e-8,"Below-ground column deducted");
 o.columns[0].state=ColumnState::Rejected;r=calc(points);require(std::abs(r.totalVolume-2)<1e-8,"Rejected column deducted");
 o.columns[0].state=ColumnState::Confirmed;o.columns[0].top=10;o.columns[0].shape="polygon";o.columns[0].footprint={{0,0},{1,0},{1,.2},{.6,.2},{.6,.8},{1,.8},{1,1},{0,1},{0,.8},{.4,.8},{.4,.2},{0,.2}};
 r=calc(points);require(std::abs(r.totalVolume-.96)<1e-8,"I-section treated as solid bounding box");
 // 凹箱體與高柱相接，各自頂高且不可重複扣除。 / Notched box and taller pier share edges, retaining independent tops.
 ColumnSemantic box=col;box.id=2;box.shape="polygon";box.top=1;
 box.footprint={{0,0},{1,0},{1,1},{.6,1},{.6,.5},{.4,.5},{.4,1},{0,1}};
 auto pier=col;pier.id=3;pier.top=3;pier.footprint=columnSection("rectangle",.4,.6,.5,1);
 o.columns={box,pier};r=calc(points);
 require(std::abs(r.columnExcludedVolume-1.1)<1e-7&&std::abs(r.totalVolume-.9)<1e-7,"Notched box/pier mesh deduction failed");
 // 梁底／梁頂穿越料面時解析裁切，保留架空區。 / Clip elevated beam intervals against a varying surface.
 auto elevated=col;elevated.groundRelative=true;elevated.bottomHeight=1;elevated.topHeight=1.5;elevated.footprint=site;o.columns={elevated};
 r=calc(points);require(std::abs(r.columnExcludedVolume-.5)<1e-8&&std::abs(r.totalVolume-1.5)<1e-8&&std::abs(r.totalArea-1)<1e-8,"Mesh beam removed clearance/projected area");
 auto varying=points;for(auto &p:varying) p.z=.5+2*p.x;r=calc(varying);require(std::abs(r.columnExcludedVolume-.3125)<1e-7,"Beam bottom/top crossing surface failed");
 auto translated=calculateMeshVolume(points,{},site,0,{.1,.2,-3},o);require(translated.error.isEmpty()&&std::abs(translated.columnExcludedVolume-.5)<1e-8,"Relative mesh beam failed on translated/sloped ground");
 // 重疊柱、箱、梁聯集與斜料面解析結果。 / Analytic overlapping column/box/beam union under a sloping surface.
 auto post=col;post.id=900;post.footprint=columnSection("rectangle",0,.6,0,1);post.top=2;
 auto boxSolid=post;boxSolid.id=901;boxSolid.footprint=columnSection("rectangle",.4,1,0,1);boxSolid.top=3;
 auto beamSolid=post;beamSolid.id=902;beamSolid.groundRelative=true;beamSolid.bottomHeight=1;beamSolid.topHeight=2.5;beamSolid.footprint=site;
 o.columns={post,boxSolid,beamSolid};auto tall=points;for(auto &p:tall) p.z=4;
 r=calc(tall);require(std::abs(r.totalVolume-1.2)<1e-7&&std::abs(r.columnExcludedVolume-2.8)<1e-7&&std::abs(r.totalArea-1)<1e-7,"Mesh triple union failed");
 for(auto &p:tall) p.z=1+3*p.x;r=calc(tall);
 require(std::abs(r.totalVolume-1./6)<1e-7&&std::abs(r.totalArea-1./3)<1e-7,"Union clipping against sloped material surface failed");
 o.columns={beamSolid,boxSolid,post};auto reordered=calc(tall);require(std::abs(reordered.totalVolume-r.totalVolume)<1e-9,"Mesh union depends on input order");
 o.columns.clear();r=calculateMeshVolume({{.1f,.1f,2},{.9f,.1f,2},{.9f,.9f,2},{.1f,.9f,2}}, {},site,.1,{0,0,0},o);require(r.error.isEmpty()&&std::abs(r.totalVolume-1.28)<1e-6,"Inset footprint failed");
 auto duplicates=points;duplicates.push_back({0,0,99});r=calc(duplicates);require(std::abs(r.totalVolume-2)<1e-8,"Minimum measured support failed");
 r=calculateMeshVolume(points,{{0,0,-10}},site,0,{0,0,0},o);require(r.error.isEmpty()&&std::abs(r.totalVolume-2)<1e-8,"Synthetic support replaced measured support");
 auto islands=points;for(auto p:points) {p.x+=10;islands.push_back(p);}Footprint big{{0,0},{11,0},{11,1},{0,1}};o.minimumArea=1.1;
 r=calculateMeshVolume(islands,{},big,0,{0,0,0},o);require(r.error.isEmpty()&&r.components==0&&r.discardedComponents==2&&std::abs(r.discardedSmallVolume-4)<1e-8,"Long gap or small component filter failed");o.minimumArea=0;
 r=calc(points);require(std::abs(r.grossVolume-r.columnExcludedVolume-r.discardedSmallVolume-r.totalVolume)<1e-8,"Mesh volume accounting failed");
 QTemporaryDir tmp;const auto path=tmp.filePath("surface.ply");exportMeshSurface(path,r);QFile file(path);require(file.open(QIODevice::ReadOnly)&&file.readAll().contains("element face 2\n"),"Mesh export missing faces");
 bool refused=false;try {exportMeshSurface(path,r);}catch(const std::exception &){refused=true;}require(refused,"Mesh export overwrote file");
 require(!calculateMeshVolume(points,{}, {},0,{0,0,0},o).error.isEmpty(),"Empty site accepted");o.mesh.enabled=false;require(!calculateMeshVolume(points,{},site,0,{0,0,0},o).enabled,"Disabled mesh ran");
 std::cout<<"PASS: analytic mesh volumes, ground, concavity, thresholds, columns, support priority, gaps, small components and mesh export\n";
 }catch(const std::exception &e){std::cerr<<e.what()<<"\n";return 1;}return 0;}
