/** @file mesh_volume.cpp
 * @brief 受限 TIN 的解析體積及柱體交集；不使用體積柵格。 / Analytic bounded-TIN volume and column intersections, without volume rasters.
 */
#include "cloud.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <QFile>
#include <QSaveFile>
#include <QtEndian>
#include <map>
#include <numeric>
#include <sstream>
#include <cstring>
namespace {
namespace bg=boost::geometry;
using P=bg::model::d2::point_xy<double>;using Polygon=bg::model::polygon<P>;using Polygons=bg::model::multi_polygon<Polygon>;using Plane=std::array<double,3>;
struct XY {double x,y;};
double value(const Plane &p,double x,double y){return p[0]*x+p[1]*y+p[2];}
Polygon polygon(const std::vector<XY> &points) {Polygon p;for(auto v:points) p.outer().push_back(P(v.x,v.y));bg::correct(p);return p;}
Polygon polygon(const Footprint &points) {std::vector<XY> v;for(auto p:points) v.push_back({p[0],p[1]});return polygon(v);}
/** @brief 凸多邊形對線性高度半平面裁切。 / Clip a convex polygon to a linear-height half-plane. */
std::vector<XY> halfplane(const std::vector<XY> &p,const Plane &h) {
 std::vector<XY> out;if(p.empty()) return out;auto a=p.back();double va=value(h,a.x,a.y);
 for(auto b:p) {double vb=value(h,b.x,b.y);if((va>=0)!=(vb>=0)) {double t=va/(va-vb);out.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});}if(vb>=0) out.push_back(b);a=b;va=vb;}return out;
}
struct Integral {double area=0,volume=0;};
/** @brief 以原點平移的鞋帶矩精確積分線性高度，孔洞扣除。 / Integrate linear heights using translated shoelace moments, subtracting holes. */
template<class Ring> Integral ringIntegral(const Ring &ring,const Plane &h) {
 if(ring.size()<4) return {};
 double area2=0,mx=0,my=0;const double ox=bg::get<0>(ring[0]),oy=bg::get<1>(ring[0]);
 for(size_t i=1;i<ring.size();++i) {double ax=bg::get<0>(ring[i-1])-ox,ay=bg::get<1>(ring[i-1])-oy,bx=bg::get<0>(ring[i])-ox,by=bg::get<1>(ring[i])-oy;double c=ax*by-bx*ay;area2+=c;mx+=(ax+bx)*c;my+=(ay+by)*c;}
 if(std::abs(area2)<1e-14) return {};
 double sign=area2>0?1:-1,area=std::abs(area2)*.5;
 return {area,area*value(h,ox,oy)+sign*(h[0]*mx+h[1]*my)/6};
}
Integral integrate(const Polygon &p,const Plane &h) {auto r=ringIntegral(p.outer(),h);for(const auto &hole:p.inners()) {auto s=ringIntegral(hole,h);r.area-=s.area;r.volume-=s.volume;}return r;}
struct Face {std::array<int,3> ids;double area=0,volume=0;size_t begin=0,end=0;};
/** @brief 耳切法將已裁切上表面拆成凸三角形，供結構聯集裁切及匯出共用。
 *         Ear clipping decomposes clipped surfaces into convex triangles for shared union clipping and export.
 */
void appendSurface(const Polygon &p,const Plane &surface,std::vector<MeshSurfaceTriangle> &out) {
 if(!p.inners().empty()) throw std::runtime_error("Unexpected hole in clipped upper surface");
 std::vector<XY> v;for(size_t i=0;i+1<p.outer().size();++i) {auto a=p.outer()[i];XY q{a.x(),a.y()};if(v.empty()||std::hypot(q.x-v.back().x,q.y-v.back().y)>1e-9) v.push_back(q);}
 if(v.size()>1&&std::hypot(v.front().x-v.back().x,v.front().y-v.back().y)<1e-9) v.pop_back();
 auto cross=[](XY a,XY b,XY c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);};
 double area=0;for(size_t i=0;i<v.size();++i) area+=v[i].x*v[(i+1)%v.size()].y-v[(i+1)%v.size()].x*v[i].y;
 if(area<0) std::reverse(v.begin(),v.end());
 while(v.size()>=3) {bool found=false;for(size_t i=0;i<v.size();++i) {auto a=v[(i+v.size()-1)%v.size()],b=v[i],c=v[(i+1)%v.size()];double turn=cross(a,b,c);
    if(std::abs(turn)<1e-12) {v.erase(v.begin()+i);found=true;break;}if(turn<0) continue;
    bool contains=false;for(size_t j=0;j<v.size();++j) if(j!=i&&j!=(i+1)%v.size()&&j!=(i+v.size()-1)%v.size()) {auto q=v[j];if(cross(a,b,q)>=-1e-12&&cross(b,c,q)>=-1e-12&&cross(c,a,q)>=-1e-12) {contains=true;break;}}
    if(contains) continue;
    MeshSurfaceTriangle t;const XY points[]{a,b,c};for(int j=0;j<3;++j) t.vertices[j]={points[j].x,points[j].y,value(surface,points[j].x,points[j].y)};out.push_back(t);v.erase(v.begin()+i);found=true;break;
 }if(!found) throw std::runtime_error("Cannot triangulate clipped surface for review");}
}
}
MeshVolumeReport calculateMeshVolume(const std::vector<Point> &measured,const std::vector<Point> &estimated,const Footprint &footprint,double margin,const std::array<double,3> &ground,const VolumeOptions &o) {
 MeshVolumeReport r;r.enabled=o.enabled&&o.mesh.enabled;r.ground=ground;r.supportSpacing=o.mesh.supportSpacing;if(!r.enabled) return r;
 try {
    validateFootprint(footprint);validateColumns(o.columns);
    if(footprint.empty()||!std::isfinite(margin)||margin<0||!std::isfinite(o.mesh.supportSpacing)||o.mesh.supportSpacing<.05||!std::isfinite(o.mesh.maxEdge)||o.mesh.maxEdge<o.mesh.supportSpacing||!std::isfinite(o.mesh.maxSlope)||o.mesh.maxSlope<=0||!std::isfinite(o.pileThreshold)||o.pileThreshold<0||!std::isfinite(o.minimumArea)||o.minimumArea<0) throw std::runtime_error("Invalid mesh/footprint parameters");
    for(auto v:ground) if(!std::isfinite(v)) throw std::runtime_error("Invalid mesh ground plane");
    /** @brief 1：獨立選取實際 XY 的最低離地支點，實測優先；上限超出即報錯，不偷偷降精度。
     *         1: independently select actual-XY minimum-height supports, measurements first; fail rather than silently coarsen.
     */
    struct Support {Point p;double h;bool synthetic;};std::map<std::pair<int,int>,Support> bins;
    const auto add=[&](const std::vector<Point> &points,bool synthetic) {for(const auto &p:points) {
        if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||!insideFootprint(footprint,p.x,p.y,margin)) continue;
        if(std::abs(p.x)>1e6||std::abs(p.y)>1e6) throw std::runtime_error("Mesh coordinates exceed range");
        const std::pair<int,int> k{int(std::floor(p.x/o.mesh.supportSpacing)),int(std::floor(p.y/o.mesh.supportSpacing))};const double h=p.z-value(ground,p.x,p.y);auto it=bins.find(k);
        if(it==bins.end()) bins.emplace(k,Support{p,h,synthetic});else if(!(synthetic&&!it->second.synthetic)&&h<it->second.h) it->second={p,h,synthetic};
        if(bins.size()>o.mesh.maxVertices) throw std::runtime_error("Mesh support cap exceeded; increase --mesh-support-spacing");
    }};add(measured,false);add(estimated,true);r.vertices=bins.size();if(r.vertices<3) throw std::runtime_error("Insufficient mesh supports");
    std::vector<Point> points;for(const auto &b:bins) points.push_back(b.second.p);
    /** @brief 2：場域內縮使用幾何 buffer；每片三角形精確與凹輪廓裁切，不以凸包替代。
     *         2: geometrically inset the footprint, then exactly clip triangles against the concave region, never its convex hull.
     */
    Polygons region;auto site=polygon(footprint);
    if(margin>0) {bg::strategy::buffer::distance_symmetric<double> distance(-margin);bg::strategy::buffer::side_straight side;bg::strategy::buffer::join_round join(32);bg::strategy::buffer::end_round end(32);bg::strategy::buffer::point_circle circle(32);bg::buffer(site,region,distance,side,join,end,circle);}else region.push_back(site);
    if(region.empty()) throw std::runtime_error("Empty inset mesh footprint");
    for(const auto &p:region) r.footprintArea+=bg::area(p);
    const ColumnVolumeMask columnMask(o.columns);
    std::vector<double> xy;const auto origin=points.front();for(auto p:points) {xy.push_back(double(p.x)-origin.x);xy.push_back(double(p.y)-origin.y);}
    orgQhull::Qhull hull;std::ostringstream errors;hull.setErrorStream(&errors);hull.disableOutputStream();hull.runQhull("volume supports",2,int(points.size()),xy.data(),"d Qbb Qt Qc Qz");
    std::vector<Face> faces;std::vector<MeshSurfaceTriangle> surface;
    /** @brief 3：三角形內為線性料面；精確裁切高度門檻後，以投影面積乘平均頂點離地高度積分。
     *         3: integrate linear triangular surfaces after exact threshold clipping: projected area times mean vertex height.
     */
    for(auto f=hull.beginFacet();f!=hull.endFacet();f=f.next()) {
        if(f.isUpperDelaunay()) continue;
        auto vertices=f.vertices().toStdVector();if(vertices.size()!=3) continue;std::array<int,3> ids;bool valid=true;
        for(int j=0;j<3;++j) {ids[j]=vertices[j].point().id();if(ids[j]<0||size_t(ids[j])>=points.size()) valid=false;}if(!valid) continue;++r.inputTriangles;
        const auto a=points[ids[0]],b=points[ids[1]],c=points[ids[2]];double ux=b.x-a.x,uy=b.y-a.y,vx=c.x-a.x,vy=c.y-a.y,det=ux*vy-uy*vx;if(std::abs(det)<1e-12) continue;
        double ab=std::hypot(ux,uy),ac=std::hypot(vx,vy),bc=std::hypot(double(b.x)-c.x,double(b.y)-c.y);
        if(std::max({ab,ac,bc})>o.mesh.maxEdge||ab*ac*bc/(2*std::abs(det))>o.mesh.maxEdge*.6) {++r.rejectedDistance;continue;}
        double gx=((b.z-a.z)*vy-(c.z-a.z)*uy)/det,gy=(ux*(c.z-a.z)-vx*(b.z-a.z))/det;
        if(std::hypot(gx,gy)>o.mesh.maxSlope) {++r.rejectedSlope;continue;}
        Plane top{gx,gy,a.z-gx*a.x-gy*a.y},height{gx-ground[0],gy-ground[1],top[2]-ground[2]};
        Polygon tri=polygon(std::vector<XY>{{a.x,a.y},{b.x,b.y},{c.x,c.y}});Polygons domain;bg::intersection(tri,region,domain);for(const auto &p:domain) r.coveredArea+=integrate(p,{0,0,0}).area;
        auto threshold=height;threshold[2]-=o.pileThreshold;auto clipped=halfplane({{a.x,a.y},{b.x,b.y},{c.x,c.y}},threshold);if(clipped.size()<3) continue;
        Polygons material;bg::intersection(polygon(clipped),region,material);Face face;face.ids=ids;face.begin=surface.size();
        for(const auto &p:material) {
            auto gross=integrate(p,height);if(gross.area<1e-12||gross.volume<=0) continue;double deduction=0,excludedArea=0;
            /** @brief 使用共用結構聯集裁切，凸片上直接積分；沒有轉成體積柵格。
             *         Clip the shared structural union and integrate convex pieces directly, without rasterization.
             */
            std::vector<MeshSurfaceTriangle> convex;appendSurface(p,top,convex);
            for(const auto &triangle:convex) {
                Footprint domain;for(auto v:triangle.vertices) domain.push_back({v[0],v[1]});
                const auto excluded=columnMask.intersectSurface(domain,top,ground);deduction+=excluded.volume;excludedArea+=excluded.fullHeightArea;
            }
            deduction=std::clamp(deduction,0.,gross.volume);r.grossVolume+=gross.volume;r.columnExcludedVolume+=deduction;r.materialArea+=gross.area;face.volume+=gross.volume-deduction;face.area+=std::max(0.,gross.area-excludedArea);
            if(gross.volume-deduction>1e-12) appendSurface(p,top,surface);
        }
        face.end=surface.size();if(face.volume>1e-12) faces.push_back(face);
    }
    /** @brief 5：有效三角形依共邊連通分群、排除過小面積；與柵格分群獨立。
     *         5: group valid faces by shared edges and reject small areas independently of raster components.
     */
    std::vector<size_t> parent(faces.size());std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](size_t n) {while(parent[n]!=n) {parent[n]=parent[parent[n]];n=parent[n];}return n;};
    std::map<std::pair<int,int>,size_t> edges;
    for(size_t i=0;i<faces.size();++i) for(int j=0;j<3;++j) {
        int a=faces[i].ids[j],b=faces[i].ids[(j+1)%3];if(a>b) std::swap(a,b);auto found=edges.find({a,b});
        if(found==edges.end()) {edges[{a,b}]=i;continue;}
        // 共邊須有高於門檻的場內段；跨地面或場外的鄰接不合併。
        // Join only through an interior segment above threshold, not across ground/exterior edges.
        const auto &p=points[a],&q=points[b];double hp=p.z-value(ground,p.x,p.y)-o.pileThreshold,hq=q.z-value(ground,q.x,q.y)-o.pileThreshold;
        if(std::max(hp,hq)<=0) continue;
        double lo=0,hi=1;if(hp<0) lo=hp/(hp-hq);if(hq<0) hi=hp/(hp-hq);
        double t=(lo+hi)*.5,x=p.x+t*(q.x-p.x),y=p.y+t*(q.y-p.y);if(!insideFootprint(footprint,x,y,margin)) continue;
        auto u=root(i),v=root(found->second);if(u!=v) parent[u]=v;
    }
    struct Component {double area=0,volume=0;};std::map<size_t,Component> components;for(size_t i=0;i<faces.size();++i) {auto &c=components[root(i)];c.area+=faces[i].area;c.volume+=faces[i].volume;}
    for(const auto &c:components) {if(c.second.area+1e-9<o.minimumArea) {++r.discardedComponents;r.discardedSmallVolume+=c.second.volume;}else {++r.components;r.totalVolume+=c.second.volume;r.totalArea+=c.second.area;}}
    for(size_t i=0;i<faces.size();++i) if(components[root(i)].area+1e-9>=o.minimumArea) r.surface.insert(r.surface.end(),surface.begin()+faces[i].begin,surface.begin()+faces[i].end);
    if(r.inputTriangles==0||r.coveredArea<=0) throw std::runtime_error("No supported mesh triangles inside site");
 }catch(const std::exception &e) {r.error=QString::fromUtf8(e.what());r.surface.clear();}
 return r;
}
void exportMeshSurface(const QString &path,const MeshVolumeReport &r) {
 if(!r.enabled||!r.error.isEmpty()||r.surface.empty()) throw std::runtime_error("No valid mesh surface to export");
 if(QFile::exists(path)) throw std::runtime_error("Output mesh exists; choose a new name");
 QSaveFile file(path);if(!file.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot create mesh output");
 QByteArray header="ply\nformat binary_little_endian 1.0\ncomment calibrated building frame; meters; origin at reference roof corner\ncomment OPEN upper surface: not watertight; column solids subtracted analytically, not cut out of this review surface\n";
 header+="comment ground z="+QByteArray::number(r.ground[0],'g',17)+"*x+"+QByteArray::number(r.ground[1],'g',17)+"*y+"+QByteArray::number(r.ground[2],'g',17)+"\n";
 header+="element vertex "+QByteArray::number(qulonglong(r.surface.size()*3))+"\nproperty double x\nproperty double y\nproperty double z\nelement face "+QByteArray::number(qulonglong(r.surface.size()))+"\nproperty list uchar int vertex_indices\nend_header\n";
 auto write=[&](const char *p,qint64 n){if(file.write(p,n)!=n) throw std::runtime_error("Mesh write failed");};write(header.constData(),header.size());
 for(const auto &t:r.surface) for(const auto &p:t.vertices) for(double value:p) {quint64 bits;std::memcpy(&bits,&value,8);bits=qToLittleEndian(bits);write(reinterpret_cast<const char *>(&bits),8);}
 for(size_t i=0;i<r.surface.size();++i) {const char count=3;write(&count,1);for(int j=0;j<3;++j) {qint32 id=qToLittleEndian(qint32(i*3+j));write(reinterpret_cast<const char *>(&id),4);}}
 if(!file.commit()) throw std::runtime_error("Mesh commit failed");
}
