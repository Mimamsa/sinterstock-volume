/** @file columns.cpp
 * @brief 幾何柱體提議、人工審核檔及柱料體積交集。 / Geometric proposals, reviewed models, and column/material volume intersections.
 */
#include "cloud.h"
#include <QFile>
#include <QSaveFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

namespace {
using SectionTriangle=std::array<std::array<double,2>,3>;
/** @brief 將凹柱截面耳切為逆時針三角形，供重疊驗證與積分共用。
 *         Ear-clip a concave section into CCW triangles shared by overlap validation and integration.
 */
std::vector<SectionTriangle> triangulateSection(Footprint p) {
    std::vector<SectionTriangle> result;
    auto turn=[](const auto &a,const auto &b,const auto &c) {return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);};
    double area=0;for(size_t i=0;i<p.size();++i) area+=p[i][0]*p[(i+1)%p.size()][1]-p[(i+1)%p.size()][0]*p[i][1];
    if(area<0) std::reverse(p.begin(),p.end());
    while(p.size()>3) {
        bool found=false;
        for(size_t i=0;i<p.size();++i) {
            const size_t previous=(i+p.size()-1)%p.size(),next=(i+1)%p.size();
            const auto a=p[previous],b=p[i],c=p[next];if(turn(a,b,c)<=1e-12) continue;
            bool contains=false;
            for(size_t j=0;j<p.size();++j) if(j!=previous&&j!=i&&j!=next&&turn(a,b,p[j])>=-1e-12&&turn(b,c,p[j])>=-1e-12&&turn(c,a,p[j])>=-1e-12) {contains=true;break;}
            if(contains) continue;
            result.push_back({a,b,c});p.erase(p.begin()+i);found=true;break;
        }
        if(!found) throw std::runtime_error("Cannot triangulate column section; remove redundant/collinear polygon vertices");
    }
    if(p.size()==3&&turn(p[0],p[1],p[2])>1e-12) result.push_back({p[0],p[1],p[2]});
    return result;
}

}

std::array<double,4> columnBounds(const ColumnSemantic &c) {
    std::array<double,4> b{1e30,-1e30,1e30,-1e30};
    for(const auto &p:c.footprint) {b[0]=std::min(b[0],p[0]);b[1]=std::max(b[1],p[0]);b[2]=std::min(b[2],p[1]);b[3]=std::max(b[3],p[1]);}return b;
}
Footprint columnSection(const QString &shape,double x0,double x1,double y0,double y1) {
    if(!std::isfinite(x0)||!std::isfinite(x1)||!std::isfinite(y0)||!std::isfinite(y1)||x1<=x0||y1<=y0) throw std::runtime_error("Invalid column bounds");
    if(shape=="rectangle") return {{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
    if(shape!="circle") throw std::runtime_error("Custom polygons must be edited in the semantic JSON");
    if(std::abs((x1-x0)-(y1-y0))>1e-5) throw std::runtime_error("Circle requires equal X/Y diameters");
    Footprint p;for(int i=0;i<64;++i) {const double a=i*2*3.141592653589793/64;p.push_back({(x0+x1)*.5+(x1-x0)*.5*std::cos(a),(y0+y1)*.5+(y1-y0)*.5*std::sin(a)});}return p;
}
double columnBottomZ(const ColumnSemantic &c,double x,double y) {return c.referenceGround[0]*x+c.referenceGround[1]*y+c.referenceGround[2]+(c.groundRelative?c.bottomHeight:0);}
double columnTopZ(const ColumnSemantic &c,double x,double y) {return c.groundRelative?c.referenceGround[0]*x+c.referenceGround[1]*y+c.referenceGround[2]+c.topHeight:c.top;}
void setColumnGround(std::vector<ColumnSemantic> &columns,const std::array<double,3> &ground) {
    for(double v:ground) if(!std::isfinite(v)) throw std::runtime_error("Invalid structural ground reference");
    for(auto &c:columns) c.referenceGround=ground;
}
void validateColumns(const std::vector<ColumnSemantic> &columns) {
    if(columns.size()>1000) throw std::runtime_error("Too many column semantics (maximum 1000)");
    std::set<int> ids;
    for(size_t i=0;i<columns.size();++i) {
        const auto &c=columns[i];if(c.id<=0||c.id>1000000||!ids.insert(c.id).second||!std::isfinite(c.top)||std::abs(c.top)>1e5) throw std::runtime_error("Invalid/duplicate column ID or top");
        if(c.groundRelative&&(!std::isfinite(c.bottomHeight)||!std::isfinite(c.topHeight)||c.bottomHeight<0||c.topHeight<=c.bottomHeight||c.topHeight>1000)) throw std::runtime_error("Invalid ground-relative structural bounds");
        for(double v:c.referenceGround) if(!std::isfinite(v)) throw std::runtime_error("Invalid structural ground reference");
        validateFootprint(c.footprint);if(!c.evidenceFootprint.empty()) validateFootprint(c.evidenceFootprint);if(c.footprint.empty()) throw std::runtime_error("Empty column section");
        if(c.shape!="rectangle"&&c.shape!="circle"&&c.shape!="polygon") throw std::runtime_error("Unknown column shape");
        const auto a=columnBounds(c);
        if(a[1]-a[0]>20||a[3]-a[2]>20) throw std::runtime_error("Column section exceeds 20 m");

    }
}
std::vector<ColumnSemantic> detectColumnCandidates(const std::vector<Point> &points,const ColumnOptions &o) {
    if(!o.enabled) return {};
    if(!std::isfinite(o.xyCell)||o.xyCell<.1||!std::isfinite(o.zCell)||o.zCell<.1||!std::isfinite(o.minHeight)||o.minHeight<1 ||
       !std::isfinite(o.minCoverage)||o.minCoverage<=0||o.minCoverage>1||!std::isfinite(o.maxWidth)||o.maxWidth<o.xyCell||o.maxWidth>10) throw std::runtime_error("Invalid column detection settings");
    using Key=std::pair<int,int>;
    struct Tile {std::set<int> levels;size_t count=0;double low=1e30,high=-1e30;bool visited=false;};
    std::map<Key,Tile> tiles;
    /** @brief 步驟 1：XY 窄格中的高度佔用須跨越至少 4 m 且持續，料堆斜面通常不符合此條件。
     *         Step 1: require persistent height occupancy over at least 4 m in narrow XY bins, unlike ordinary pile slopes.
     * 只用量測；牆面會形成過長連通片，於下一步拒絕，短牆／吊掛物仍可能成為待審候選。
     * Use measurements only; long connected walls are rejected next, while short walls/hanging objects may remain proposals.
     */
    for(const auto &p:points) {
        if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||std::max({std::abs(p.x),std::abs(p.y),std::abs(p.z)})>1e5) continue;
        auto &t=tiles[{int(std::floor(p.x/o.xyCell)),int(std::floor(p.y/o.xyCell))}];
        ++t.count;t.low=std::min(t.low,double(p.z));t.high=std::max(t.high,double(p.z));t.levels.insert(int(std::floor(p.z/o.zCell)));
    }
    for(auto it=tiles.begin();it!=tiles.end();) {
        const auto &t=it->second;const double occupancy=t.levels.size()/(std::floor(t.high/o.zCell)-std::floor(t.low/o.zCell)+1);
        if(t.high-t.low<o.minHeight||occupancy<o.minCoverage||t.count<20) it=tiles.erase(it);else ++it;
    }
    std::vector<ColumnSemantic> result;
    for(auto &entry:tiles) if(!entry.second.visited) {
        std::vector<Key> queue{entry.first};entry.second.visited=true;
        int xmin=entry.first.first,xmax=xmin,ymin=entry.first.second,ymax=ymin;
        ColumnSemantic c;c.observedBottom=1e30;c.observedTop=-1e30;double occupancy=0;
        for(size_t q=0;q<queue.size();++q) {
            const auto k=queue[q];auto &t=tiles.at(k);xmin=std::min(xmin,k.first);xmax=std::max(xmax,k.first);ymin=std::min(ymin,k.second);ymax=std::max(ymax,k.second);
            c.evidencePoints+=t.count;c.observedBottom=std::min(c.observedBottom,t.low);c.observedTop=std::max(c.observedTop,t.high);
            occupancy+=t.levels.size()/(std::floor(t.high/o.zCell)-std::floor(t.low/o.zCell)+1);
            for(int dy=-1;dy<=1;++dy) for(int dx=-1;dx<=1;++dx) {auto it=tiles.find({k.first+dx,k.second+dy});if(it!=tiles.end()&&!it->second.visited) {it->second.visited=true;queue.push_back(it->first);}}
        }
        const double w=(xmax-xmin+1)*o.xyCell,d=(ymax-ymin+1)*o.xyCell;
        if(w>o.maxWidth||d>o.maxWidth||std::max(w,d)/std::min(w,d)>4) continue;
        c.id=int(result.size()+1);c.top=c.observedTop;c.verticalCoverage=occupancy/queue.size();
        c.footprint=columnSection("rectangle",xmin*o.xyCell,(xmax+1)*o.xyCell,ymin*o.xyCell,(ymax+1)*o.xyCell);
        result.push_back(c);if(result.size()>=1000) break;
    }
    return result;
}
std::vector<ColumnSemantic> loadColumnSemantics(const QString &path,const QString &siteId) {
    QFile f(path);if(!f.exists()) return {};if(!f.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot read column semantics");
    QJsonParseError error;const auto doc=QJsonDocument::fromJson(f.readAll(),&error);const auto root=doc.object();
    if(error.error!=QJsonParseError::NoError||root.value("version").toInt()!=1||root.value("site_id").toString()!=siteId||root.value("base").toString("ground")!="ground"||!root.value("columns").isArray()) throw std::runtime_error("Invalid column semantics version/site ID");
    std::vector<ColumnSemantic> columns;
    for(const auto &v:root.value("columns").toArray()) {
        const auto a=v.toObject();ColumnSemantic c;c.id=a.value("id").toInt();c.shape=a.value("shape").toString("rectangle");
        const auto state=a.value("state").toString();if(state=="confirmed") c.state=ColumnState::Confirmed;else if(state=="rejected") c.state=ColumnState::Rejected;else if(state!="pending") throw std::runtime_error("Invalid column review state");
        if(!a.value("top_z").isDouble()) throw std::runtime_error("Missing column top_z");
        c.top=a.value("top_z").toDouble();c.kind=a.value("kind").toString("column");
        if(a.contains("height_above_ground_m")) {
            if(!a.value("height_above_ground_m").isDouble()||!a.value("bottom_height_above_ground_m").isDouble()) throw std::runtime_error("Invalid relative structural heights");
            c.groundRelative=true;c.topHeight=a.value("height_above_ground_m").toDouble();c.bottomHeight=a.value("bottom_height_above_ground_m").toDouble();
        }
        if(a.contains("reference_ground_plane")||root.contains("reference_ground_plane")) {
            const auto plane=(a.contains("reference_ground_plane")?a.value("reference_ground_plane"):root.value("reference_ground_plane")).toArray();if(plane.size()!=3) throw std::runtime_error("Invalid reference ground plane");
            for(int i=0;i<3;++i) {if(!plane[i].isDouble()) throw std::runtime_error("Invalid reference ground value");c.referenceGround[i]=plane[i].toDouble();}
        }
        for(const auto &value:a.value("footprint_xy_m").toArray()) {const auto p=value.toArray();if(p.size()!=2||!p[0].isDouble()||!p[1].isDouble()) throw std::runtime_error("Invalid column polygon");c.footprint.push_back({p[0].toDouble(),p[1].toDouble()});}
        c.provenance=a.value("provenance").toString();c.group=a.value("group").toString();
        for(const auto &id:a.value("legacy_ids").toArray()) {if(id.toInt()<=0) throw std::runtime_error("Invalid legacy ID");c.legacyIds.push_back(id.toInt());}
        if(root.value("finalized").toBool()&&(c.state!=ColumnState::Confirmed||c.group.isEmpty())) throw std::runtime_error("Final model requires confirmed grouped columns");
        for(const auto &value:a.value("evidence_footprint_xy_m").toArray()) {const auto p=value.toArray();if(p.size()!=2||!p[0].isDouble()||!p[1].isDouble()) throw std::runtime_error("Invalid evidence polygon");c.evidenceFootprint.push_back({p[0].toDouble(),p[1].toDouble()});}
        c.observedBottom=a.value("observed_bottom_z").toDouble(c.top);c.observedTop=a.value("observed_top_z").toDouble(c.top);
        c.verticalCoverage=a.value("vertical_coverage").toDouble();c.evidencePoints=size_t(std::max(0.,a.value("evidence_points").toDouble()));columns.push_back(c);
    }
    validateColumns(columns);return columns;
}
bool columnSemanticsFinalized(const QString &path,const QString &siteId) {
    QFile f(path);if(!f.exists()) return false;if(!f.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot read final column flag");
    const auto root=QJsonDocument::fromJson(f.readAll()).object();if(root.value("site_id").toString()!=siteId||root.value("version").toInt()!=1) throw std::runtime_error("Invalid finalized model site/version");
    return root.value("finalized").toBool(false);
}
void saveColumnSemantics(const QString &path,const QString &siteId,const std::vector<ColumnSemantic> &columns,bool finalized) {
    validateColumns(columns);const ColumnVolumeMask checked(columns);(void)checked;QJsonArray array;
    for(const auto &c:columns) {
        if(finalized&&(c.state!=ColumnState::Confirmed||c.group.isEmpty())) throw std::runtime_error("Final model requires confirmed grouped columns");
        QJsonArray legacy;for(auto id:c.legacyIds) legacy.append(id);
        QJsonArray polygon;for(const auto &p:c.footprint) polygon.append(QJsonArray{p[0],p[1]});
        QJsonArray evidence;for(const auto &p:c.evidenceFootprint) evidence.append(QJsonArray{p[0],p[1]});
        QJsonObject item{{"group",c.group},{"legacy_ids",legacy},{"provenance",c.provenance},{"evidence_footprint_xy_m",evidence},{"id",c.id},{"state",c.state==ColumnState::Confirmed?"confirmed":(c.state==ColumnState::Rejected?"rejected":"pending")},
            {"shape",c.shape},{"footprint_xy_m",polygon},{"top_z",c.top},{"observed_bottom_z",c.observedBottom},{"observed_top_z",c.observedTop},
            {"vertical_coverage",c.verticalCoverage},{"evidence_points",double(c.evidencePoints)}};
        item["kind"]=c.kind;item["reference_ground_plane"]=QJsonArray{c.referenceGround[0],c.referenceGround[1],c.referenceGround[2]};
        if(c.groundRelative) {item["height_above_ground_m"]=c.topHeight;item["bottom_height_above_ground_m"]=c.bottomHeight;}
        array.append(item);}
    QSaveFile f(path);if(!f.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot write column semantics");
    const auto bytes=QJsonDocument(QJsonObject{{"version",1},{"finalized",finalized},{"site_id",siteId},{"base","ground"},{"columns",array}}).toJson();
    if(f.write(bytes)!=bytes.size()||!f.commit()) throw std::runtime_error("Cannot commit column semantics");
}
size_t removeColumnPoints(std::vector<Point> &points,const std::vector<ColumnSemantic> &columns) {
    if(std::none_of(columns.begin(),columns.end(),[](const auto &c) {return c.state==ColumnState::Confirmed;})) return 0;
    const ColumnVolumeMask mask(columns);const auto n=points.size();
    points.erase(std::remove_if(points.begin(),points.end(),[&](const Point &p) {return mask.contains(p.x,p.y,p.z);}),points.end());return n-points.size();
}
std::vector<Point> sampleColumnEvidence(const std::vector<Point> &points,const std::vector<ColumnSemantic> &columns,size_t limit) {
    std::map<std::pair<int,int>,std::vector<size_t>> bins;
    std::vector<std::array<double,4>> bounds;
    for(size_t i=0;i<columns.size();++i) {auto b=columnBounds(columns[i]);
        for(const auto &p:columns[i].evidenceFootprint) {b[0]=std::min(b[0],p[0]);b[1]=std::max(b[1],p[0]);b[2]=std::min(b[2],p[1]);b[3]=std::max(b[3],p[1]);}
        bounds.push_back(b);
        for(int y=int(std::floor(b[2]-.05));y<=int(std::floor(b[3]+.05));++y) for(int x=int(std::floor(b[0]-.05));x<=int(std::floor(b[1]+.05));++x) bins[{x,y}].push_back(i);}
    std::vector<Point> result;size_t seen=0;std::mt19937 random(481);
    if(!limit) return result;
    for(const auto &p:points) {if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||std::abs(p.x)>1e6||std::abs(p.y)>1e6) continue;auto it=bins.find({int(std::floor(p.x)),int(std::floor(p.y))});if(it==bins.end()) continue;
        bool match=false;for(auto i:it->second) {const auto &b=bounds[i];if(p.x>=b[0]-.05&&p.x<=b[1]+.05&&p.y>=b[2]-.05&&p.y<=b[3]+.05&&p.z<=columnTopZ(columns[i],p.x,p.y)+.05&&(!columns[i].groundRelative||p.z>=columnBottomZ(columns[i],p.x,p.y)-.05)) {match=true;break;}}
        if(!match) continue;
        ++seen;
        if(result.size()<limit) result.push_back(p);else {std::uniform_int_distribution<size_t> pick(0,seen-1);const auto id=pick(random);if(id<limit) result[id]=p;}
    }return result;
}

namespace {
double cross(const std::array<double,2> &a,const std::array<double,2> &b,const std::array<double,2> &c) {return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);}

}
ColumnVolumeMask::ColumnVolumeMask(const std::vector<ColumnSemantic> &columns) {
    validateColumns(columns);
    /** @brief 耳切三角化支援工字形等凹實心截面；不以包圍盒代替實際形狀。
     *         Ear clipping supports concave solid sections such as I-beams; bounding boxes never replace actual sections.
     */
    for(const auto &c:columns) if(c.state==ColumnState::Confirmed) {
        for(const auto &t:triangulateSection(c.footprint)) triangles.push_back({t[0],t[1],t[2],c.top,c.id,c.groundRelative,c.bottomHeight,c.topHeight,c.referenceGround});
    }
    if(triangles.empty()) return;
    double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30;
    for(const auto &t:triangles) for(const auto &p:{t.a,t.b,t.c}) {xmin=std::min(xmin,p[0]);xmax=std::max(xmax,p[0]);ymin=std::min(ymin,p[1]);ymax=std::max(ymax,p[1]);}
    x0=int(std::floor(xmin));y0=int(std::floor(ymin));nx=int(std::floor(xmax))-x0+1;ny=int(std::floor(ymax))-y0+1;
    if(double(nx)*ny>1000000) throw std::runtime_error("Column index exceeds 1,000,000 square meters");
    bins.resize(size_t(nx)*ny);
    for(size_t i=0;i<triangles.size();++i) {const auto &t=triangles[i];const int xa=int(std::floor(std::min({t.a[0],t.b[0],t.c[0]})))-x0,xb=int(std::floor(std::max({t.a[0],t.b[0],t.c[0]})))-x0;
        const int ya=int(std::floor(std::min({t.a[1],t.b[1],t.c[1]})))-y0,yb=int(std::floor(std::max({t.a[1],t.b[1],t.c[1]})))-y0;
        for(int y=ya;y<=yb;++y) for(int x=xa;x<=xb;++x) bins[size_t(y)*nx+x].push_back(i);}
}
namespace {
using HeightPlane=std::array<double,3>;
struct MaterialPrism {Footprint xy;HeightPlane bottom,top;};
double planeAt(const HeightPlane &h,const std::array<double,2> &p) {return h[0]*p[0]+h[1]*p[1]+h[2];}
HeightPlane difference(const HeightPlane &a,const HeightPlane &b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
/** @brief 平移鞋帶矩：積分線性高度，避免大座標相消。 / Translated shoelace moments integrate linear heights stably. */
double planeIntegral(const Footprint &p,const HeightPlane &h) {
    if(p.size()<3) return 0;
    const auto origin=p[0];double twice=0,mx=0,my=0;
    for(size_t i=0;i<p.size();++i) {const auto &a=p[i],&b=p[(i+1)%p.size()];const double ax=a[0]-origin[0],ay=a[1]-origin[1],bx=b[0]-origin[0],by=b[1]-origin[1],k=ax*by-bx*ay;twice+=k;mx+=(ax+bx)*k;my+=(ay+by)*k;}
    const double sign=twice<0?-1:1;return std::abs(twice)*.5*planeAt(h,origin)+sign*(h[0]*mx+h[1]*my)/6;
}
double projectedArea(const Footprint &p) {return planeIntegral(p,{0,0,1});}
double prismVolume(const MaterialPrism &p) {return planeIntegral(p.xy,difference(p.top,p.bottom));}
/** @brief 用 XY 半平面分割凸多邊形；共面情況只歸一側，不製造重複片。
 *         Split a convex polygon by an XY half-plane; coplanar pieces belong to only one side.
 */
std::pair<Footprint,Footprint> splitXY(const Footprint &p,const HeightPlane &line) {
    if(p.empty()) return {};
    double low=1e300,high=-1e300;for(auto v:p) {double d=planeAt(line,v);low=std::min(low,d);high=std::max(high,d);}
    if(low>=-1e-10) return {p,{}};
    if(high<=1e-10) return {{},p};
    Footprint positive,negative;auto a=p.back();double da=planeAt(line,a);
    for(auto b:p) {double db=planeAt(line,b);if((da>=0)!=(db>=0)) {double t=da/(da-db);std::array<double,2> q{a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1])};positive.push_back(q);negative.push_back(q);}if(db>=0) positive.push_back(b);else negative.push_back(b);a=b;da=db;}
    if(projectedArea(positive)<1e-12) positive.clear();
    if(projectedArea(negative)<1e-12) negative.clear();
    return {positive,negative};
}
void keepPrism(std::vector<MaterialPrism> &out,MaterialPrism p) {if(prismVolume(p)>1e-12) out.push_back(std::move(p));}
/** @brief 以高度平面分割有界料柱，精確保留平面穿越上下界時的楔形部分。
 *         Split a bounded prism at a height plane, preserving wedges where the plane crosses its bounds.
 */
std::pair<std::vector<MaterialPrism>,std::vector<MaterialPrism>> splitHeight(const MaterialPrism &p,const HeightPlane &height) {
    std::vector<MaterialPrism> below,above;
    auto topSplit=splitXY(p.xy,difference(height,p.top));
    if(!topSplit.first.empty()) keepPrism(below,{topSplit.first,p.bottom,p.top});
    if(!topSplit.second.empty()) {
        auto bottomSplit=splitXY(topSplit.second,difference(height,p.bottom));
        if(!bottomSplit.second.empty()) keepPrism(above,{bottomSplit.second,p.bottom,p.top});
        if(!bottomSplit.first.empty()) {keepPrism(below,{bottomSplit.first,p.bottom,height});keepPrism(above,{bottomSplit.first,height,p.top});}
    }
    return {below,above};
}
/** @brief 逐半空間做三維差集；輸出互不重疊的凸投影有界料柱，不以整格上限掩蓋重扣。
 *         Subtract a convex solid by half-spaces, producing disjoint bounded prisms without hiding duplicate deductions via caps.
 */
std::vector<MaterialPrism> subtractSolid(const MaterialPrism &p,const MaterialPrism &mask) {
    std::vector<MaterialPrism> outside;Footprint inside=p.xy;
    for(size_t i=0;i<mask.xy.size()&&!inside.empty();++i) {
        const auto &a=mask.xy[i],&b=mask.xy[(i+1)%mask.xy.size()];HeightPlane side{a[1]-b[1],b[0]-a[0],(b[1]-a[1])*a[0]-(b[0]-a[0])*a[1]};
        auto parts=splitXY(inside,side);if(!parts.second.empty()) keepPrism(outside,{parts.second,p.bottom,p.top});inside=std::move(parts.first);
    }
    if(!inside.empty()) {
        auto lower=splitHeight({inside,p.bottom,p.top},mask.bottom);outside.insert(outside.end(),lower.first.begin(),lower.first.end());
        for(const auto &aboveLower:lower.second) {auto upper=splitHeight(aboveLower,mask.top);outside.insert(outside.end(),upper.second.begin(),upper.second.end());}
    }
    return outside;
}
/** @brief 剩餘料體投影取聯集，垂直堆疊殘片不重複計面積。 / Union residual projections without counting vertically stacked pieces twice. */
double residualProjectionArea(const std::vector<MaterialPrism> &prisms) {
    std::vector<Footprint> disjoint;
    for(const auto &prism:prisms) {
        std::vector<Footprint> pieces{prism.xy};
        for(const auto &mask:disjoint) {
            std::vector<Footprint> next;
            for(auto p:pieces) for(size_t i=0;i<mask.size()&&!p.empty();++i) {const auto &a=mask[i],&b=mask[(i+1)%mask.size()];auto cut=splitXY(p,{a[1]-b[1],b[0]-a[0],(b[1]-a[1])*a[0]-(b[0]-a[0])*a[1]});if(!cut.second.empty()) next.push_back(std::move(cut.second));p=std::move(cut.first);}
            pieces=std::move(next);if(pieces.empty()) break;
        }
        disjoint.insert(disjoint.end(),pieces.begin(),pieces.end());
    }
    double area=0;for(const auto &p:disjoint) area+=projectedArea(p);return area;
}
}
ColumnCellIntersection ColumnVolumeMask::intersect(double x,double y,double size,double h,const std::array<double,3> &ground) const {
    if(h<=0) return {};
    auto surface=ground;surface[2]+=h;
    return intersectSurface({{x,y},{x+size,y},{x+size,y+size},{x,y+size}},surface,ground);
}
ColumnCellIntersection ColumnVolumeMask::intersectSurface(const Footprint &domain,const std::array<double,3> &surface,const std::array<double,3> &ground) const {
    ColumnCellIntersection result;if(bins.empty()||domain.size()<3) return result;
    double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30;for(auto p:domain) {xmin=std::min(xmin,p[0]);xmax=std::max(xmax,p[0]);ymin=std::min(ymin,p[1]);ymax=std::max(ymax,p[1]);}
    std::vector<size_t> ids;
    for(int yy=std::max(y0,int(std::floor(ymin)));yy<=std::min(y0+ny-1,int(std::floor(ymax)));++yy) for(int xx=std::max(x0,int(std::floor(xmin)));xx<=std::min(x0+nx-1,int(std::floor(xmax)));++xx) {const auto &b=bins[size_t(yy-y0)*nx+xx-x0];ids.insert(ids.end(),b.begin(),b.end());}
    std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    std::stable_sort(ids.begin(),ids.end(),[&](size_t a,size_t b){return triangles[a].id<triangles[b].id;});
    if(ids.empty()) return result;
    // 僅保留地面以上且料面以下的實際料柱。 / Material exists only between ground and the supported surface.
    auto positive=splitXY(domain,difference(surface,ground)).first;if(positive.empty()) return result;
    std::vector<MaterialPrism> remaining{{positive,ground,surface}};const double gross=prismVolume(remaining.front());double net=gross;
    for(auto id:ids) {
        const auto &t=triangles[id];HeightPlane lower=ground,upper{0,0,t.top};if(t.relative) {lower[2]+=t.lower;upper=ground;upper[2]+=t.upper;}
        Footprint maskXY=splitXY({t.a,t.b,t.c},difference(upper,lower)).first;if(maskXY.empty()) continue;
        MaterialPrism mask{maskXY,lower,upper};std::vector<MaterialPrism> next;
        for(const auto &piece:remaining) {auto parts=subtractSolid(piece,mask);next.insert(next.end(),parts.begin(),parts.end());}
        if(next.size()>20000) throw std::runtime_error("Structural union local partition exceeds 20000 pieces");
        double after=0;for(const auto &piece:next) after+=prismVolume(piece);
        const double removed=std::max(0.,net-after);
        if(removed>1e-12) {auto found=std::find_if(result.contributions.begin(),result.contributions.end(),[&](const auto &part){return part.first==t.id;});if(found==result.contributions.end()) result.contributions.push_back({t.id,removed});else found->second+=removed;}
        remaining=std::move(next);net=after;if(remaining.empty()) break;
    }
    result.volume=std::max(0.,gross-net);
    if(result.volume>1e-12) result.fullHeightArea=std::max(0.,projectedArea(positive)-residualProjectionArea(remaining));
    return result;
}

bool ColumnVolumeMask::contains(double x,double y,double z) const {
    if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)||bins.empty()||x<x0||y<y0||x>=x0+nx||y>=y0+ny) return false;
    const auto &ids=bins[size_t(int(std::floor(y))-y0)*nx+int(std::floor(x))-x0];
    const std::array<double,2> p{x,y};
    for(auto id:ids) {const auto &t=triangles[id];if(z<=(t.relative?t.reference[0]*x+t.reference[1]*y+t.reference[2]+t.upper:t.top)+.001&&(!t.relative||z>=t.reference[0]*x+t.reference[1]*y+t.reference[2]+t.lower-.001)&&cross(t.a,t.b,p)>=-1e-10&&cross(t.b,t.c,p)>=-1e-10&&cross(t.c,t.a,p)>=-1e-10) return true;}return false;
}

std::vector<int> completeColumnGroup(std::vector<ColumnSemantic> &columns,const std::vector<int> &ids,int axis,const QString &group) {
    if(ids.size()<2||ids.size()>4||group.isEmpty()) throw std::runtime_error("A group requires 2..4 reviewed members and a name");
    auto work=columns;std::vector<size_t> members;std::set<int> unique;int next=0;
    for(const auto &c:work) next=std::max(next,c.id);
    std::vector<std::array<double,2>> centers;double width=0,depth=0;
    /** @brief 步驟 1：只使用指定的確認柱，取最大截面尺寸以避免縮掉已確認柱段。
     *         Step 1: use explicitly confirmed members and maximum section dimensions to avoid shrinking reviewed sections.
     */
    for(auto id:ids) {
        if(!unique.insert(id).second) throw std::runtime_error("Repeated group ID");
        auto it=std::find_if(work.begin(),work.end(),[&](const auto &c){return c.id==id;});
        if(it==work.end()||it->state!=ColumnState::Confirmed||it->shape!="rectangle") throw std::runtime_error("Group members must be confirmed rectangular columns");
        members.push_back(size_t(it-work.begin()));auto b=columnBounds(*it);centers.push_back({(b[0]+b[1])*.5,(b[2]+b[3])*.5});width=std::max(width,b[1]-b[0]);depth=std::max(depth,b[3]-b[2]);
    }
    std::vector<std::array<double,2>> additions;
    /** @brief 步驟 2：兩柱視為相鄰角，以其距離為正方形邊長，沿指定正方向補排。
     *         Step 2: treat two members as adjacent square corners and offset a row by their separation toward the specified positive axis.
     * 保留掃描定位的小幅偏差；位移向量垂直既有邊，並非強制重排既有中心。
     * Preserve small localization errors: the offset is perpendicular to the existing edge, without rearranging observed centers.
     */
    if(ids.size()==2) {
        if(axis!=0&&axis!=1) throw std::runtime_error("Pair completion requires +X or +Y direction");
        const double dx=centers[1][0]-centers[0][0],dy=centers[1][1]-centers[0][1];double vx=-dy,vy=dx;
        if((axis==0?vx:vy)<0) {vx=-vx;vy=-vy;}
        if(std::hypot(dx,dy)<1||std::abs(axis==0?vx:vy)<.9*std::hypot(vx,vy)) throw std::runtime_error("Pair is not an adjacent edge perpendicular to the requested direction");
        for(auto p:centers) additions.push_back({p[0]+vx,p[1]+vy});
    } else {
        /** @brief 步驟 3：三柱由近似直角的共同角及兩鄰角求平行四邊形第四角。
         *         Step 3: infer the fourth corner from the approximately right-angle vertex and its two adjacent corners.
         * 允許矩形跨距；不把既有非精確正方形強制變形成正方形。
         * Allow rectangular spans rather than deforming already observed non-square layouts.
         */
        if(ids.size()==3) {
            size_t corner=0;double best=1e30;
            for(size_t i=0;i<3;++i) {const auto &a=centers[i],&b=centers[(i+1)%3],&c=centers[(i+2)%3];const double bx=b[0]-a[0],by=b[1]-a[1],cx=c[0]-a[0],cy=c[1]-a[1];const double norm=std::hypot(bx,by)*std::hypot(cx,cy);if(norm<1) continue;const double error=std::abs(bx*cx+by*cy)/norm;if(error<best) {best=error;corner=i;}}
            if(best>.2) throw std::runtime_error("Three members do not define a credible right-angle group");
            const auto &a=centers[corner],&b=centers[(corner+1)%3],&c=centers[(corner+2)%3];additions.push_back({b[0]+c[0]-a[0],b[1]+c[1]-a[1]});
        }
        double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30;auto all=centers;all.insert(all.end(),additions.begin(),additions.end());
        for(auto p:all) {xmin=std::min(xmin,p[0]);xmax=std::max(xmax,p[0]);ymin=std::min(ymin,p[1]);ymax=std::max(ymax,p[1]);}
        if(xmax-xmin<1||ymax-ymin<1) throw std::runtime_error("Degenerate four-corner group");
        std::set<int> corners;for(auto p:all) {if(std::min(p[0]-xmin,xmax-p[0])>.2*(xmax-xmin)||std::min(p[1]-ymin,ymax-p[1])>.2*(ymax-ymin)) throw std::runtime_error("Member is not near a group corner");corners.insert((p[0]>(xmin+xmax)*.5?1:0)+(p[1]>(ymin+ymax)*.5?2:0));}
        if(corners.size()!=4) throw std::runtime_error("Duplicate or missing group corners");
    }
    /** @brief 步驟 4：保留既有實測範圍及頂高，統一同組矩形；新增柱採最近鄰柱頂推估。
     *         Step 4: retain evidence/top elevations and equalize rectangles; estimate new tops from nearest members.
     * 新增柱記錄推估來源、零觀測點數，不冒充實測；校驗整場重疊後才提交。
     * Inferred columns retain provenance and zero measured evidence; commit only after site-wide overlap validation.
     */
    std::vector<int> result=ids;
    for(size_t n=0;n<members.size();++n) {auto &c=work[members[n]];if(c.evidenceFootprint.empty()) c.evidenceFootprint=c.footprint;c.footprint=columnSection("rectangle",centers[n][0]-width*.5,centers[n][0]+width*.5,centers[n][1]-depth*.5,centers[n][1]+depth*.5);c.provenance="reviewed_group_"+group;}
    for(auto p:additions) {
        ColumnSemantic c;c.id=++next;c.state=ColumnState::Confirmed;c.provenance="pattern_inferred_group_"+group;
        c.footprint=columnSection("rectangle",p[0]-width*.5,p[0]+width*.5,p[1]-depth*.5,p[1]+depth*.5);
        double best=1e30;size_t nearest=0;for(size_t n=0;n<centers.size();++n) {double d=std::hypot(p[0]-centers[n][0],p[1]-centers[n][1]);if(d<best) {best=d;nearest=n;}}
        c.top=work[members[nearest]].top;c.observedBottom=c.observedTop=c.top;result.push_back(c.id);work.push_back(c);
    }
    const ColumnVolumeMask checked(work);(void)checked;columns=std::move(work);return result;
}

void fitMeasuredColumnGroup(std::vector<ColumnSemantic> &columns,const std::vector<int> &ids,const std::vector<int> &anchors,double ns,double ew,double section,const QString &name) {
    if(ids.size()!=4||anchors.size()<2||!std::isfinite(ns)||!std::isfinite(ew)||!std::isfinite(section)||section<=0||ns<=section||ew<=section||ns>100||ew>100||name.isEmpty()) throw std::runtime_error("Invalid measured group dimensions/anchors");
    auto work=columns;std::vector<size_t> members;std::vector<std::array<double,2>> centers;
    std::set<int> unique,anchorSet(anchors.begin(),anchors.end());
    if(anchorSet.size()!=anchors.size()) throw std::runtime_error("Duplicate measured anchors");
    double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30;
    /** @brief 步驟 1：由既有四個柱中心辨認四角，保留各柱 ID 與角位的對應。
     *         Step 1: identify the four existing corners while retaining the ID-to-corner correspondence.
     */
    for(auto id:ids) {
        if(!unique.insert(id).second) throw std::runtime_error("Duplicate group member");
        auto it=std::find_if(work.begin(),work.end(),[&](const auto &c){return c.id==id;});
        if(it==work.end()||it->state!=ColumnState::Confirmed) throw std::runtime_error("Missing/unconfirmed group member");
        auto b=columnBounds(*it);double x=(b[0]+b[1])*.5,y=(b[2]+b[3])*.5;members.push_back(size_t(it-work.begin()));centers.push_back({x,y});
        xmin=std::min(xmin,x);xmax=std::max(xmax,x);ymin=std::min(ymin,y);ymax=std::max(ymax,y);
    }
    for(auto id:anchors) if(!unique.count(id)) throw std::runtime_error("Anchor outside group");
    if(xmax-xmin<1||ymax-ymin<1) throw std::runtime_error("Degenerate group layout");
    std::set<int> corners;std::vector<std::array<int,2>> signs;
    double originX=0,originY=0;
    /** @brief 步驟 2：固定量測跨距及建築軸方向，以實測錨柱做等權平移最小平方擬合。
     *         Step 2: hold measured spacings/building axes fixed and fit translation by equal-weight least squares over observed anchors.
     * 每個錨柱對組中心的估計是 observed_center − signed_half_spacing，取平均即最小平方解。
     * Each anchor estimates the group center as observed_center minus signed_half_spacing; their mean is the least-squares solution.
     * 原證據截面可避免多次調整累積漂移；推估柱不參與定位。
     * Original evidence sections prevent cumulative drift; inferred members do not anchor placement.
     */
    for(size_t n=0;n<members.size();++n) {
        const int sx=centers[n][0]>(xmin+xmax)*.5?1:-1,sy=centers[n][1]>(ymin+ymax)*.5?1:-1;
        signs.push_back({sx,sy});corners.insert((sx>0?1:0)+(sy>0?2:0));
        const auto &c=work[members[n]];
        if(anchorSet.count(c.id)) {
            if(c.provenance.startsWith("pattern_inferred")) throw std::runtime_error("Inferred column cannot be a measured anchor");
            auto observed=c;if(!c.evidenceFootprint.empty()) observed.footprint=c.evidenceFootprint;auto b=columnBounds(observed);
            originX+=(b[0]+b[1])*.5-sx*ns*.5;originY+=(b[2]+b[3])*.5-sy*ew*.5;
        }
    }
    if(corners.size()!=4) throw std::runtime_error("Ambiguous group corners");
    originX/=anchors.size();originY/=anchors.size();
    /** @brief 步驟 3：依量測跨距重建精確四角，所有柱採指定正方形截面；保留推估來源及柱頂。
     *         Step 3: rebuild exact measured corners with identical square sections, retaining inference provenance and tops.
     */
    for(size_t n=0;n<members.size();++n) {
        auto &c=work[members[n]];if(c.evidenceFootprint.empty()&&!c.provenance.startsWith("pattern_inferred")) c.evidenceFootprint=c.footprint;
        const double x=originX+signs[n][0]*ns*.5,y=originY+signs[n][1]*ew*.5;
        c.footprint=columnSection("rectangle",x-section*.5,x+section*.5,y-section*.5,y+section*.5);c.shape="rectangle";
        c.provenance=(c.provenance.startsWith("pattern_inferred")?"pattern_inferred_group_":"reviewed_group_")+name;
    }
    const ColumnVolumeMask checked(work);(void)checked;columns=std::move(work);
}
