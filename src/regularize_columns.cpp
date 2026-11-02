/** @file regularize_columns.cpp
 * @brief 依人工確認的跨距補柱並將等尺寸壁柱延伸到牆面。 / Complete reviewed structural bays and extend equal-sized wall piers to walls.
 * 只處理由審查資料支持的規律柱列，不推廣到設備柱。 / Operate only on regular reviewed rows, never irregular equipment supports.
 */
#include "columns.h"
#include <QCoreApplication>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace {
double median(std::vector<double> a) {std::sort(a.begin(),a.end());return a[a.size()/2];}
double cx(const ColumnSemantic &c) {auto b=columnBounds(c);return (b[0]+b[1])*.5;}
double cy(const ColumnSemantic &c) {auto b=columnBounds(c);return (b[2]+b[3])*.5;}
struct Row {double x=0;std::vector<size_t> members;std::vector<std::vector<size_t>> bays;};
}
/** @brief 輸入審核檔、場域邊界與輸出路徑；輸出需與輸入分開，保存原始審查。
 *         Read reviewed semantics, site boundary and a distinct output path, preserving the original review.
 * @return 0 成功；1 為資料或幾何驗證失敗。 / Zero on success, one on invalid data/geometry.
 */
int main(int argc,char **argv) {
 QCoreApplication app(argc,argv);QTextStream out(stdout),err(stderr);
 try {
    const auto args=app.arguments();
    if(args.size()!=4) throw std::runtime_error("Usage: regularize_columns reviewed.json processing_region.json output.json");
    if(args[1]==args[3]||QFile::exists(args[3])) throw std::runtime_error("Output must be new; preserve your original review");
    QFile f(args[2]);if(!f.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot read region");
    const auto region=QJsonDocument::fromJson(f.readAll()).object();const auto site=region.value("site_id").toString();
    if(site.isEmpty()) throw std::runtime_error("Missing site ID");
    Footprint footprint;double left=1e30,right=-1e30;
    for(auto v:region.value("footprint_xy_m").toArray()) {auto p=v.toArray();if(p.size()!=2) throw std::runtime_error("Invalid boundary");footprint.push_back({p[0].toDouble(),p[1].toDouble()});left=std::min(left,p[0].toDouble());right=std::max(right,p[0].toDouble());}
    validateFootprint(footprint);
    auto columns=loadColumnSemantics(args[1],site);const auto original=columns;int next=0;
    std::vector<size_t> confirmed;for(size_t i=0;i<columns.size();++i) {next=std::max(next,columns[i].id);if(columns[i].state==ColumnState::Confirmed) confirmed.push_back(i);}
    std::sort(confirmed.begin(),confirmed.end(),[&](auto a,auto b){return cx(columns[a])<cx(columns[b]);});
    /** @brief 步驟 1：以橫向位置分柱列，再合併同一 Y 跨的多個量測碎片。
     *         Step 1: group by transverse position, then merge multiple measured fragments in the same longitudinal bay.
     * 2 m 橫向容差、1 m 同跨容差只用於歸組，不直接當作柱寬。
     * The 2 m transverse and 1 m bay tolerances group observations, not physical dimensions.
     */
    std::vector<Row> rows;
    for(auto i:confirmed) {auto it=std::find_if(rows.begin(),rows.end(),[&](auto &r){return std::abs(cx(columns[i])-r.x)<=2.;});
        if(it==rows.end()) {rows.push_back({cx(columns[i]),{i},{}});} else {it->members.push_back(i);std::vector<double> xs;for(auto j:it->members) xs.push_back(cx(columns[j]));it->x=median(xs);}}
    for(auto &r:rows) {
        std::sort(r.members.begin(),r.members.end(),[&](auto a,auto b){return cy(columns[a])<cy(columns[b]);});
        for(auto i:r.members) {if(r.bays.empty()||cy(columns[i])-cy(columns[r.bays.back().front()])>1.) r.bays.push_back({i});else r.bays.back().push_back(i);}
    }
    /** @brief 步驟 2：至少十跨、六成相鄰跨距為 10–16 m 才認為是主結構柱列。
     *         Step 2: retain rows with ten or more bays and at least 60% of gaps in 10–16 m.
     * 此範圍由本場域已審查的 12/15 m 跨距建立，不適用其他建築。
     * This site-specific range follows reviewed 12/15 m bays and is not a universal building prior.
     */
    rows.erase(std::remove_if(rows.begin(),rows.end(),[&](auto &r){size_t good=0;for(size_t j=1;j<r.bays.size();++j) {double d=cy(columns[r.bays[j].front()])-cy(columns[r.bays[j-1].front()]);if(d>=10&&d<=16) ++good;}return r.bays.size()<10||double(good)/double(r.bays.size()-1)<.6;}),rows.end());
    if(rows.size()<3) throw std::runtime_error("Insufficient regular reviewed rows; no extrapolation performed");
    struct Station {std::vector<double> ys;std::set<size_t> rows;};std::vector<std::pair<double,size_t>> centers;
    for(size_t r=0;r<rows.size();++r) for(auto &bay:rows[r].bays) {std::vector<double> ys;for(auto i:bay) ys.push_back(cy(columns[i]));centers.push_back({median(ys),r});}
    std::sort(centers.begin(),centers.end());std::vector<Station> stations;
    for(auto p:centers) {if(stations.empty()||p.first-median(stations.back().ys)>1.) stations.push_back({{p.first},{p.second}});else {stations.back().ys.push_back(p.first);stations.back().rows.insert(p.second);}}
    /** @brief 步驟 3：左右壁柱使用同一實心矩形，沿牆寬取中位數，向內深度取最大已確認伸入量。
     *         Step 3: use one solid rectangle for wall piers: median along-wall width and maximum reviewed inward reach.
     * 保留所有已確認柱段的內側位置，以延伸到牆方式貼牆，並非把量測柱段移走。
     * Preserve every reviewed inward face and extend to the wall, rather than translating measured columns.
     */
    std::vector<double> widths;double depth=0;
    for(const auto &r:rows) if(r.x-left<4||right-r.x<4) for(auto i:r.members) {auto b=columnBounds(columns[i]);widths.push_back(b[3]-b[2]);depth=std::max(depth,r.x-left<4?b[1]-left:right-b[0]);}
    if(widths.empty()) throw std::runtime_error("No wall rows");
    const double width=median(widths);
    int merged=0,normalized=0,added=0;
    for(auto &r:rows) {
        const bool wall=r.x-left<4||right-r.x<4;const bool west=r.x-left<4;
        std::vector<size_t> representatives;
        for(auto bay:r.bays) {
            auto master=*std::max_element(bay.begin(),bay.end(),[&](auto a,auto b){auto aa=columnBounds(columns[a]),bb=columnBounds(columns[b]);return (aa[1]-aa[0])*(aa[3]-aa[2])<(bb[1]-bb[0])*(bb[3]-bb[2]);});
            representatives.push_back(master);auto &c=columns[master];auto bounds=columnBounds(c);
            c.evidenceFootprint=c.footprint;
            for(auto i:bay) if(i!=master) {auto b=columnBounds(columns[i]);bounds={std::min(bounds[0],b[0]),std::max(bounds[1],b[1]),std::min(bounds[2],b[2]),std::max(bounds[3],b[3])};c.top=std::max(c.top,columns[i].top);columns[i].state=ColumnState::Rejected;columns[i].provenance=QString("merged_into_C%1").arg(c.id);++merged;}
            c.evidenceFootprint=columnSection("rectangle",bounds[0],bounds[1],bounds[2],bounds[3]);
            if(wall) {const double y=cy(original[master]);c.footprint=columnSection("rectangle",west?left:right-depth,west?left+depth:right,y-width*.5,y+width*.5);c.shape="rectangle";c.provenance="reviewed_wall_extended";++normalized;}
            else if(bay.size()>1) {c.footprint=c.evidenceFootprint;c.provenance="reviewed_fragments_merged";}
        }
        std::vector<double> rowY;for(auto i:representatives) rowY.push_back(cy(columns[i]));
        const double low=*std::min_element(rowY.begin(),rowY.end()),high=*std::max_element(rowY.begin(),rowY.end());
        /** @brief 步驟 4：僅在已有柱列的兩端之間、且至少兩條其他柱列支持的跨位補缺柱。
         *         Step 4: fill only interior gaps supported by at least two other rows; never extrapolate row ends.
         * 拒絕項在同一跨位則不補。新增柱保留推估來源、零量測數，避免誤認為掃描實測。
         * Explicitly rejected stations block completion. Inferred columns retain provenance and zero observed evidence.
         */
        for(const auto &s:stations) {
            const double y=median(s.ys);if(s.rows.size()<2||y<=low||y>=high) continue;
            if(std::any_of(representatives.begin(),representatives.end(),[&](auto i){return std::abs(cy(columns[i])-y)<1.;})) continue;
            bool rejected=false;for(const auto &c:original) if(c.state==ColumnState::Rejected&&std::abs(cx(c)-r.x)<1.&&std::abs(cy(c)-y)<1.) rejected=true;
            if(rejected) {out<<"Skipped explicitly rejected row "<<r.x<<" Y "<<y<<"\n";continue;}
            std::vector<double> dx,dy;size_t lower=representatives.front(),upper=representatives.back();
            for(auto i:representatives) {auto b=columnBounds(columns[i]);dx.push_back(b[1]-b[0]);dy.push_back(b[3]-b[2]);if(cy(columns[i])<y&&cy(columns[i])>cy(columns[lower])) lower=i;if(cy(columns[i])>y&&cy(columns[i])<cy(columns[upper])) upper=i;}
            const double t=(y-cy(columns[lower]))/(cy(columns[upper])-cy(columns[lower]));
            ColumnSemantic c;c.id=++next;c.state=ColumnState::Confirmed;c.top=columns[lower].top*(1-t)+columns[upper].top*t;c.observedBottom=c.observedTop=c.top;
            c.provenance=QString("pattern_inferred_between_C%1_C%2").arg(columns[lower].id).arg(columns[upper].id);
            const double d=wall?depth:median(dx),w=wall?width:median(dy),x=wall?(west?left+d*.5:right-d*.5):r.x;
            c.footprint=columnSection("rectangle",x-d*.5,x+d*.5,y-w*.5,y+w*.5);
            // 新增幾何須在真實場域內；不跨越倒凸字凹角。 / New geometry must remain within the actual concave site.
            if(!std::all_of(c.footprint.begin(),c.footprint.end(),[&](const auto &p){return insideFootprint(footprint,p[0],p[1]);})) throw std::runtime_error("Inferred column outside site");
            columns.push_back(c);++added;out<<"Added C"<<c.id<<" row X "<<r.x<<" Y "<<y<<" top "<<c.top<<"\n";
        }
    }
    saveColumnSemantics(args[3],site,columns,columnSemanticsFinalized(args[1],site));
    out<<"Regular rows: "<<rows.size()<<"; wall size "<<depth<<" x "<<width<<" m; normalized "<<normalized<<"; merged fragments "<<merged<<"; added "<<added<<"\n";
    return 0;
 } catch(const std::exception &e) {err<<e.what()<<"\n";return 1;}
}
