/** @file complete_column_groups.cpp
 * @brief 依人工柱組補角或套用量測長方形，另存模型。 / Complete reviewed corners or fit measured rectangles into a new model.
 */
#include "columns.h"
#include <QCoreApplication>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QTextStream>
#include <set>
#include <algorithm>
#include <cmath>
#include <stdexcept>
int main(int argc,char **argv) {
 QCoreApplication app(argc,argv);QTextStream out(stdout),err(stderr);
 try {
    const auto a=app.arguments();if(a.size()!=4) throw std::runtime_error("Usage: complete_column_groups input.json groups.json NEW-output.json");
    if(QFile::exists(a[3])) throw std::runtime_error("Output exists; refusing overwrite");
    QFile f(a[2]);if(!f.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot open group definitions");
    QJsonParseError error;const auto root=QJsonDocument::fromJson(f.readAll(),&error).object();
    if(error.error!=QJsonParseError::NoError||root.value("site_id").toString().isEmpty()||!root.value("groups").isArray()) throw std::runtime_error("Invalid group definitions");
    const auto site=root.value("site_id").toString();auto columns=loadColumnSemantics(a[1],site);std::set<int> used;
    for(auto value:root.value("groups").toArray()) {const auto g=value.toObject();std::vector<int> ids;
        for(auto id:g.value("member_ids").toArray()) {if(!id.isDouble()||id.toInt()<=0||!used.insert(id.toInt()).second) throw std::runtime_error("Invalid or overlapping group membership");ids.push_back(id.toInt());}
        const auto direction=g.value("missing_direction").toString();int axis=direction=="+X"?0:direction=="+Y"?1:-1;
        std::vector<int> members;
        if(g.contains("north_south_m")) {
            std::vector<int> anchors;for(auto id:g.value("anchor_ids").toArray()) anchors.push_back(id.toInt());
            fitMeasuredColumnGroup(columns,ids,anchors,g.value("north_south_m").toDouble(),g.value("east_west_m").toDouble(),g.value("section_m").toDouble(),g.value("name").toString());members=ids;
        } else members=completeColumnGroup(columns,ids,axis,g.value("name").toString());out<<g.value("name").toString()<<":";for(auto id:members) out<<" C"<<id;out<<"\n";
    }
    /** @brief 依已審查標準柱統一同列截面與 X 中心，保留各柱 Y、柱頂及原證據。
     *         Standardize row sections and X alignment from reviewed references while retaining Y, tops, and evidence.
     * 標準柱尺寸及 X 必須一致，避免將錯誤參考平均成新幾何；所有修改仍先在記憶體驗證再另存。
     * Reference dimensions/X must agree rather than averaging incompatible references; validate before saving a new model.
     */
    for(auto value:root.value("rows").toArray()) {
        const auto row=value.toObject();const auto name=row.value("name").toString();
        if(name.isEmpty()||row.value("reference_ids").toArray().size()<2) throw std::runtime_error("Row requires a name and at least two references");
        std::vector<int> targets,references;std::set<int> rowIds,refIds;
        for(auto id:row.value("member_ids").toArray()) {if(id.toInt()<=0||!rowIds.insert(id.toInt()).second||!used.insert(id.toInt()).second) throw std::runtime_error("Invalid/overlapping row membership");targets.push_back(id.toInt());}
        for(auto id:row.value("reference_ids").toArray()) {if(!rowIds.count(id.toInt())||!refIds.insert(id.toInt()).second) throw std::runtime_error("Invalid row reference");references.push_back(id.toInt());}
        double x=0,w=0,d=0;
        for(auto id:references) {auto it=std::find_if(columns.begin(),columns.end(),[&](const auto &c){return c.id==id;});
            if(it==columns.end()||it->state!=ColumnState::Confirmed||it->shape!="rectangle") throw std::runtime_error("Invalid reference column");
            const auto b=columnBounds(*it);const double xx=(b[0]+b[1])*.5,ww=b[1]-b[0],dd=b[3]-b[2];
            if(id==references.front()) {x=xx;w=ww;d=dd;} else if(std::abs(x-xx)>1e-6||std::abs(w-ww)>1e-6||std::abs(d-dd)>1e-6) throw std::runtime_error("Row references disagree in section or alignment");
        }
        for(auto id:targets) {auto it=std::find_if(columns.begin(),columns.end(),[&](const auto &c){return c.id==id;});
            if(it==columns.end()||it->state!=ColumnState::Confirmed||it->shape!="rectangle") throw std::runtime_error("Invalid row column");
            const auto b=columnBounds(*it);const double y=(b[2]+b[3])*.5;
            if(it->evidenceFootprint.empty()&&!it->provenance.startsWith("pattern_inferred")) it->evidenceFootprint=it->footprint;
            it->footprint=columnSection("rectangle",x-w*.5,x+w*.5,y-d*.5,y+d*.5);
            const auto tag="_row_"+name;if(!it->provenance.contains(tag)) it->provenance=(it->provenance.isEmpty()?QString("reviewed"):it->provenance)+tag;
        }
        out<<name<<": "<<targets.size()<<" pillars, X="<<x<<", section="<<w<<" x "<<d<<" m\n";
    }
    saveColumnSemantics(a[3],site,columns,columnSemanticsFinalized(a[1],site));return 0;
 } catch(const std::exception &e) {err<<e.what()<<"\n";return 1;}
}
