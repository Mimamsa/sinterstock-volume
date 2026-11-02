/** @file columns_test.cpp
 * @brief 語意提議、人工狀態、真實截面及體積排除驗證。 / Proposals, review states, actual sections, and volume subtraction.
 */
#include "cloud.h"
#include <QTemporaryDir>
#include <QCoreApplication>
#include <QProcess>
#include <QFile>
#include <cmath>
#include <iostream>
#include <random>
#include <algorithm>
#include <stdexcept>
void require(bool value,const char *message) {if(!value) throw std::runtime_error(message);}
int main(int argc,char **argv) {
 QCoreApplication app(argc,argv);
 try {
    std::vector<Point> scene;
    for(float center:{2.f,8.f}) for(int iz=0;iz<=40;++iz) for(int side=0;side<4;++side) for(int n=0;n<=6;++n) {
        const float t=-.3f+n*.1f;scene.push_back({center+(side<2?(side==0?-.3f:.3f):t),3+(side>=2?(side==2?-.3f:.3f):t),iz*.2f,42});
    }
    for(int iy=0;iy<=120;++iy) for(int iz=0;iz<=40;++iz) scene.push_back({0,iy*.1f,iz*.2f,42});
    for(int x=0;x<50;++x) for(int y=0;y<50;++y) scene.push_back({x*.1f,y*.1f,1+x*.02f+y*.01f,42});
    ColumnOptions options;options.enabled=true;auto proposals=detectColumnCandidates(scene,options);
    require(proposals.size()==2,"Expected two narrow vertical columns, excluding long wall and pile slope");
    for(const auto &c:proposals) require(c.state==ColumnState::Pending&&c.verticalCoverage>=.65,"Automatic proposal was confirmed");
    auto untouched=scene;require(removeColumnPoints(untouched,proposals)==0,"Pending columns removed measurements");
    ColumnSemantic c;c.id=1;c.state=ColumnState::Confirmed;c.footprint=columnSection("rectangle",1,2,1,2);c.top=4;
    ColumnVolumeMask mask({c});auto partial=mask.intersect(1.95,1.95,.1,3,{0,0,0});
    require(std::abs(partial.volume-.0075)<1e-10&&std::abs(partial.fullHeightArea-.0025)<1e-10,"Partial-cell area intersection failed");
    require(mask.contains(1.5,1.5,3)&&!mask.contains(1.5,1.5,5)&&!mask.contains(3,3,2),"Solid point classification failed");
    c.top=2;ColumnVolumeMask shortMask({c});partial=shortMask.intersect(1,1,.1,3,{0,0,0});
    require(std::abs(partial.volume-.02)<1e-10&&partial.fullHeightArea==0,"Column height clipping failed");
    c.shape="circle";c.footprint=columnSection("circle",0,1,0,1);c.top=4;ColumnVolumeMask circle({c});double v=0;
    for(int y=0;y<10;++y) for(int x=0;x<10;++x) v+=circle.intersect(x*.1,y*.1,.1,2,{0,0,0}).volume;
    require(std::abs(v-3.141592653589793*.25*2)<.003,"Circle deducted its bounding box");
    c.shape="polygon";c.footprint={{0,0},{1,0},{1,.2},{.6,.2},{.6,.8},{1,.8},{1,1},{0,1},{0,.8},{.4,.8},{.4,.2},{0,.2}};
    ColumnVolumeMask beam({c});v=0;
    for(int y=0;y<10;++y) for(int x=0;x<10;++x) v+=beam.intersect(x*.1,y*.1,.1,2,{0,0,0}).volume;
    require(std::abs(v-1.04)<1e-8,"Concave I-section volume failed");
    QTemporaryDir dir;const auto path=dir.filePath("column_semantics.json");saveColumnSemantics(path,"site-test",{c});auto read=loadColumnSemantics(path,"site-test");
    require(read.size()==1&&read[0].state==ColumnState::Confirmed&&read[0].shape=="polygon","Review persistence failed");
    bool bad=false;try {loadColumnSemantics(path,"foreign-site");} catch(const std::runtime_error &) {bad=true;}require(bad,"Foreign site accepted");
    auto duplicate=c;duplicate.id=2;ColumnVolumeMask overlapping({c,duplicate});require(std::abs(overlapping.intersect(0,0,1,2,{0,0,0}).volume-1.04)<1e-8,"Duplicate I-sections were deducted twice");
    // 箱體避讓既有柱：包圍盒交疊但凹截面不重疊，共享邊不可重複扣體積。
    // A box notched around an existing pier may share bounds/edges without sharing solid area.
    auto notched=c;notched.id=3;notched.shape="polygon";notched.top=4;
    notched.footprint={{0,0},{2,0},{2,2},{1.2,2},{1.2,1},{.8,1},{.8,2},{0,2}};
    auto pier=c;pier.id=4;pier.shape="rectangle";pier.top=8;pier.footprint=columnSection("rectangle",.8,1.2,1,2);
    ColumnVolumeMask joined({notched,pier});double joinedVolume=0;
    for(int y=0;y<20;++y) for(int x=0;x<20;++x) joinedVolume+=joined.intersect(x*.1,y*.1,.1,6,{0,0,0}).volume;
    require(std::abs(joinedVolume-16.8)<1e-8,"Notched box/pier overlap was rejected or double-counted");
    require(joined.contains(.5,1.5,3)&&!joined.contains(.5,1.5,5)&&joined.contains(1,1.5,5),"Notched solids lost separate tops");
    auto intrusive=pier;intrusive.footprint=columnSection("rectangle",.79,1.2,1,2);
    ColumnVolumeMask intrusiveUnion({notched,intrusive});
    require(std::abs(intrusiveUnion.intersect(0,0,2,6,{0,0,0}).volume-16.82)<1e-8,"Partial notch overlap was double-counted");
    saveColumnSemantics(dir.filePath("notched.json"),"site-test",{notched,pier});
    require(loadColumnSemantics(dir.filePath("notched.json"),"site-test").size()==2,"Notched semantics did not persist");
    c.provenance="pattern_test";c.evidenceFootprint=columnSection("rectangle",2,3,2,3);
    saveColumnSemantics(path,"site-test",{c});read=loadColumnSemantics(path,"site-test");
    require(read[0].provenance==c.provenance&&read[0].evidenceFootprint==c.evidenceFootprint,"Provenance/evidence lost on save");
    const auto evidence=sampleColumnEvidence({{2.5f,2.5f,1,42}},read);
    require(evidence.size()==1,"Extended section lost original measured evidence");
    c.group="S1";c.legacyIds={53,111};saveColumnSemantics(path,"site-test",{c},true);
    read=loadColumnSemantics(path,"site-test");require(columnSemanticsFinalized(path,"site-test")&&read[0].group=="S1"&&read[0].legacyIds==c.legacyIds,"Finalization/group/legacy IDs lost");
    c.evidenceFootprint.clear();
    c.shape="rectangle";c.footprint=columnSection("rectangle",0,.5,0,.5);c.top=10;
    VolumeOptions volume;volume.enabled=true;volume.automaticGround=false;volume.minimumArea=0;volume.columns={c};
    std::vector<Point> surface;
    for(int y=0;y<10;++y) for(int x=0;x<10;++x) surface.push_back({(x+.5f)*.1f,(y+.5f)*.1f,2,42});
    const auto report=calculateStockpileVolumes(surface,{},{{0,0},{1,0},{1,1},{0,1}},0,volume);
    require(report.error.isEmpty()&&std::abs(report.totalVolume-1.5)<1e-7&&std::abs(report.columnExcludedVolume-.5)<1e-7,"Net pile volume failed");
    require(std::abs(report.grossRasterVolume-report.totalVolume-report.columnExcludedVolume-report.discardedSmallVolume)<1e-8,"Gross/net accounting failed");
    volume.columns[0].state=ColumnState::Rejected;const auto rejected=calculateStockpileVolumes(surface,{},{{0,0},{1,0},{1,1},{0,1}},0,volume);
    require(rejected.error.isEmpty()&&std::abs(rejected.totalVolume-2)<1e-7&&rejected.columnExcludedVolume==0,"Rejected semantics changed volume");
    // 梁底以下不得移除／扣體積，兩層共用 XY 但有不同高度。 / Preserve beam clearance and stack disjoint levels in XY.
    auto lowerBeam=c;lowerBeam.id=801;lowerBeam.kind="beam";lowerBeam.groundRelative=true;lowerBeam.bottomHeight=5;lowerBeam.topHeight=5.9;lowerBeam.referenceGround={0,0,-24};lowerBeam.footprint=columnSection("rectangle",0,1,0,1);
    auto upperBeam=lowerBeam;upperBeam.id=802;upperBeam.bottomHeight=8.1;upperBeam.topHeight=9;
    ColumnVolumeMask levels({lowerBeam,upperBeam});
    require(!levels.contains(.5,.5,-20)&&levels.contains(.5,.5,-18.5)&&!levels.contains(.5,.5,-17),"Elevated beam clearance classification failed");
    auto intersection=levels.intersect(0,0,1,10,{0,0,-24});require(std::abs(intersection.volume-1.8)<1e-8&&intersection.fullHeightArea==0,"Elevated beams deducted the space below");
    intersection=levels.intersect(0,0,1,5.4,{0,0,-24});require(std::abs(intersection.volume-.4)<1e-8,"Partially buried beam failed");
    intersection=levels.intersect(0,0,1,4,{.02,.01,-23});require(intersection.volume==0,"Ground-relative beam did not follow calculation ground");
    saveColumnSemantics(dir.filePath("beams.json"),"site-test",{lowerBeam,upperBeam});
    auto savedBeams=loadColumnSemantics(dir.filePath("beams.json"),"site-test");require(savedBeams[0].groundRelative&&savedBeams[0].bottomHeight==5&&savedBeams[0].topHeight==5.9&&savedBeams[0].kind=="beam","Beam bounds lost on save");
    auto overlapBeam=upperBeam;overlapBeam.bottomHeight=5.8;ColumnVolumeMask beamUnion({lowerBeam,overlapBeam});require(std::abs(beamUnion.intersect(0,0,1,10,{0,0,-24}).volume-4)<1e-8,"Overlapping beam intervals were deducted twice");
    // 柱、箱與梁的三重交集，不能只以整格上限截斷掩飾重複扣除。
    // Triple column/box/beam overlap must use the union, not a full-cell clamp.
    auto post=c;post.id=900;post.footprint=columnSection("rectangle",0,.6,0,1);post.top=2;
    auto box=post;box.id=901;box.footprint=columnSection("rectangle",.4,1,0,1);box.top=3;
    auto crossBeam=lowerBeam;crossBeam.id=902;crossBeam.bottomHeight=1;crossBeam.topHeight=2.5;crossBeam.referenceGround={0,0,0};
    ColumnVolumeMask tripleUnion({post,box,crossBeam}),reverseUnion({crossBeam,box,post});
    auto unionAmount=tripleUnion.intersect(0,0,1,4,{0,0,0});
    require(std::abs(unionAmount.volume-2.8)<1e-8&&unionAmount.fullHeightArea<1e-8,"Triple union volume or projected area failed");
    double contributionSum=0;for(auto part:unionAmount.contributions) contributionSum+=part.second;
    require(std::abs(contributionSum-unionAmount.volume)<1e-8&&reverseUnion.intersect(0,0,1,4,{0,0,0}).contributions==unionAmount.contributions,"Union attribution changed with input order");
    unionAmount=tripleUnion.intersect(0,0,1,2.7,{0,0,0});
    require(std::abs(unionAmount.volume-2.62)<1e-8&&std::abs(unionAmount.fullHeightArea-.6)<1e-8,"Partially buried union/stacked full-height area failed");
    unionAmount=tripleUnion.intersect(0,0,1,2.2,{0,0,0});
    require(std::abs(unionAmount.volume-2.2)<1e-8&&std::abs(unionAmount.fullHeightArea-1)<1e-8,"Adjacent vertical intervals did not fully exclude material");
    // 與獨立座標分割解析解比較：隨機長方體、部分出界及多重重疊。
    // Compare with an independent coordinate-partition oracle for overlapping and partially exterior boxes.
    std::mt19937 unionRandom(721);std::uniform_int_distribution<int> digit(0,10);
    for(int trial=0;trial<24;++trial) {
        std::vector<ColumnSemantic> solids;std::vector<double> xs{0,1},ys{0,1},zs{0,2.7};
        for(int i=0;i<6;++i) {
            auto solid=lowerBeam;solid.id=1000+i;solid.referenceGround={0,0,0};
            const double x=digit(unionRandom)*.1-.2,y=digit(unionRandom)*.1-.2,w=.2+digit(unionRandom)*.1;
            solid.footprint=columnSection("rectangle",x,x+w,y,y+w);solid.bottomHeight=digit(unionRandom)*.2;solid.topHeight=solid.bottomHeight+.2+digit(unionRandom)*.2;
            solids.push_back(solid);for(double v:{x,x+w}) if(v>0&&v<1) xs.push_back(v);for(double v:{y,y+w}) if(v>0&&v<1) ys.push_back(v);for(double v:{solid.bottomHeight,solid.topHeight}) if(v>0&&v<2.7) zs.push_back(v);
        }
        for(auto *v:{&xs,&ys,&zs}) {std::sort(v->begin(),v->end());v->erase(std::unique(v->begin(),v->end()),v->end());}
        double expectedVolume=0,expectedFullArea=0;
        for(size_t x=1;x<xs.size();++x) for(size_t y=1;y<ys.size();++y) {
            const double cx=(xs[x]+xs[x-1])*.5,cy=(ys[y]+ys[y-1])*.5,area=(xs[x]-xs[x-1])*(ys[y]-ys[y-1]);bool full=true;
            for(size_t z=1;z<zs.size();++z) {
                const double cz=(zs[z]+zs[z-1])*.5;bool covered=false;
                for(const auto &solid:solids) {const auto b=columnBounds(solid);if(cx>b[0]&&cx<b[1]&&cy>b[2]&&cy<b[3]&&cz>solid.bottomHeight&&cz<solid.topHeight) covered=true;}
                if(covered) expectedVolume+=area*(zs[z]-zs[z-1]);else full=false;
            }
            if(full) expectedFullArea+=area;
        }
        auto actual=ColumnVolumeMask(solids).intersect(0,0,1,2.7,{0,0,0});
        require(std::abs(actual.volume-expectedVolume)<1e-8&&std::abs(actual.fullHeightArea-expectedFullArea)<1e-8,"Structural union disagrees with independent partition oracle");
    }
    auto belowFloor=c;belowFloor.id=1100;belowFloor.footprint=columnSection("rectangle",0,1,0,1);belowFloor.top=0;
    require(std::abs(ColumnVolumeMask({belowFloor}).intersect(0,0,1,3,{2,0,-1}).volume-.25)<1e-8,"Sloped floor crossing absolute structure top failed");
    // 四角補齊：三角補第四角、兩柱向指定方向補排、同組尺寸及不移動既有中心。
    // Four-corner completion: three-to-four inference, directed pair completion, equal sections, stable observed centers.
    auto makeMember=[](int id,double x,double y,double w) {ColumnSemantic p;p.id=id;p.state=ColumnState::Confirmed;p.top=5;p.footprint=columnSection("rectangle",x-w*.5,x+w*.5,y-.25,y+.25);return p;};
    auto triple=std::vector<ColumnSemantic>{makeMember(1,2,2,.5),makeMember(2,8,2,1),makeMember(3,8,8,.5)};
    auto groupIds=completeColumnGroup(triple,{1,2,3},-1,"G-test");require(groupIds.size()==4&&triple.size()==4,"Triple completion count failed");
    auto corner=columnBounds(triple.back());require(std::abs((corner[0]+corner[1])*.5-2)<1e-8&&std::abs((corner[2]+corner[3])*.5-8)<1e-8,"Wrong missing corner");
    require(triple.back().evidencePoints==0&&triple.back().provenance.startsWith("pattern_inferred_group"),"Inference presented as observation");
    for(const auto &p:triple) {auto b=columnBounds(p);require(std::abs(b[1]-b[0]-1)<1e-8,"Group section mismatch");}
    auto pair=std::vector<ColumnSemantic>{makeMember(1,2,2,.5),makeMember(2,8,2,.5)};
    completeColumnGroup(pair,{1,2},1,"west");require(pair.size()==4&&columnBounds(pair[2])[2]>7,"+Y pair direction failed");
    auto north=std::vector<ColumnSemantic>{makeMember(1,2,2,.5),makeMember(2,2,8,.5)};
    completeColumnGroup(north,{1,2},0,"north");require(columnBounds(north[2])[0]>7,"+X pair direction failed");
    const auto beforeFailure=north.size();bad=false;try {completeColumnGroup(north,{1,1},0,"invalid");} catch(const std::runtime_error &) {bad=true;}
    require(bad&&north.size()==beforeFailure,"Failed group operation changed model");
    // 量測長方形：不受推估柱位置影響，精確套用 6.4 × 5.0 m 中心距及 0.8 m 截面。
    // Measured rectangle: fit only observed anchors, enforcing exact 6.4 x 5.0 m spacings and 0.8 m sections.
    fitMeasuredColumnGroup(triple,{1,2,3,4},{1,2,3},6.4,5.,.8,"S1");
    std::vector<std::array<double,4>> measuredBounds;for(const auto &p:triple) {auto b=columnBounds(p);measuredBounds.push_back(b);require(std::abs(b[1]-b[0]-.8)<1e-9&&std::abs(b[3]-b[2]-.8)<1e-9,"Measured section mismatch");}
    require(std::abs(measuredBounds[1][0]-measuredBounds[0][0]-6.4)<1e-9&&std::abs(measuredBounds[3][2]-measuredBounds[0][2]-5)<1e-9,"Measured center spacing mismatch");
    require(triple.back().provenance=="pattern_inferred_group_S1","Lost inference provenance");
    fitMeasuredColumnGroup(triple,{1,2,3,4},{1,2,3},6.4,5.,.8,"S1");
    for(size_t n=0;n<triple.size();++n) require(columnBounds(triple[n])==measuredBounds[n],"Repeated fit drifted");
    bad=false;try {fitMeasuredColumnGroup(triple,{1,2,3,4},{1,4},6.4,5,.8,"bad");} catch(const std::runtime_error &) {bad=true;}require(bad,"Inferred anchor accepted");
    if(argc>1) {
        // 規律補齊整合驗證：只補主列內部缺柱，尊重否決，且保持兩壁同尺寸。
        // Completion integration: interior bays only, honor rejection, and equal wall dimensions.
        std::vector<ColumnSemantic> reviewed;int id=0;
        for(double x:{1.5,20.,52.}) for(int bay=0;bay<12;++bay) {
            ColumnSemantic p;p.id=++id;p.state=ColumnState::Confirmed;p.top=5;p.footprint=columnSection("rectangle",x-.4,x+.4,10+bay*12-.3,10+bay*12+.3);
            if(x==1.5&&bay==3) continue;
            if(x==20.&&bay==5) continue;
            if(x==20.&&bay==7) p.state=ColumnState::Rejected;
            reviewed.push_back(p);
        }
        const auto input=dir.filePath("reviewed.json"),region=dir.filePath("region.json"),output=dir.filePath("completed.json");
        saveColumnSemantics(input,"site-test",reviewed);QFile file(region);require(file.open(QIODevice::WriteOnly),"Region write failed");
        file.write(R"({"site_id":"site-test","footprint_xy_m":[[0,0],[54,0],[54,160],[0,160]]})");file.close();
        QProcess process;process.start(QString::fromLocal8Bit(argv[1]),{input,region,output});require(process.waitForFinished(30000)&&process.exitCode()==0,"Pattern tool failed");
        const auto completed=loadColumnSemantics(output,"site-test");size_t inferred=0;double sectionDepth=-1;
        for(const auto &p:completed) {
            const auto b=columnBounds(p);
            if(p.provenance.startsWith("pattern_inferred")) {++inferred;require(p.evidencePoints==0,"Inference fabricated measurements");require(std::abs((b[2]+b[3])*.5-94)>1,"Overrode explicitly rejected bay");}
            if(p.provenance=="reviewed_wall_extended") {
                require(std::abs(b[0])<1e-8||std::abs(b[1]-54)<1e-8,"Column not attached to wall");
                if(sectionDepth<0) sectionDepth=b[1]-b[0];else require(std::abs(sectionDepth-(b[1]-b[0]))<1e-8,"Wall column sizes differ");
            }
        }
        require(inferred==2,"Unexpected inferred bay count");
        process.start(QString::fromLocal8Bit(argv[1]),{input,region,output});require(process.waitForFinished(30000)&&process.exitCode()!=0,"Tool overwrote existing result");
    }
    std::cout<<"PASS: proposals, no auto acceptance, partial raster intersections, finite height, round/I sections, persistence, buried structural unions, gross/net volumes\n";
 } catch(const std::exception &e) {std::cerr<<e.what()<<"\n";return 1;}
}
