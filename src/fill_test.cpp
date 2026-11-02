/**
 * @file fill_test.cpp
 * @brief 合成堆面補洞與保守拒絕條件測試。 / Synthetic pile hole-filling and conservative rejection tests.
 */
#include "cloud.h"
#include <QFile>
#include <QDir>
#include <QtEndian>
#include <cstring>
#include <QTemporaryDir>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

/** @brief 建立已知斜平面與中央 3x3 缺口。 / Build a known sloping plane with a central 3x3 hole.
 * @param hole 是否省略中央量測點。 / Whether to omit central measurements.
 * @param slope X 向坡度。 / X-direction slope.
 * @return 量測點集。 / Measured point set.
 */
std::vector<Point> plane(bool hole=true, float slope=.2f) {
    std::vector<Point> points;
    for (int y=0;y<9;++y) for (int x=0;x<9;++x) {
        if (hole && x>=3 && x<=5 && y>=3 && y<=5) continue;
        float xx=float((x+.5)*.2),yy=float((y+.5)*.2);
        points.push_back({xx,yy,1+slope*xx+.3f*yy,42});
    }
    return points;
}

/** @brief 驗證內插準確性、禁止外插、斷層拒絕與顯示預算。 / Verify interpolation, no extrapolation, discontinuity rejection, and display budgets.
 * @return 全部通過為 0，失敗為 1。 / Zero on success, one on failure.
 */
int main() {
    FillOptions o{true,0,1.8,0,1.8,0,10,.2,1,.05,2};
    size_t holes=0;
    const auto measured=plane();
    const auto result=fillStockpileHoles(measured,o,holes);
    if (result.size()!=9 || holes!=1 || measured.size()!=72) return 1;
    for (const auto &p : result) {
        if (std::abs(p.z-(1+.2*p.x+.3*p.y))>1e-5 || !std::isnan(p.intensity) ||
            p.x<.6 || p.x>1.2 || p.y<.6 || p.y>1.2) return 1;
    }
    // 量測 XY 不在格心時，仍應重建相同斜面。 / Reconstruct the same slope when measured XY is off-center.
    auto offset=measured;
    for (auto &p : offset) {p.x+=.025f; p.y-=.03f; p.z=1+.2f*p.x+.3f*p.y;}
    auto offsetFilled=fillStockpileHoles(offset,o,holes);
    if (offsetFilled.size()!=9) return 1;
    for (const auto &p : offsetFilled) if (std::abs(p.z-(1+.2*p.x+.3*p.y))>1e-5) return 1;
    if (!fillStockpileHoles(plane(false),o,holes).empty()) return 1;
    auto small=o; small.maxHole=.4;
    if (!fillStockpileHoles(measured,small,holes).empty()) return 1;
    auto open=measured;
    open.erase(std::remove_if(open.begin(),open.end(),[](Point p) {return p.x<.6 && p.y>.8 && p.y<1.;}),open.end());
    if (!fillStockpileHoles(open,o,holes).empty()) return 1;
    auto layers=measured;
    for (auto p : measured) {p.z+=3; layers.push_back(p);}
    if (!fillStockpileHoles(layers,o,holes).empty()) return 1;
    if (!fillStockpileHoles(plane(true,3),o,holes).empty()) return 1;
    auto step=measured;
    for (auto &p : step) if (p.x>1) p.z+=1;
    if (!fillStockpileHoles(step,o,holes).empty()) return 1;
    auto outside=o; outside.zMin=8;
    if (!fillStockpileHoles(measured,outside,holes).empty()) return 1;
    auto disabled=o; disabled.enabled=false;
    if (!fillStockpileHoles(measured,disabled,holes).empty()) return 1;
    auto invalid=o; invalid.xMax=100000;
    bool rejected=false;
    try { fillStockpileHoles(measured,invalid,holes); } catch (const std::runtime_error &) { rejected=true; }
    if (!rejected) return 1;
    // 使用完整資料補洞，合成／量測顯示樣本共用一個上限。 / Fill from full data; synthetic/measured samples share one display cap.
    QTemporaryDir dir; if (!dir.isValid()) return 1;
    const auto path=dir.filePath("pile.ply"); QFile file(path);
    if (!file.open(QIODevice::WriteOnly)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 72\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n");
    for (const auto &p : measured)
        file.write(QByteArray::number(p.x,'g',9)+" "+QByteArray::number(p.y,'g',9)+" "+QByteArray::number(p.z,'g',9)+" 42\n");
    file.close();
    const auto cloud=processCloud(path,1,{},o);
    if (!cloud.error.isEmpty() || cloud.filledCount!=9 || cloud.filledHoles!=1 || cloud.reducedCount!=72 ||
        cloud.after.size()+cloud.filled.size()!=1 || cloud.before.size()!=1) return 1;
    // 匯出在顯示抽樣前執行：1 點預覽仍須保存 72 量測 + 9 補點。
    // Export precedes display sampling: a one-point preview must still save 72 measurements and 9 estimates.
    const auto output=dir.filePath("filled.ply");
    const auto exported=processCloud(path,1,{},o,true,{}, {},output);
    if(!exported.error.isEmpty() || !exported.exportError.isEmpty() || exported.exportedPath!=output) return 1;
    QFile saved(output);if(!saved.open(QIODevice::ReadOnly)) return 1;
    const auto bytes=saved.readAll();saved.close();const auto end=bytes.indexOf("end_header\n");
    if(end<0 || !bytes.left(end).contains("element vertex 81\n") || !bytes.left(end).contains("property uchar is_estimated")) return 1;
    const auto payload=bytes.mid(end+11);if(payload.size()!=81*17) return 1;
    for(int i=0;i<81;++i) {
        const auto *record=reinterpret_cast<const uchar *>(payload.constData()+i*17);
        const quint32 bits=qFromLittleEndian<quint32>(record+12);float intensity;std::memcpy(&intensity,&bits,4);
        if(record[16]!=(i>=72?1:0) || (i>=72?!std::isnan(intensity):intensity!=42)) return 1;
    }
    const auto reread=processCloud(output,1,{}, {},false);
    if(!reread.error.isEmpty() || reread.originalCount!=81) return 1;
    const auto duplicate=processCloud(path,1,{},o,true,{}, {},output);
    if(duplicate.exportError.isEmpty() || !duplicate.error.isEmpty()) return 1;
    if(!saved.open(QIODevice::ReadOnly) || saved.readAll()!=bytes) return 1;saved.close();
    if(!QDir(dir.path()).entryList({".pointcloud-export-*"},QDir::Files|QDir::Hidden).isEmpty()) return 1;
    const auto calibratedPath=dir.filePath("calibrated_export.ply");
    exportProcessedPly(calibratedPath,{{2,3,-4,42}},{{4,5,-3}},true);
    const auto calibratedReload=processCloud(calibratedPath,0);
    if(!calibratedReload.error.isEmpty() || calibratedReload.coordinateFrame.isEmpty() || calibratedReload.originalCount!=2 ||
       calibratedReload.before[0].x!=2 || calibratedReload.before[0].z!=-4) return 1;
    bool badParent=false;try {exportProcessedPly(dir.filePath("missing/failed.ply"),{}, {},false);} catch(const std::runtime_error &) {badParent=true;}
    if(!badParent) return 1;
    auto global=o; global.fullDomain=true;
    global.xMin=global.xMax=global.yMin=global.yMax=global.zMin=global.zMax=0;
    auto fullDomain=processCloud(path,100,{},global);
    if (!fullDomain.error.isEmpty() || !fullDomain.filledCount || !fullDomain.appliedFill.fullDomain ||
        fullDomain.appliedFill.xMin>=fullDomain.cleanLow.x || fullDomain.appliedFill.xMax<=fullDomain.cleanHigh.x ||
        fullDomain.filled.size()!=fullDomain.filledCount || fullDomain.after.size()+fullDomain.filled.size()>100) return 1;
    // 流程回歸：先裁切牆邊再補洞；開放邊界不得被當作封閉缺口補回。
    // Pipeline regression: crop walls before filling; an opened boundary must not be filled as an enclosed gap.
    QFile tagged(path); if(!tagged.open(QIODevice::ReadOnly)) return 1;
    auto data=tagged.readAll();tagged.close();
    data.replace("format ascii 1.0\n","format ascii 1.0\ncomment calibrated building frame; meters; origin at reference roof corner\n");
    if(!tagged.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    tagged.write(data);tagged.close();
    QFile model(dir.filePath("site_model.json"));if(!model.open(QIODevice::WriteOnly)) return 1;
    model.write(R"({"id":"synthetic-site","planes":[[1,0,0,0],[1,0,0,-10],[0,1,0,0],[0,1,0,-20],[0,0,1,-10],[0,0,1,-10]]})");model.close();
    StructureOptions filtering;filtering.surfaceFilter.enabled=false;filtering.enabled=true;filtering.wallMargin=.9;filtering.roofMargin=.5;filtering.cellSize=.2;
    const auto filtered=processCloud(path,0,{},o,true,filtering);
    if(!filtered.error.isEmpty() || filtered.filledCount!=0 || filtered.structureRemoved==0 || !filtered.floatingFiltered ||
       filtered.after.size()!=filtered.reducedCount || filtered.filled.size()!=filtered.filledCount) {
        std::cerr<<"Structure-before-fill regression: "<<filtered.error.toStdString()<<"\n";return 1;
    }
    for(const auto &p:filtered.filled) if(!std::isnan(p.intensity)) return 1;
    for(const auto &p:filtered.after) if(p.intensity!=42) return 1;
    // 整合驗證：與料堆 XY 不重疊的密集懸浮群，必須在補洞之前被移除。
    // Integration: a dense airborne cluster with no underlying pile observations must be removed before filling.
    data.replace("element vertex 72", "element vertex 97");
    for(int x=0;x<5;++x) for(int y=0;y<5;++y)
        data+=QByteArray::number(3.05+x*.1,'g',9)+" "+QByteArray::number(3.05+y*.1,'g',9)+" 3 99\n";
    if(!tagged.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    tagged.write(data);tagged.close();filtering.wallMargin=.01;
    // 此測試只隔離連通性；密度門檻由下方獨立測試驗證。 / Isolate connectivity here; density is tested separately below.
    filtering.floating.minPoints=1;
    const auto cleaned=processCloud(path,0,{},o,true,filtering);
    if(!cleaned.error.isEmpty() || cleaned.floatingRemoved!=25 || cleaned.reducedCount!=72 || cleaned.filledCount!=9) {
        std::cerr<<"Pre-fill floating cleanup regression: "<<cleaned.error.toStdString()<<" removed="<<cleaned.floatingRemoved<<" retained="<<cleaned.reducedCount<<" filled="<<cleaned.filledCount<<"\n";return 1;
    }
    for(const auto &p:cleaned.after) if(p.intensity!=42) return 1;
    // 固定場域旁的處理邊界可延伸到舊短牆以外，並須拒絕其他場域的邊界。
    // A processing profile beside the fixed model restores wings beyond old end walls and rejects a foreign site profile.
    QFile region(dir.filePath("processing_region.json"));if(!region.open(QIODevice::WriteOnly)) return 1;
    const QByteArray profile=R"({"version":1,"site_id":"synthetic-site","footprint_xy_m":[[0,0],[4,0],[4,-4],[10,-4],[10,26],[6,26],[6,24],[4,24],[4,20],[0,20]]})";
    region.write(profile);region.close();
    data.replace("element vertex 97","element vertex 99");data+="8.05 -2 1 42\n2 -2 1 99\n";
    if(!tagged.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    tagged.write(data);tagged.close();filtering.floating.enabled=false;
    const auto wings=processCloud(path,0,{},o,true,filtering);
    if(!wings.error.isEmpty() || wings.footprintSource.isEmpty() || wings.restoredEnds[0]!=1 || wings.restoredFootprintPoints!=1) {
        std::cerr<<"Persisted footprint regression: "<<wings.error.toStdString()<<"\n";return 1;
    }
    // 顯示預算不得影響體積；相同完整量測以 1 點與全點預覽必須得到同一數值。
    // The display cap must not affect volume: one displayed point and a full preview must use identical measurements.
    VolumeOptions volumes;volumes.enabled=true;volumes.automaticGround=false;volumes.ground={0,0,0};volumes.minimumArea=0;
    const auto volumeFull=processCloud(path,0,{},o,true,filtering,volumes);
    const auto volumeLimited=processCloud(path,1,{},o,true,filtering,volumes);
    if(!volumeFull.volume.error.isEmpty() || !volumeLimited.volume.error.isEmpty() ||
       std::abs(volumeFull.volume.totalVolume-volumeLimited.volume.totalVolume)>1e-12 || volumeFull.volume.totalVolume<=0) return 1;
    if(!volumeFull.meshVolume.error.isEmpty()||!volumeLimited.meshVolume.error.isEmpty()||std::abs(volumeFull.meshVolume.totalVolume-volumeLimited.meshVolume.totalVolume)>1e-12) return 1;
    auto rasterOnly=volumes;rasterOnly.mesh.enabled=false;const auto withoutMesh=processCloud(path,1,{},o,true,filtering,rasterOnly);
    if(!withoutMesh.error.isEmpty()||withoutMesh.meshVolume.enabled||std::abs(withoutMesh.volume.totalVolume-volumeFull.volume.totalVolume)>1e-12) return 1;
    // 確認柱體先移除量測，再由周邊支點建立虛擬料面計算交集，最後刪除柱內補點。
    // Confirmed columns remove observations first; virtual interpolated material supports volume intersection before estimates are removed.
    ColumnSemantic column;column.id=1;column.state=ColumnState::Confirmed;column.top=10;column.footprint=columnSection("rectangle",.65,1.15,.65,1.15);
    saveColumnSemantics(dir.filePath("column_semantics.json"),"synthetic-site",{column});filtering.columns.enabled=true;
    const auto columnResult=processCloud(path,0,{},o,true,filtering,volumes);
    if(!columnResult.error.isEmpty() || !columnResult.volume.error.isEmpty() || columnResult.volume.columnExcludedVolume<=0 ||
       columnResult.columns.removedEstimated!=9 || columnResult.filledCount!=0) {
        std::cerr<<"Column pipeline regression: "<<columnResult.error.toStdString()<<" "<<columnResult.volume.error.toStdString()<<"\n";return 1;
    }
    // 回歸：與料面分離的稀疏細柱會被去雜點及結構裁切刪除，但必須先參與柱體偵測。
    // Regression: a disconnected thin column removed by denoising/shell filtering must still supply semantic evidence first.
    auto withColumn=data;withColumn.replace("element vertex 99","element vertex 140");
    for(int z=0;z<=40;++z) withColumn+="5.1 5.1 "+QByteArray::number(z*.2,'g',9)+" 77\n";
    if(!tagged.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    tagged.write(withColumn);tagged.close();
    DenoiseOptions sparseCleanup{true,.5,4,true};
    const auto rawColumn=processCloud(path,1,sparseCleanup,{},true,filtering);
    bool detectedRawColumn=false;
    for(const auto &c:rawColumn.columns.items) {const auto b=columnBounds(c);
        if(b[0]<=5.1&&b[1]>=5.1&&b[2]<=5.1&&b[3]>=5.1&&c.observedTop-c.observedBottom>7.9) detectedRawColumn=true;}
    bool retainedEvidence=false;for(const auto &p:rawColumn.columnEvidence) if(p.intensity==77&&p.z>7) retainedEvidence=true;
    if(!rawColumn.error.isEmpty()||!detectedRawColumn||!retainedEvidence||!rawColumn.denoised) {
        std::cerr<<"Pre-denoise column evidence regression: "<<rawColumn.error.toStdString()<<"\n";return 1;
    }
    // 定案檔不再從同一掃描加入新候選，原審核幾何保留。
    // A finalized model must never repopulate fresh candidates from the same scan.
    column.group="R1";column.legacyIds={77};saveColumnSemantics(dir.filePath("column_semantics.json"),"synthetic-site",{column},true);
    const auto finalized=processCloud(path,1,sparseCleanup,{},true,filtering);
    if(!finalized.error.isEmpty()||!finalized.columns.finalized||finalized.columns.items.size()!=1||finalized.columns.items[0].group!="R1") {
        std::cerr<<"Final model repopulated candidates: "<<finalized.error.toStdString()<<"\n";return 1;
    }
    filtering.columns.enabled=false;
    if(!region.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    auto wrong=profile;wrong.replace("synthetic-site","foreign-site");region.write(wrong);region.close();
    if(processCloud(path,0,{},o,true,filtering).error.isEmpty()) return 1;
    // 懸浮點群即使本身密集，只要與主要料面分離，也應在補洞前被移除。
    // Even dense floating clusters should be removed before filling when detached from the main surface.
    std::vector<Point> clusters;
    for(int x=0;x<20;++x) for(int y=0;y<20;++y) clusters.push_back({x*.1f,y*.1f,-5,42});
    for(int x=0;x<5;++x) for(int y=0;y<5;++y) clusters.push_back({x*.1f,y*.1f,0,99});
    clusters.push_back({8,8,-5,99});
    const auto originalClusters=clusters;
    DenoiseOptions cleanup{true,.5,2,true};
    if(removeFloatingNoise(clusters,cleanup)!=26 || clusters.size()!=400) return 1;
    for(const auto &p:clusters) if(p.intensity!=42 || p.z!=-5) return 1;
    auto allDense=originalClusters;cleanup.largestComponent=false;
    if(removeFloatingNoise(allDense,cleanup)!=1 || allDense.size()!=425) return 1;
    auto disabledClusters=originalClusters;cleanup.enabled=false;
    if(removeFloatingNoise(disabledClusters,cleanup)!=0 || disabledClusters.size()!=426) return 1;
    std::vector<Point> empty;cleanup.enabled=true;
    if(removeFloatingNoise(empty,cleanup)!=0) return 1;
    auto sparse=std::vector<Point>{{0,0,0},{10,10,10}};
    if(removeFloatingNoise(sparse,cleanup)!=2 || !sparse.empty()) return 1;
    std::cout << "PASS: slope reconstruction, closed holes only, size/slope/residual limits, layers, ROI, resource cap, display budget\n";
    return 0;
}
