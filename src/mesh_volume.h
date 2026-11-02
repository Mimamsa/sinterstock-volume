/** @file mesh_volume.h
 * @brief 獨立料面 TIN 積分、幾何交集及可檢查表面輸出。 / Independent TIN integration, solid intersections and inspectable surface export.
 */
#pragma once
#include "volume.h"
struct MeshSurfaceTriangle {std::array<std::array<double,3>,3> vertices;};
struct MeshVolumeReport {
    bool enabled=false;
    QString error;
    size_t vertices=0,inputTriangles=0,rejectedDistance=0,rejectedSlope=0,components=0,discardedComponents=0;
    double supportSpacing=0,footprintArea=0,coveredArea=0,materialArea=0,grossVolume=0,columnExcludedVolume=0,discardedSmallVolume=0,totalVolume=0,totalArea=0;
    std::array<double,3> ground{};
    std::vector<MeshSurfaceTriangle> surface; ///< 通過面積篩選的上表面，柱體仍以獨立交集扣除。 / Accepted upper surface, with columns subtracted analytically.
};
/** @brief 使用固定地平面積分 XY Delaunay TIN；不呼叫柵格積分或由柵格體積反推。
 *         Integrate an XY Delaunay TIN on a fixed ground, never rasterizing or deriving volume from raster results.
 * @details 量測支點優先；支點格只選原點，不平均。凹邊界／高度門檻先裁切，柱體交集分段精確積分。
 *          Prefer measured supports without averaging; clip concave boundaries/height thresholds and integrate piecewise column intersections.
 */
MeshVolumeReport calculateMeshVolume(const std::vector<Point> &measured,const std::vector<Point> &estimated,const Footprint &footprint,double margin,const std::array<double,3> &ground,const VolumeOptions &options);
/** @brief 匯出開放上表面供檢查，非已封底實體；不覆寫檔案。 / Export the open upper surface for inspection, not a capped solid; never overwrite files. */
void exportMeshSurface(const QString &path,const MeshVolumeReport &report);
