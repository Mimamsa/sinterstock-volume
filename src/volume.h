/** @file volume.h
 * @brief 地平面、料堆邊界及固定 10 cm 柵格體積。 / Ground plane, pile boundaries, and fixed 10 cm raster volumes.
 */
#pragma once
#include "footprint.h"
#include "columns.h"
#include <QString>
#include <array>
#include <vector>
#include <cstddef>
#include <optional>
struct Point;

/** @brief 實驗性料面網格；支點間距不是體積積分格距。 / Experimental surface mesh; support spacing is not a volume raster step. */
struct MeshVolumeOptions {bool enabled=true;double supportSpacing=.25,maxEdge=2.,maxSlope=3.;size_t maxVertices=250000;};
/** @brief 體積參數；所有距離以公尺表示。 / Volume options; all distances are meters. */
struct VolumeOptions {
    bool enabled=false; ///< API 預設停用，應用程式預設啟用。 / Disabled in the API, enabled in the application.
    bool automaticGround=true; ///< 只由量測估計地面，不使用補點。 / Estimate ground from measurements only, never estimates.
    std::array<double,3> ground{0,0,0}; ///< 手動地面 z=a*x+b*y+c。 / Manual ground z=a*x+b*y+c.
    std::optional<double> groundDistance; ///< 原點至地面法線距離 d（m，帶正負號）；保留擬合 a、b，依 d 正負號決定地面位於原點上方或下方。 / Signed origin-to-ground normal distance d (m); retain fitted a,b and use the sign of d to place ground above/below origin.
    double groundTolerance=.15; ///< 地面內點垂直容差。 / Vertical ground inlier tolerance.
    double pileThreshold=.2; ///< 高於地面此高度才列入料堆。 / Minimum above-ground height classified as material.
    std::vector<ColumnSemantic> columns{}; ///< 只扣除確認截面內與料柱相交的實心柱體。 / Deduct only confirmed solid intersections.
    MeshVolumeOptions mesh{};
    double minimumArea=1.; ///< 最小料堆投影面積，平方公尺。 / Minimum pile projected area, square meters.
};
/** @brief 一個 8 鄰接料堆的柵格積分。 / Raster integral of one 8-connected pile. */
struct PileVolume {
    size_t id=0, measuredCells=0, estimatedCells=0, enclosedUnknownCells=0;
    /// meanHeight 為排柱後淨體積 / 淨投影面積，表示等效料高。 / meanHeight is net volume/net projected area after column exclusion.
    double area=0, volume=0, measuredVolume=0, estimatedVolume=0, meanHeight=0, maxHeight=0;
    double centerX=0,centerY=0,boundaryLength=0;
};
/** @brief 地面上的柵格邊界線，供右側視覺化；不參與數值積分。 / Raster boundary edge on ground for right-pane rendering; never used for integration. */
struct VolumeBoundary { float x1,y1,x2,y2;size_t pile=0; };
/** @brief 體積與資料覆蓋統計；未知格不當作零高度或完整體積。 / Volumes and coverage; unknown cells are neither zero height nor complete volume. */
struct VolumeReport {
    bool enabled=false, groundValid=false, manualGround=false;
    std::array<double,3> ground{0,0,0};
    std::array<double,3> fittedGround{0,0,0}; ///< 原始擬合平面；RMSE 與內點統計屬於此平面。 / Original fit underlying RMSE and inlier statistics.
    bool adjustedGroundDistance=false;
    double groundRMSE=0,groundSpanX=0,groundSpanY=0;
    size_t groundSamples=0,groundInliers=0;
    static constexpr double spacing=.1; ///< 體積網格固定 10 cm，獨立於補點與顯示格距。 / Fixed 10 cm volume grid, independent of fill/display spacing.
    size_t footprintCells=0,coveredCells=0,unknownCells=0,groundCells=0,belowGroundCells=0,discardedSmallCells=0;
    double totalVolume=0,totalArea=0;
    std::vector<std::pair<int,double>> columnVolumes; ///< 逐柱扣除量，與全域扣除量同一積分。 / Per-column deductions from the same integral.
    double grossRasterVolume=0,columnExcludedVolume=0,columnExcludedArea=0,discardedSmallVolume=0;
    std::vector<PileVolume> piles;
    std::vector<VolumeBoundary> boundaries;
    bool boundaryDisplayLimited=false;
    QString error;
    double groundZ(double x,double y) const {return ground[0]*x+ground[1]*y+ground[2];}
};
/** @brief 估計地面，再以最低離地高度建立 10x10 cm XY 柵格、分割料堆並積分體積。
 *         Estimate ground, rasterize minimum above-ground heights on a 10x10 cm XY grid, segment piles, and integrate volumes.
 * @param measured 完整清理後量測，非顯示樣本。 / Full cleaned measurements, never display samples.
 * @param estimated 補點只支撐沒有量測的格子，不與量測重複加總。 / Estimates support only cells lacking measurements, avoiding double counting.
 * @param footprint 明確場域凹多邊形，不接受空邊界。 / Explicit concave site footprint; empty is rejected.
 * @param margin 真實邊界內縮距離。 / Perimeter setback.
 * @param options 地面及分割設定；groundDistance 在擬合驗證後，以 c = d*sqrt(1+a²+b²) 保留坡度並依 d 正負號設定地面。
 *                Ground/segmentation settings; groundDistance sets c=d*sqrt(1+a²+b²) after fit validation, preserving slopes and the sign of d.
 * @return 每格 V=0.01*min(z-ground(x,y)) 的累加與每堆報告；失敗不影響點雲顯示。
 *         Sum V=0.01*min(z-ground(x,y)) per cell and report per pile; errors do not prevent point rendering.
 * @warning 接觸且連通的料堆視為同一堆，地面被完全遮住時須提供手動地平面。
 *          Touching connected piles form one component; a fully occluded floor requires a manually supplied ground plane.
 * @note 調整 c 後，報告內 RMSE／內點統計仍描述 fittedGround；不可與完整手動地面併用。
 *       With c overridden, RMSE/inliers still describe fittedGround; a full manual plane cannot be combined with this override.
 */
VolumeReport calculateStockpileVolumes(const std::vector<Point> &measured,const std::vector<Point> &estimated,
                                      const Footprint &footprint,double margin,const VolumeOptions &options);
