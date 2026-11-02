/**
 * @file cloud.h
 * @brief 點雲資料結構與處理介面。 / Point-cloud data structures and processing API.
 */
#pragma once
#include <QString>
#include "footprint.h"
#include "volume.h"
#include "mesh_volume.h"
#include "columns.h"
#include <array>
#include <vector>
#include <cstddef>
#include <limits>

/** @brief 公尺制 XYZ 與強度；NaN 強度表示缺值。 / XYZ in meters and intensity; NaN intensity denotes missing data. */
struct Point { float x, y, z; float intensity = std::numeric_limits<float>::quiet_NaN(); };
/** @brief 空間密度與連通性去雜點設定。 / Spatial density and connectivity denoising settings. */
struct DenoiseOptions {
    bool enabled = false; ///< API 預設不去雜點，GUI 預設啟用。 / Disabled by default in the API, enabled by default in the GUI.
    double cellSize = .5; ///< 密度網格邊長（公尺）。 / Density grid cell size in meters.
    size_t minPoints = 3; ///< 每格至少需多少個原始量測點。 / Minimum measured points per cell.
    bool largestComponent = true; ///< 只保留點數最多的 26 鄰接連通區。 / Keep only the largest 26-connected component by point count.
};
/** @brief 體積前的保守局部料面離群清理。 / Conservative local-surface filtering before volume integration. */
struct SurfaceFilterOptions {bool enabled=true;double cell=.5,radius=1.5,tolerance=.6;};
struct SurfaceFilterReport {size_t removedMeasured=0,removedEstimated=0;bool comparisonValid=false;double volumeBefore=0,volumeAfter=0;};
/** @brief 用實測鄰格的穩健平面檢查原始點及補點，不移動保留點；刪除點另存供同地面體積比較。
 *         Check measurements/estimates against robust measured neighborhoods without moving survivors; return rejects for same-ground volume comparison.
 */
SurfaceFilterReport filterSurfaceOutliers(std::vector<Point> &measured,std::vector<Point> &estimated,const SurfaceFilterOptions &options,std::vector<Point> &rejectedMeasured,std::vector<Point> &rejectedEstimated);
/** @brief 建築外殼及料堆上方結構物的移除設定。 / Building-shell and overhead-structure removal settings. */
struct StructureOptions {
    bool enabled = false; ///< API 預設關閉；GUI 預設啟用。 / Disabled in the API; enabled in the GUI.
    double wallMargin = .4; ///< 牆面內縮距離，公尺。 / Interior wall setback in meters.
    double roofMargin = 2.; ///< 屋頂下方排除距離，公尺。 / Roof exclusion distance in meters.
    double cellSize = .5; ///< XY 網格邊長，公尺。 / XY cell width in meters.
    double surfaceBand = .8; ///< 同格最低點上方保留厚度，公尺。 / Retained thickness above the cell minimum in meters.
    Footprint footprint{}; ///< 真實料場凹多邊形；空值沿用四牆裁切。 / Actual concave site footprint; empty uses the four walls.
    ColumnOptions columns{}; ///< 柱體提議與人工確認語意。 / Column proposals and reviewed semantics.
    SurfaceFilterOptions surfaceFilter{}; ///< 補洞後、估地面前的末段清理。 / Final cleanup after filling and before fitting ground.
    DenoiseOptions floating{true,.5,2,true}; ///< 結構移除後的獨立清理設定。 / Independent cleanup after structure removal.
};
/** @brief 在已校準座標內清除建築外殼，再保留每個 XY 格的下層表面。
 *         Remove the calibrated building shell, then retain the lower measured surface in each XY cell.
 * @param points 原地壓縮，保留點的座標與強度不變。 / Compacted in place; retained coordinates and intensity are unchanged.
 * @param planes 固定場域模型的四面牆及兩面屋頂。 / Four wall and two roof planes from the fixed site model.
 * @param options 邊界及下層表面參數。 / Boundary and lower-surface parameters.
 * @return 移除點數。 / Number of removed points.
 * @warning 同格無下層量測時無法辨識所有懸空物件；貼牆或貼頂料堆也可能被裁除。
 *          Overhead objects without lower observations can remain; material near walls or roofs can be cropped.
 */
size_t removeBuildingStructures(std::vector<Point> &points, const std::vector<std::array<double,4>> &planes, const StructureOptions &options);
/** @brief 結構移除後，清除稀疏格及與主體分離的點群，原地保留座標與強度。
 *         After structure removal, discard sparse cells and detached clusters without changing retained coordinates/intensities.
 * @param points 完整量測資料，不包含補點。 / Full measured points, excluding estimates.
 * @param options 獨立的密度與 26 鄰接連通區設定。 / Independent density and 26-connected component settings.
 * @return 移除點數；停用時為零。 / Removed count; zero when disabled.
 * @warning 最大連通區模式也會移除分離的真實料堆，必要時停用此模式。
 *          Largest-component mode also removes detached real piles; disable this mode when needed.
 */
size_t removeFloatingNoise(std::vector<Point> &points, const DenoiseOptions &options);
/** @brief 局部或全場域的高度表面補洞設定。 / Local or full-domain height-surface hole-filling settings. */
enum class FillMethod { LocalPlane, LayeredTin };
struct FillOptions {
    bool enabled = false; ///< API 預設停用；應用程式預設啟用全場域。 / Disabled by default in the API; the application enables full-domain filling by default.
    double xMin = 0, xMax = 0, yMin = 0, yMax = 0, zMin = 0, zMax = 0;
    double spacing = .1; ///< 補點 XY 間距（公尺）。 / Generated XY spacing in meters.
    double maxHole = 8; ///< 缺口包圍盒對角線上限（公尺）。 / Maximum hole bounding-box diagonal in meters.
    double maxError = .75; ///< 邊界與擬合平面的最大高度誤差（公尺）。 / Maximum boundary-to-plane height residual in meters.
    double maxSlope = 3; ///< 最大坡度 dz / 水平距離。 / Maximum dz per horizontal distance.
    bool fullDomain = false; ///< 由清理後點雲自動決定全場域範圍。 / Derive full-domain bounds from the cleaned cloud.
    FillMethod method = FillMethod::LocalPlane; ///< 保留舊 API 行為；應用程式使用分層 TIN。 / Preserve legacy API behavior; the application uses layered TIN.
    double supportSpacing = .5; ///< 三角網格量測支點間距，不是輸出點距。 / Measured TIN support spacing, distinct from output spacing.
    double layerHeight = 8; ///< 重疊高度帶厚度（公尺）。 / Overlapping height-band thickness in meters.
    Footprint footprint{}; ///< 補點必須位於實際料場內。 / Estimates must lie inside the actual footprint.
    double footprintMargin = 0; ///< 真實邊界內縮距離，公尺。 / Actual perimeter setback in meters.
    bool validate = false; ///< 遮去已知區域驗證內插，不將測試補點當量測。 / Validate by withholding known patches, never treating predictions as measurements.
};
/** @brief 單一人工遮擋區域的量測對照誤差。 / Measured-reference errors for one artificial occlusion patch. */
struct ValidationPatch {
    Point center;
    double radius=0;
    size_t targets=0,predicted=0;
    double rmse=0,maxError=0;
};
/** @brief 補洞的幾何限制及已知區域遮擋驗證統計。 / Geometric limits and known-patch holdout validation statistics. */
struct FillReport {
    size_t rejectedFootprint = 0; ///< 跨越真實料場凹角的三角形。 / Triangles crossing actual site recesses.
    size_t bands = 0, supports = 0, triangles = 0, patches = 0, rejectedDistance = 0, rejectedSlope = 0, rejectedConsistency = 0;
    size_t validationPatches = 0, validationTargets = 0, validationPredicted = 0;
    double validationSquaredError = 0, validationAbsoluteError = 0, validationMaxError = 0;
    bool outputLimited = false;
    std::vector<ValidationPatch> validation;
};
/**
 * @brief 完整處理的統計與有限大小的顯示樣本。 / Full-processing statistics and bounded display samples.
 * @note before/after 不一定包含所有點；完整工作陣列在處理結束後釋放。
 *       before/after may omit points; the full working array is released after processing.
 */
struct Cloud {
    QString exportedPath, exportError; ///< 完整結果匯出路徑與獨立錯誤。 / Full-result export path and independent error.
    QString coordinateFrame; ///< 已套用的共用座標系，空值表示原始座標。 / Applied shared frame; empty means source coordinates.
    /// 清理前／後的顯示樣本。 / Display samples before/after cleaning.
    std::vector<Point> before, after;
    /// 與量測樣本分開保存的合成補點；強度為 NaN，不捏造量測值。 / Separate synthetic samples; intensity is NaN rather than an invented measurement.
    std::vector<Point> filled;
    size_t filledCount = 0, filledHoles = 0;
    FillOptions appliedFill; ///< 實際使用的補點範圍與參數。 / Actual fill bounds and parameters used.
    FillReport fillReport;
    ColumnReport columns;
    std::vector<Point> columnEvidence; ///< 結構移除前的柱體證據疊圖。 / Pre-removal evidence overlay, excluded from export/volume.
    bool hasRoofPlanes=false; ///< 僅供屋頂顯示裁切，不改變計算資料。 / Display-only roof clipping availability.
    std::array<std::array<double,4>,2> roofPlanes{}; ///< 朝上的單位法向屋頂平面。 / Roof planes with upward unit normals.
    SurfaceFilterReport surfaceFilter;
    MeshVolumeReport meshVolume; ///< 與柵格共用地面的獨立 TIN 積分。 / Independent TIN integration on the same ground.
    VolumeReport volume; ///< 完整清理量測與補點的體積報告。 / Volume report from full cleaned measurements and estimates.
    /// 去雜點後的幾何邊界，供選取補點範圍。 / Denoised geometry bounds for selecting a fill region.
    Point cleanLow{0,0,0}, cleanHigh{0,0,0};
    /// 有效輸入點數／清理後量測數，不受顯示上限影響。 / Valid input/cleaned measurement counts, independent of the display cap.
    size_t originalCount = 0, reducedCount = 0;
    size_t floatingRemoved = 0; ///< 結構移除後清除的懸浮／稀疏候選點數。 / Detached/sparse candidates removed after structure filtering.
    bool floatingFiltered = false; ///< 是否執行結構移除後的清理。 / Whether post-structure cleanup ran.
    size_t structureRemoved = 0; ///< 移除的結構候選點數。 / Removed structure candidate count.
    std::array<size_t,2> restoredEnds{0,0}; ///< 舊兩端短牆之外的保留點數。 / Retained counts beyond each original end wall.
    size_t restoredFootprintPoints = 0; ///< 清理後保留在舊矩形之外的量測點數。 / Cleaned measurements retained outside the old rectangle.
    QString footprintSource; ///< 獨立處理邊界檔案。 / Independent processing-footprint file.
    bool structuresFiltered = false; ///< 是否移除結構物。 / Whether structure filtering was applied.
    bool denoised = false; ///< 是否套用去雜點。 / Whether denoising was applied.
    /// 包含有效輸入點與原點的邊界。 / Bounds enclosing valid input points and the origin.
    Point low{0, 0, 0}, high{0, 0, 0};
    /// 完整有效輸入資料的著色範圍，兩側共用；高度範圍不強制包含原點。 / Shared full-input color ranges; height does not force inclusion of the origin.
    float zMin = 0, zMax = 0, intensityMin = 0, intensityMax = 0;
    /// 完整有效輸入的 P5/P95，線性內插；自動色階使用此範圍。 / Full valid-input P5/P95 with linear interpolation, used for automatic color limits.
    float zP05 = 0, zP95 = 0, intensityP05 = 0, intensityP95 = 0;
    /// 是否至少包含一個有效強度值。 / Whether at least one finite intensity value exists.
    bool hasIntensity = false;
    /// 空字串表示成功；失敗時其他欄位可能只完成部分填入。 / Empty on success; other fields may be partial on failure.
    QString error;
};
/**
 * @brief 讀取校準後原始量測，清理結構、補面、估地面並積分體積。
 *         Load calibrated raw measurements, clean structures, reconstruct surfaces, fit ground, and integrate volume.
 * @param path 單份 PLY 路徑。 / Path of one PLY scan.
 * @param displayLimit 每側顯示上限，0 為全部；不影響清理、體積或匯出。
 *                     Per-pane preview cap; zero means all, independent of cleaning, volume and export.
 * @param denoise 原始量測的密度／連通區清理。 / Density/connectivity filtering of raw measurements.
 * @param fill 清理後實測料面的補洞設定。 / Reconstruction options for cleaned measured material surfaces.
 * @param applyCalibration 在所有清理之前套用共同座標。 / Apply the common frame before any cleaning.
 * @param structures 建築及確認柱體清理；體積須有明確場域邊界。
 *                   Building/confirmed-column filtering; volume requires an explicit site footprint.
 * @param volume 10 cm XY 格最低離地高度積分設定。 / Minimum-height integration on a 10 cm XY grid.
 * @param exportPath 可選完整處理點雲輸出路徑。 / Optional full processed-cloud export path.
 * @return 原量測及處理結果的顯示樣本、統計與錯誤。 / Display samples, statistics and errors.
 * @note 不進行三維體素降取樣或平均座標，保留實測最低高度。
 *       Never voxel-average/downsample measurements, preserving measured lower heights.
 */
Cloud processCloud(const QString &path, size_t displayLimit, const DenoiseOptions &denoise = {}, const FillOptions &fill = {}, bool applyCalibration = true, const StructureOptions &structures = {}, const VolumeOptions &volume = {}, const QString &exportPath = {});

/**
 * @brief 以封閉缺口邊界擬合局部平面產生估計表面點。 / Fit local planes to enclosed hole boundaries and generate estimated surface points.
 * @param points 完整去雜點結果，保持不變。 / Full denoised points, left unchanged.
 * @param options 使用者指定 ROI 與保守的缺口限制。 / User-selected ROI and conservative hole limits.
 * @param holes 回傳成功補上的缺口數。 / Receives the number of successfully filled holes.
 * @return 合成點，強度為 NaN；不保證還原真實遮擋表面。 / Synthetic points with NaN intensity; actual occluded geometry is not guaranteed.
 * @throws std::runtime_error 無效參數或超過網格資源上限。 / Invalid parameters or grid resource limits exceeded.
 */
std::vector<Point> fillStockpileHoles(const std::vector<Point> &points, const FillOptions &options, size_t &holes);

/**
 * @brief 在重疊高度帶內建立受距離與坡度限制的 Delaunay 三角網格，以重心座標補點。
 *        Build distance/slope-bounded Delaunay TINs in overlapping height bands and fill using barycentric interpolation.
 * @param points 完整清理後量測資料，不修改。 / Full cleaned measurements, left unchanged.
 * @param options 範圍、網格、分層與驗證設定。 / Bounds, grid, layering, and validation settings.
 * @param report 回傳限制與遮擋驗證統計。 / Receives limit and holdout-validation statistics.
 * @return 估計點，強度為 NaN；只補三角形內部，不在凸包外外插。
 *         Estimated points with NaN intensity; interpolate only inside accepted triangles, never outside the convex hull.
 * @warning 三角網格不會辨識原料、建築或真實空洞；大遮擋內部的真實峰谷仍不可由單次掃描確定。
 *          TINs do not recognize material, buildings, or real voids; hidden peaks/valleys remain indeterminate from one scan.
 */
std::vector<Point> reconstructStockpileSurface(const std::vector<Point> &points, const FillOptions &options, FillReport &report);

/** @brief 匯出完整量測與補點為 binary little-endian PLY；不覆寫既有檔。
 *         Export full measurements and estimates as binary little-endian PLY without overwriting existing files.
 * @param path 新 PLY 路徑，父目錄須存在。 / New PLY path; its parent directory must exist.
 * @param measured 完整清理量測，不是顯示樣本。 / Full cleaned measurements, not display samples.
 * @param estimated 完整合成補點。 / Full synthetic estimates.
 * @param calibrated 是否已為建築座標；標頭可防止再次校準。 / Whether coordinates are calibrated; the header prevents re-registration.
 * @throws std::runtime_error 輸出已存在或寫入失敗。 / Existing output or write failure.
 */
void exportProcessedPly(const QString &path,const std::vector<Point> &measured,const std::vector<Point> &estimated,bool calibrated);
