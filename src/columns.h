/** @file columns.h
 * @brief 可審核的柱體幾何語意與排除模型。 / Reviewable column geometry semantics and exclusion models.
 */
#pragma once
#include "footprint.h"
#include <QString>
#include <vector>
#include <cstddef>
#include <utility>
struct Point;

/** @brief 候選須由使用者確認才排除；非柱體可記錄以避免重複提議。 / Candidates require confirmation; rejected candidates are remembered. */
enum class ColumnState { Pending, Confirmed, Rejected };
/** @brief 柱／梁截面與上下界；舊柱預設從地面至頂高。 / Column/beam sections and bounds; legacy columns extend from ground to top. */
struct ColumnSemantic {
    int id=0;
    ColumnState state=ColumnState::Pending;
    QString shape="rectangle"; ///< rectangle / circle / polygon；工字鋼須提供實際實心截面。 / I-beams require their actual solid polygon.
    Footprint footprint{};
    double top=0, observedBottom=0, observedTop=0, verticalCoverage=0;
    bool groundRelative=false; ///< true 時上下界皆為離地高度。 / Bounds are ground-relative when true.
    double bottomHeight=0,topHeight=0;
    std::array<double,3> referenceGround{0,0,0}; ///< 清理／預覽參考；積分使用當次地面。 / Cleanup/preview reference; integration uses current ground.
    QString kind="column"; ///< column / beam / box.
    size_t evidencePoints=0;
    QString group; ///< 定案分群：W_S、R1、S1–S8、E、W_N。 / Final site groups.
    std::vector<int> legacyIds; ///< 舊審查 ID，只供追溯。 / Previous review IDs, for audit only.
    QString provenance; ///< 人工調整、碎片合併或規律推估來源。 / Review edits, merged fragments, or pattern inference origin.
    Footprint evidenceFootprint; ///< 原始柱段範圍；貼牆延伸後仍可匹配與檢查實測。 / Original measured bounds for matching/review after wall extension.
    double excludedVolume=0; ///< 本次有料面支撐的扣除量，不寫入永久語意。 / Current supported deduction, never persisted as geometry.
};
struct ColumnOptions {
    bool enabled=false;
    double xyCell=.25,zCell=.5,minHeight=4.,minCoverage=.65,maxWidth=2.5;
};
struct ColumnReport {
    std::vector<ColumnSemantic> items;
    QString modelPath,siteId,error;
    bool finalized=false; ///< 定案模型停用自動提議。 / Final models disable automatic proposals.
    size_t removedMeasured=0,removedEstimated=0;
};
/** @brief 從完整結構量測找垂直持續、窄截面的候選，不把自動候選視為確認語意。
 *         Detect persistent narrow vertical structures from full measurements; detections are never confirmed automatically.
 */
std::vector<ColumnSemantic> detectColumnCandidates(const std::vector<Point> &points,const ColumnOptions &options);
/** @brief 讀取固定場域旁的人工審核語意，驗證場域 ID 與確認截面。 / Load reviewed semantics and validate the site ID and confirmed sections. */
std::vector<ColumnSemantic> loadColumnSemantics(const QString &path,const QString &siteId);
/** @brief 原子保存人工審核結果；此檔案不修改校準模型。 / Atomically save reviewed semantics independently of calibration. */
void saveColumnSemantics(const QString &path,const QString &siteId,const std::vector<ColumnSemantic> &columns,bool finalized=false);
/** @brief 載入場域綁定的定案旗標。 / Read the site-bound finalization flag. */
bool columnSemanticsFinalized(const QString &path,const QString &siteId);
/** @brief 核對 ID、截面及上下界；幾何重疊由交集聯集積分處理。 / Validate IDs, sections and bounds; geometric overlaps are handled by union integration. */
void validateColumns(const std::vector<ColumnSemantic> &columns);
/** @brief 指定地面上的柱／梁上下界；公尺制，離地模式隨地面平移。 / Bounds in meters; relative solids follow the ground plane. */
double columnBottomZ(const ColumnSemantic &c,double x,double y);
double columnTopZ(const ColumnSemantic &c,double x,double y);
/** @brief 更新離地實體的清理與顯示參考，不改永久幾何。 / Update cleanup/display ground without altering persistent geometry. */
void setColumnGround(std::vector<ColumnSemantic> &columns,const std::array<double,3> &ground);
/** @brief 取得截面的 xmin,xmax,ymin,ymax。 / Return cross-section xmin,xmax,ymin,ymax. */
std::array<double,4> columnBounds(const ColumnSemantic &column);
/** @brief 建立矩形或圓形截面；圓以 64 邊形近似。 / Construct rectangular/circular sections; circles use 64 sides. */
Footprint columnSection(const QString &shape,double xMin,double xMax,double yMin,double yMax);
/** @brief 只刪除確認實體上下界內的點，保留梁下空間。 / Remove points only within confirmed solid bounds, preserving beam clearance. */
size_t removeColumnPoints(std::vector<Point> &points,const std::vector<ColumnSemantic> &columns);
/** @brief 一格內柱體與料柱的交集，未知料面不在此函數推估。 / Intersection of columns and one material cell; this function never infers missing surfaces. */
struct ColumnCellIntersection {double volume=0,fullHeightArea=0;std::vector<std::pair<int,double>> contributions;};
/** @brief 預先三角化確認截面並以 1 m 網格索引，避免每格遍歷所有柱。
 *         Pre-triangulate confirmed sections and index by 1 m bins, avoiding all-column searches per raster cell.
 */
class ColumnVolumeMask {
    struct Triangle {std::array<double,2> a,b,c;double top;int id;bool relative;double lower,upper;std::array<double,3> reference;};
    std::vector<Triangle> triangles;
    std::vector<std::vector<size_t>> bins;
    int x0=0,y0=0,nx=0,ny=0;
public:
    explicit ColumnVolumeMask(const std::vector<ColumnSemantic> &columns);
    /** @brief 扣除凸投影料柱內的結構聯集；所有單位為公尺，surface/ground 為 Z=aX+bY+c。
     *         Subtract the structural union inside a convex material prism; planes use Z=aX+bY+c in meters.
     * @details 逐實體從剩餘料體裁切，因此重疊只扣一次；逐柱歸屬依較小 ID 優先，總量不受輸入順序影響。
     *          Clip remaining material successively, counting overlaps once; smaller IDs own shared deductions.
     */
    ColumnCellIntersection intersectSurface(const Footprint &domain,const std::array<double,3> &surface,const std::array<double,3> &ground) const;
    bool contains(double x,double y,double z) const; ///< 確認柱體的量測點分類。 / Classification against confirmed solids.
    /** @brief 對 10 cm 料柱求截面交集 × 垂直交集高度；相交截面只扣一次。
     *         Integrate cross-section intersection times vertical overlap for a 10 cm material prism, without duplicate subtraction.
     */
    ColumnCellIntersection intersect(double x,double y,double size,double materialHeight,const std::array<double,3> &ground) const;
};

/** @brief 抽樣保留柱體候選的原始量測證據，不參與料堆計算。 / Sample measured candidate evidence independently of material calculations. */
std::vector<Point> sampleColumnEvidence(const std::vector<Point> &points,const std::vector<ColumnSemantic> &columns,size_t limit=100000);

/** @brief 依指定人工柱組補齊四角並統一截面；既有中心不移動。
 *         Complete an explicitly reviewed four-corner group and equalize sections without moving observed centers.
 * @param columns 完整場域語意，成功後原子替換；失敗不修改。 / Full site semantics, replaced only on success.
 * @param ids 二至四個確認柱 ID。 / Two to four confirmed member IDs.
 * @param positiveAxis 二柱時指定缺排朝 +X=0 或 +Y=1；三／四柱忽略。 / Missing row toward +X=0 or +Y=1 for pairs.
 * @param group 可追溯的組名。 / Auditable group name.
 * @return 含新增柱的四個組員 ID。 / Four member IDs including inferred additions.
 */
std::vector<int> completeColumnGroup(std::vector<ColumnSemantic> &columns,const std::vector<int> &ids,int positiveAxis,const QString &group);

/** @brief 將四柱組擬合到使用者量測的長方形；只用指定實測柱定位，不以推估柱當錨點。
 *         Fit four members to a measured axis-aligned rectangle, using only designated observed anchors.
 * @param northSouth 建築 X 方向中心距（公尺）。 / Center spacing along building X (meters).
 * @param eastWest 建築 Y 方向中心距（公尺）。 / Center spacing along building Y (meters).
 * @param section 正方形柱截面邊長（公尺）。 / Square section side length (meters).
 * @details 保留 ID、狀態、柱頂與原量測範圍；先驗證幾何，失敗不修改。
 *          Preserve IDs, states, tops and original evidence; validate before committing, with no mutation on failure.
 */
void fitMeasuredColumnGroup(std::vector<ColumnSemantic> &columns,const std::vector<int> &ids,const std::vector<int> &anchorIds,double northSouth,double eastWest,double section,const QString &name);
