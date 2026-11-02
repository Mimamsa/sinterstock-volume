/** @file calibration_transform.h
 * @brief 共用座標系的剛體轉換載入與套用。 / Loading and applying rigid transforms into a shared coordinate frame.
 */
#pragma once
#include <QString>
#include <array>

/** @brief 列優先 4x4 剛體矩陣；公尺制、無縮放。 / Row-major 4x4 rigid transform, meters, without scaling. */
struct CalibrationTransform {
    std::array<double,16> matrix{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    bool applied=false;
    QString frame;
};
/** @brief 讀取既有矩陣；場域模型存在時自動註冊新掃描，否則保留原始座標。
 *         Load a cached matrix, automatically register new scans if a site model exists, or retain source coordinates otherwise.
 * @param plyPath 原始 PLY 路徑。 / Source PLY path.
 * @return 已驗證的剛體矩陣與座標系名稱。 / Validated rigid matrix and frame label.
 * @throws std::runtime_error 側錄無效、來源已變更且無模型、註冊品質不足或方向不明。
 *                           Invalid sidecar, stale source without a model, inadequate registration quality, or ambiguous orientation.
 */
CalibrationTransform loadCalibrationTransform(const QString &plyPath);
