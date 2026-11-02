/** @file config_paths.h
 * @brief 固定場域定義檔的集中查找。 / Central lookup of predefined site configuration.
 */
#pragma once
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QStringList>
/** @brief 優先使用側錄提示及來源設定，再找工作目錄／執行檔旁的專案設定。
 *         Prefer manifest/source configuration, then project configuration near the cwd/executable.
 * @param source 掃描目錄。 / Scan directory.
 * @param hint 側錄中的模型路徑；失效時繼續查找。 / Manifest hint; stale paths fall through.
 * @param projectFallback 是否查找專案設定；原始資料 API 可停用。 / Enable project fallback; raw-data APIs may disable it.
 * @return 存在的絕對路徑，否則空字串。 / Existing absolute path, or an empty string.
 * @details 保留來源舊版 site_model.json 相容性；相關結構／邊界檔與模型同目錄。
 *          Legacy source-side site_model.json remains supported; companion models reside beside the site file.
 */
inline QString resolveSiteModel(const QDir &source,const QString &hint={},bool projectFallback=true) {
    QStringList candidates;
    if(!hint.isEmpty()) candidates<<hint;
    candidates<<source.filePath("config/site_model.json")<<source.filePath("site_model.json");
    if(projectFallback) candidates<<QDir::current().filePath("config/site_model.json")<<QDir::current().filePath("site_model.json")
        <<QDir(QCoreApplication::applicationDirPath()).filePath("../config/site_model.json")<<QDir(QCoreApplication::applicationDirPath()).filePath("../site_model.json");
    for(const auto &path:candidates) if(QFileInfo(path).isFile()) return QFileInfo(path).absoluteFilePath();
    return {};
}
