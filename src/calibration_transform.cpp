/** @file calibration_transform.cpp
 * @brief 校準側錄格式與剛體矩陣驗證。 / Calibration sidecar parsing and rigid-matrix validation.
 */
#include "calibration_transform.h"
#include "config_paths.h"
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <cmath>
#include <stdexcept>
#include <QCoreApplication>
#include <QProcess>
#include <QCryptographicHash>

CalibrationTransform loadCalibrationTransform(const QString &plyPath) {
    CalibrationTransform result;
    const QFileInfo info(plyPath);
    /** @brief 已輸出的校準 PLY 帶有明確標記，不得因找到場域模型而再次註冊。
     *         Exported calibrated PLYs carry a marker and must not be registered again merely because a site model is available.
     */
    QFile tagged(plyPath);
    if(tagged.open(QIODevice::ReadOnly)) {
        const auto prefix=tagged.read(65536);const auto end=prefix.indexOf("end_header");
        if(end>=0 && prefix.left(end).contains("comment calibrated building frame; meters; origin at reference roof corner")) {
            result.applied=true;result.frame="Pre-calibrated building coordinates (identity transform)";return result;
        }
    }
    const auto manifestPath=info.dir().filePath("calibration.json");
    QJsonObject object;
    auto readManifest=[&] {
        QFile file(manifestPath); if (!file.exists()) return QJsonObject{};
        if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot read calibration.json");
        QJsonParseError error; const auto doc=QJsonDocument::fromJson(file.readAll(),&error);
        if (error.error!=QJsonParseError::NoError || !doc.isObject() || doc.object().value("version").toInt()!=1)
            throw std::runtime_error("Invalid calibration.json format/version");
        return doc.object();
    };
    object=readManifest();
    auto record=object.value("files").toObject().value(info.fileName()).toObject();
    const QString site=resolveSiteModel(info.dir(),object.value("site_model").toString(),QCoreApplication::applicationName()=="pointcloud_compare");
    bool stale=!record.isEmpty() && record.contains("source_bytes") && qint64(record.value("source_bytes").toDouble())!=info.size();
    QByteArray digest;
    if (record.contains("source_sha256")) {
        QFile input(plyPath);if(!input.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot hash source PLY");
        QCryptographicHash hash(QCryptographicHash::Sha256);if(!hash.addData(&input)) throw std::runtime_error("Cannot verify source SHA-256");
        digest=hash.result().toHex();stale=stale || QString::fromLatin1(digest)!=record.value("source_sha256").toString();
    }
    /** @brief 新檔案或來源變更時，用固定場域在背景註冊；失敗不得悄悄顯示成同一座標。
     *         Register new/changed sources against the fixed site in the background; never silently label a failed scan as aligned.
     */
    if (record.isEmpty() || stale) {
        if (site.isEmpty()) {if(stale) throw std::runtime_error("Source changed; provide its immutable site model to recalibrate");return result;}
        const auto executable=QDir(QCoreApplication::applicationDirPath()).filePath("calibrate_clouds");
        if(!QFile::exists(executable)) throw std::runtime_error("Build calibrate_clouds to register this new scan");
        QProcess process;
        process.start(executable,{QFileInfo(plyPath).absoluteFilePath(),"--site",QFileInfo(site).absoluteFilePath(),"--matrix-only"});
        if(!process.waitForStarted(5000) || !process.waitForFinished(120000)) {process.kill();process.waitForFinished();throw std::runtime_error("Automatic fixed-site registration timed out");}
        if(process.exitStatus()!=QProcess::NormalExit || process.exitCode()!=0)
            throw std::runtime_error(("Automatic calibration rejected scan: "+QString::fromUtf8(process.readAllStandardError())).toStdString());
        object=readManifest();record=object.value("files").toObject().value(info.fileName()).toObject();
        if(record.isEmpty()) throw std::runtime_error("Automatic registration produced no transform");
    }
    const auto m=record.value("source_to_building").toArray();
    if (m.size()!=16) throw std::runtime_error("Calibration matrix requires 16 row-major values");
    for (int i=0;i<16;++i) {
        if (!m[i].isDouble() || !std::isfinite(m[i].toDouble())) throw std::runtime_error("Invalid calibration coefficient");
        result.matrix[i]=m[i].toDouble();
    }
    const auto &a=result.matrix;
    if (std::abs(a[12])+std::abs(a[13])+std::abs(a[14])+std::abs(a[15]-1)>1e-8) throw std::runtime_error("Invalid homogeneous matrix row");
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        double dot=0; for (int k=0;k<3;++k) dot+=a[4*i+k]*a[4*j+k];
        if (std::abs(dot-(i==j ? 1.:0.))>1e-5) throw std::runtime_error("Calibration rotation is not orthonormal");
    }
    const double determinant=a[0]*(a[5]*a[10]-a[6]*a[9])-a[1]*(a[4]*a[10]-a[6]*a[8])+a[2]*(a[4]*a[9]-a[5]*a[8]);
    if (std::abs(determinant-1)>1e-5) throw std::runtime_error("Calibration must preserve right-handed coordinates");
    result.applied=true;
    result.frame=object.value("frame").toString("Building roof corner");
    return result;
}
