/** @file calibration_test.cpp
 * @brief 含料堆變動的雙掃描端對端校準測試。 / End-to-end two-scan calibration tests with changing piles.
 */
#include "cloud.h"
#include "config_paths.h"
#include "calibration_transform.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <QCoreApplication>
#include <QTemporaryDir>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QJsonDocument>
#include <QJsonObject>
#include <QCryptographicHash>
#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

struct TestRecord {float x,y,z,intensity;};

/** @brief 建立已知建築物、不同料堆與不同掃描姿態，檢查角點、剛體對齊及無損強度輸出。
 *         Create known architecture, changed piles, and different scan poses; check corners, rigid alignment, and lossless intensity.
 * @param argc 參數數量。 / Argument count.
 * @param argv 第二參數是校準執行檔路徑。 / The second argument is the calibration executable path.
 * @return 成功為 0。 / Zero on success.
 */
int main(int argc,char **argv) {
    QCoreApplication app(argc,argv); if (argc<2) return 1;
    QTemporaryDir temp; if (!temp.isValid()) return 1;
    // 集中設定查找、失效側錄與舊目錄相容。 / Config lookup, stale hints and legacy source layout.
    const QDir lookupDir(QDir(temp.path()).filePath("lookup"));
    if(!QDir().mkpath(lookupDir.filePath("config"))) return 1;
    QFile configSite(lookupDir.filePath("config/site_model.json"));if(!configSite.open(QIODevice::WriteOnly)) return 1;configSite.write("{}");configSite.close();
    if(resolveSiteModel(lookupDir,"/missing/model.json",false)!=QFileInfo(configSite).absoluteFilePath()) return 1;
    QFile legacySite(lookupDir.filePath("site_model.json"));if(!legacySite.open(QIODevice::WriteOnly)) return 1;legacySite.write("{}");legacySite.close();
    if(resolveSiteModel(lookupDir,legacySite.fileName(),false)!=QFileInfo(legacySite).absoluteFilePath()) return 1;
    if(!resolveSiteModel(QDir(lookupDir.filePath("empty")),{},false).isEmpty()) return 1;
    std::vector<Eigen::Vector3d> structure;
    for (int y=0;y<200;++y) for (int x=0;x<50;++x) {
        const double xx=(x+.5)*.4,yy=(y+.5)*.4;
        structure.emplace_back(xx,yy,xx<=13 ? 10+.2*xx : 12.6-.2*(xx-13));
    }
    for (int y=0;y<200;++y) for (int z=0;z<25;++z) {
        structure.emplace_back(0,(y+.5)*.4,(z+.5)*.4); structure.emplace_back(20,(y+.5)*.4,(z+.5)*.4);
    }
    for (int x=0;x<50;++x) for (int z=0;z<25;++z) {
        structure.emplace_back((x+.5)*.4,0,(z+.5)*.4); structure.emplace_back((x+.5)*.4,80,(z+.5)*.4);
    }
    const Eigen::Matrix3d rotations[]{(Eigen::AngleAxisd(.03,Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(.07,Eigen::Vector3d::UnitY())).toRotationMatrix(),
        (Eigen::AngleAxisd(3.08,Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(-.04,Eigen::Vector3d::UnitX())).toRotationMatrix()};
    const Eigen::Vector3d translations[]{{1,-4,2},{-3,2,1}};
    std::vector<std::vector<TestRecord>> records(2);
    std::array<QByteArray,2> hashes;
    std::mt19937 random(91); std::uniform_real_distribution<double> coordinate(0,1);
    for (int scan=0;scan<2;++scan) {
        auto world=structure;
        for (int i=0;i<700;++i) world.emplace_back(3+coordinate(random)*14,5+coordinate(random)*70,1+coordinate(random)*(scan+2));
        for (size_t i=0;i<world.size();++i) {
            const Eigen::Vector3d raw=rotations[scan]*world[i]+translations[scan];
            records[scan].push_back({float(raw.x()),float(raw.y()),float(raw.z()),float(i%1024+1)});
        }
        QFile file(QDir(temp.path()).filePath(QString("scan%1.ply").arg(scan)));
        if (!file.open(QIODevice::WriteOnly)) return 1;
        file.write(QString("ply\nformat binary_little_endian 1.0\nelement vertex %1\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n").arg(qulonglong(world.size())).toUtf8());
        file.write(reinterpret_cast<const char*>(records[scan].data()),qint64(records[scan].size()*sizeof(TestRecord))); file.close();
        if (!file.open(QIODevice::ReadOnly)) return 1; hashes[scan]=QCryptographicHash::hash(file.readAll(),QCryptographicHash::Sha256);
    }
    const auto sitePath=QDir(temp.path()).filePath("site_model.json");
    QProcess initialize;
    initialize.start(QString::fromLocal8Bit(argv[1]),{"--init-site",sitePath,"--reference",QDir(temp.path()).filePath("scan0.ply")});
    if (!initialize.waitForFinished(60000) || initialize.exitCode()) {std::cerr<<initialize.readAllStandardError().toStdString();return 1;}
    QFile siteBefore(sitePath);siteBefore.open(QIODevice::ReadOnly);const auto frozen=siteBefore.readAll();siteBefore.close();
    QProcess process; process.start(QString::fromLocal8Bit(argv[1]),{temp.path(),"--site",sitePath});
    if (!process.waitForFinished(60000) || process.exitCode()) {std::cerr<<process.readAllStandardOutput().toStdString()<<process.readAllStandardError().toStdString();return 1;}
    siteBefore.open(QIODevice::ReadOnly);if(siteBefore.readAll()!=frozen)return 1;siteBefore.close();
    for (int scan=0;scan<2;++scan) {
        const auto path=QDir(temp.path()).filePath(QString("scan%1.ply").arg(scan));
        QFile raw(path); raw.open(QIODevice::ReadOnly);
        if (QCryptographicHash::hash(raw.readAll(),QCryptographicHash::Sha256)!=hashes[scan]) return 1;
        const auto transform=loadCalibrationTransform(path); if (!transform.applied) return 1;
        Eigen::Matrix4d actual; for (int i=0;i<4;++i) for (int j=0;j<4;++j) actual(i,j)=transform.matrix[i*4+j];
        for (size_t i=0;i<structure.size();i+=61) {
            const auto &p=records[scan][i]; const Eigen::Vector3d aligned=actual.block<3,3>(0,0)*Eigen::Vector3d(p.x,p.y,p.z)+actual.block<3,1>(0,3);
            const Eigen::Vector3d expected=structure[i]-Eigen::Vector3d(0,0,10);
            if ((aligned-expected).norm()>.03) {std::cerr<<"Calibration geometry mismatch "<<(aligned-expected).norm()<<"\n";return 1;}
        }
        QFile out(QDir(temp.path()).filePath(QString("calibrated/scan%1.ply").arg(scan))); if (!out.open(QIODevice::ReadOnly)) return 1;
        while (!out.atEnd()) if (out.readLine().startsWith("end_header")) break;
        const auto bytes=out.readAll(); if (size_t(bytes.size())!=records[scan].size()*sizeof(TestRecord)) return 1;
        for (size_t i=0;i<records[scan].size();++i) {
            TestRecord p; std::memcpy(&p,bytes.constData()+i*sizeof(TestRecord),sizeof(p));
            if (p.intensity!=records[scan][i].intensity) return 1;
        }
        const auto calibrated=processCloud(path,100,{},{}),original=processCloud(path,100,{},{},false);
        const auto exported=processCloud(QDir(temp.path()).filePath(QString("calibrated/scan%1.ply").arg(scan)),100,{},{});
        if (!calibrated.error.isEmpty() || calibrated.coordinateFrame.isEmpty() || !original.coordinateFrame.isEmpty() ||
            calibrated.before.size()!=exported.before.size()) return 1;
        for (size_t i=0;i<calibrated.before.size();++i) {
            const auto a=calibrated.before[i],b=exported.before[i];
            if (std::abs(a.x-b.x)+std::abs(a.y-b.y)+std::abs(a.z-b.z)>1e-4) return 1;
        }
    }
    // 後來新增一份排序在舊檔之前、掃描方向不同且料堆改變的資料，座標系不得重設。
    // Add a later scan sorting before old files, with a different pose/pile; the frame must remain unchanged.
    const auto futurePath=QDir(temp.path()).filePath("a_future.ply");
    const Eigen::Matrix3d futureRotation=(Eigen::AngleAxisd(1.2,Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(.05,Eigen::Vector3d::UnitY())).toRotationMatrix();
    const Eigen::Vector3d futureTranslation(-8,7,-3);
    std::vector<TestRecord> future;
    for (size_t i=0;i<records[0].size();++i) {
        const auto &p=records[0][i];Eigen::Vector3d world=rotations[0].transpose()*(Eigen::Vector3d(p.x,p.y,p.z)-translations[0]);
        if(i>=structure.size())world.z()+=1;
        const Eigen::Vector3d raw=futureRotation*world+futureTranslation;
        future.push_back({float(raw.x()),float(raw.y()),float(raw.z()),p.intensity});
    }
    QFile futureFile(futurePath);if(!futureFile.open(QIODevice::WriteOnly))return 1;
    futureFile.write(QString("ply\nformat binary_little_endian 1.0\nelement vertex %1\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n").arg(qulonglong(future.size())).toUtf8());
    const auto futureHeader=futureFile.pos();
    futureFile.write(reinterpret_cast<const char*>(future.data()),qint64(future.size()*sizeof(TestRecord)));futureFile.close();
    const auto newTransform=loadCalibrationTransform(futurePath);if(!newTransform.applied)return 1;
    for(size_t i=0;i<structure.size();i+=113) {
        const auto &p=future[i];const auto &a=newTransform.matrix;
        const Eigen::Vector3d q(a[0]*p.x+a[1]*p.y+a[2]*p.z+a[3],a[4]*p.x+a[5]*p.y+a[6]*p.z+a[7],a[8]*p.x+a[9]*p.y+a[10]*p.z+a[11]);
        if((q-(structure[i]-Eigen::Vector3d(0,0,10))).norm()>.03)return 1;
    }
    // 延伸料區端牆比主端牆密，仍須以永久跨度選回主建築牆，不移動共同原點。
    // A denser extension wall must not replace the main end wall or shift the frozen origin.
    auto extension=future;
    for(int x=0;x<200;++x) for(int z=0;z<50;++z) {
        const Eigen::Vector3d raw=futureRotation*Eigen::Vector3d((x+.5)*.1,-13,(z+.5)*.2)+futureTranslation;
        extension.push_back({float(raw.x()),float(raw.y()),float(raw.z()),42});
    }
    const auto extensionPath=QDir(temp.path()).filePath("extended_scan.ply");QFile extensionFile(extensionPath);
    if(!extensionFile.open(QIODevice::WriteOnly)) return 1;
    extensionFile.write(QString("ply\nformat binary_little_endian 1.0\nelement vertex %1\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n").arg(qulonglong(extension.size())).toUtf8());
    extensionFile.write(reinterpret_cast<const char*>(extension.data()),qint64(extension.size()*sizeof(TestRecord)));extensionFile.close();
    const auto extensionTransform=loadCalibrationTransform(extensionPath);if(!extensionTransform.applied) return 1;
    for(size_t i=0;i<structure.size();i+=113) {
        const auto &p=extension[i];const auto &a=extensionTransform.matrix;
        const Eigen::Vector3d q(a[0]*p.x+a[1]*p.y+a[2]*p.z+a[3],a[4]*p.x+a[5]*p.y+a[6]*p.z+a[7],a[8]*p.x+a[9]*p.y+a[10]*p.z+a[11]);
        if((q-(structure[i]-Eigen::Vector3d(0,0,10))).norm()>.03) {std::cerr<<"Extension wall shifted shared coordinates: "<<q.transpose()<<" expected "<<(structure[i]-Eigen::Vector3d(0,0,10)).transpose()<<"\n";return 1;}
    }
    siteBefore.open(QIODevice::ReadOnly);if(siteBefore.readAll()!=frozen)return 1;siteBefore.close();
    if(QFile::exists(QDir(temp.path()).filePath("calibrated/a_future.ply")))return 1; // 自動註冊只寫矩陣。 / Automatic registration writes matrices only.
    QFile cache(QDir(temp.path()).filePath("calibration.json"));cache.open(QIODevice::ReadOnly);
    const auto oldCache=QJsonDocument::fromJson(cache.readAll()).object();cache.close();
    // 同名、同大小檔案的內容變更也必須偵測到，而不重用舊 SHA。 / Detect same-name/same-size content changes rather than reuse a stale SHA.
    future[0].intensity+=1;
    if(!futureFile.open(QIODevice::ReadWrite)||!futureFile.seek(futureHeader))return 1;
    futureFile.write(reinterpret_cast<const char*>(&future[0]),sizeof(TestRecord));futureFile.close();
    if(!loadCalibrationTransform(futurePath).applied)return 1;
    cache.open(QIODevice::ReadOnly);const auto newCache=QJsonDocument::fromJson(cache.readAll()).object();cache.close();
    const auto previous=oldCache.value("files").toObject(),current=newCache.value("files").toObject();
    if(previous.value("scan0.ply")!=current.value("scan0.ply") || previous.value("scan1.ply")!=current.value("scan1.ply") ||
        previous.value("a_future.ply").toObject().value("source_sha256")==current.value("a_future.ply").toObject().value("source_sha256"))return 1;
    siteBefore.open(QIODevice::ReadOnly);if(siteBefore.readAll()!=frozen)return 1;siteBefore.close();
    // 只有料堆、看不到固定結構的掃描應失敗，且不能新增有效校準紀錄。
    // A pile-only scan must fail without adding a valid calibration record.
    const auto badPath=QDir(temp.path()).filePath("pile_only.ply");QFile bad(badPath);if(!bad.open(QIODevice::WriteOnly))return 1;
    bad.write("ply\nformat binary_little_endian 1.0\nelement vertex 20\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n");
    bad.write(reinterpret_cast<const char*>(future.data()),20*sizeof(TestRecord));bad.close();
    bool failed=false;try{loadCalibrationTransform(badPath);}catch(const std::exception &){failed=true;}if(!failed)return 1;
    cache.open(QIODevice::ReadOnly);const auto rejectedCache=QJsonDocument::fromJson(cache.readAll()).object();cache.close();
    if(rejectedCache.value("files").toObject().contains("pile_only.ply"))return 1;
    // 不允許鏡射／縮放矩陣。 / Reject reflection/scale matrices.
    QFile sidecar(QDir(temp.path()).filePath("calibration.json"));
    if (!sidecar.open(QIODevice::WriteOnly|QIODevice::Truncate)) return 1;
    sidecar.write("{\"version\":1,\"files\":{\"scan0.ply\":{\"source_to_building\":[-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}}}");sidecar.close();
    bool rejected=false; try {loadCalibrationTransform(QDir(temp.path()).filePath("scan0.ply"));} catch (...) {rejected=true;}
    if (!rejected) return 1;
    std::cout<<"PASS: frozen site, future scans, 180-degree pose, changed piles, incremental records, full SHA invalidation, rejection of insufficient structure, preserved intensity/counts\n";
    return 0;
}
