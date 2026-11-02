/**
 * @file cloud_test.cpp
 * @brief 原始量測保留、資料驗證與顯示上限測試。 / Tests for raw measurements, input validation, and display limits.
 */
#include "cloud.h"
#include <QTemporaryDir>
#include <QFile>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstring>

/**
 * @brief 使用暫存 PLY 驗證處理流程。 / Verify processing using temporary PLY fixtures.
 * @return 0 表示全部通過，1 表示檔案建立或檢查失敗。 / Zero if all checks pass; one on fixture or verification failure.
 */
int main() {
    /** @brief 建立含同格多點、負座標及 NaN 的測試輸入。
     *         Create a fixture with nearby distinct measurements, negative coordinates, and NaN.
     */
    QTemporaryDir dir;
    if (!dir.isValid()) return 1;
    const QString path = dir.filePath("test.ply");
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 7\nproperty float x\nproperty float y\nproperty float z\nend_header\n"
               "0 0 0\n0.01 0 0\n0.20 0 0\n0.21 0 0\n-0.01 0 0\n-0.02 0 0\nnan 0 0\n");
    file.close();
    /** @brief 不限制樣本，檢查有效點數與每個原始座標。
     *         Keep all samples and verify the valid input count and three raw measurements.
     */
    auto full = processCloud(path,0);
    if (!full.error.isEmpty() || full.originalCount != 6 || full.reducedCount != 6) return 1;
    if (full.hasIntensity || !std::isnan(full.before.front().intensity)) return 1;
    if (full.intensityP05!=0 || full.intensityP95!=0 || full.zP05!=0 || full.zP95!=0) return 1;
    std::sort(full.after.begin(),full.after.end(),[](Point a,Point b) {return a.x < b.x;});
    const double expected[]{-.02,-.01,0,.01,.20,.21};
    for (int i=0;i<6;++i) {
        if (std::abs(full.after[i].x-expected[i])>1e-6 || full.after[i].y || full.after[i].z) return 1;
    }
    /** @brief 限制顯示點數，確認完整量測統計維持相同。
     *         Cap display samples and confirm full measurement statistics remain unchanged.
     */
    auto limited = processCloud(path,2);
    if (!limited.error.isEmpty() || limited.originalCount != 6 || limited.reducedCount != 6 ||
        limited.before.size()!=2 || limited.after.size()!=2) return 1;
    /** @brief 故意省略標頭宣告的第二個點，確認截斷資料回報錯誤。
     *         Omit the declared second vertex to verify truncated input reports an error.
     */
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 2\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n");
    file.close();
    if (processCloud(path,2).error.isEmpty()) return 1;
    /** @brief 驗證三種編碼的強度讀取、NaN/Inf 強度處理、原始強度保留與共用範圍。
     *         Verify intensity decoding in three encodings, NaN/Inf handling, original intensities, and shared ranges.
     */
    const double intensities[]{10,std::numeric_limits<double>::quiet_NaN(),30,
                               std::numeric_limits<double>::infinity(),-5};
    const float xs[]{0,.01f,.2f,.21f,-.01f};
    for (int encoding=0;encoding<3;++encoding) {
        if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
        const QByteArray format = encoding == 0 ? "ascii" : encoding == 1 ? "binary_little_endian" : "binary_big_endian";
        // 強度放在 XYZ 前以驗證欄位位移。 / Place intensity before XYZ to verify property offsets.
        file.write("ply\nformat " + format + " 1.0\nelement vertex 5\nproperty double intensity\nproperty float z\nproperty float x\nproperty float y\nend_header\n");
        for (int i=0;i<5;++i) {
            float z = i==4 ? 3.f : 2.f;
            if (!encoding) {
                file.write(QByteArray::number(intensities[i],'g',17) + " " + QByteArray::number(z) + " " + QByteArray::number(xs[i],'g',9) + " 0\n");
            } else {
                auto writeBits = [&](auto bits) {
                    auto ordered = encoding==1 ? qToLittleEndian(bits) : qToBigEndian(bits);
                    file.write(reinterpret_cast<const char *>(&ordered),sizeof(ordered));
                };
                quint64 bits; std::memcpy(&bits,&intensities[i],sizeof(bits)); writeBits(bits);
                for (float value : {z,xs[i],0.f}) {
                    quint32 b; std::memcpy(&b,&value,sizeof(b)); writeBits(b);
                }
            }
        }
        file.close();
        auto colored = processCloud(path,0);
        if (!colored.error.isEmpty() || colored.originalCount!=5 || colored.reducedCount!=5 || !colored.hasIntensity ||
            colored.intensityMin!=-5 || colored.intensityMax!=30 || colored.zMin!=2 || colored.zMax!=3) return 1;
        if (std::abs(colored.intensityP05+3.5f)>1e-5 || std::abs(colored.intensityP95-28.f)>1e-5 ||
            colored.zP05!=2 || std::abs(colored.zP95-2.8f)>1e-5) return 1;
        std::sort(colored.after.begin(),colored.after.end(),[](Point a,Point b) {return a.x<b.x;});
        const float expectedIntensity[]{-5,10,std::numeric_limits<float>::quiet_NaN(),30,std::numeric_limits<float>::infinity()};
        for (int i=0;i<5;++i) if (std::isfinite(expectedIntensity[i]) ? colored.after[i].intensity!=expectedIntensity[i] : std::isfinite(colored.after[i].intensity)) return 1;
        auto capped = processCloud(path,1);
        if (capped.intensityMin!=-5 || capped.intensityMax!=30 || capped.zMin!=2 || capped.zMax!=3) return 1;
        if (capped.zP05!=colored.zP05 || capped.zP95!=colored.zP95 ||
            capped.intensityP05!=colored.intensityP05 || capped.intensityP95!=colored.intensityP95) return 1;
    }
    // 同一 10 cm 空間內仍保留不同量測與強度。 / Preserve distinct measurements/intensities even within one 10 cm cell.
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 2\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n0 0 0 10\n0.01 0 0 30\n");
    file.close();
    auto mean = processCloud(path,0);
    if (!mean.error.isEmpty() || mean.reducedCount!=2) return 1;
    std::sort(mean.after.begin(),mean.after.end(),[](Point a,Point b){return a.x<b.x;});
    if(mean.after[0].intensity!=10 || mean.after[1].intensity!=30 || std::abs(mean.after[1].x-.01f)>1e-7) return 1;
    if (mean.intensityP05!=11 || mean.intensityP95!=29) return 1;
    /** @brief 驗證離群值不會成為 P5/P95，且顯示上限不影響百分位數。
     *         Verify outliers do not become P5/P95 and display caps do not affect percentiles.
     */
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 102\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n");
    for (int i=0;i<=100;++i) {
        const int value = i==100 ? 10000 : i;
        file.write("0 0 " + QByteArray::number(value) + " " + QByteArray::number(value) + "\n");
    }
    file.write("nan 0 0 -10000\n"); file.close();
    auto robust = processCloud(path,1);
    if (!robust.error.isEmpty() || robust.originalCount!=101 || robust.before.size()!=1 ||
        robust.zP05!=5 || robust.zP95!=95 || robust.intensityP05!=5 || robust.intensityP95!=95 ||
        robust.zMax!=10000 || robust.intensityMax!=10000) return 1;
    // 單點與強度全缺值的邊界案例。 / Singleton and all-missing-intensity edge cases.
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n0 0 7 nan\n");
    file.close();
    auto singleton = processCloud(path,0);
    if (!singleton.error.isEmpty() || singleton.zP05!=7 || singleton.zP95!=7 || singleton.hasIntensity ||
        singleton.intensityP05!=0 || singleton.intensityP95!=0) return 1;
    /** @brief 建立相鄰主體、分離密集群及孤立點，驗證密度與連通性過濾。
     *         Build an adjacent main structure, detached dense cluster, and singleton to verify density and connectivity filtering.
     */
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) return 1;
    file.write("ply\nformat ascii 1.0\nelement vertex 10\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nend_header\n");
    for (float x : {-.49f,-.38f,-.27f,.01f,.12f,.23f,10.01f,10.12f,10.23f,20.f}) {
        file.write(QByteArray::number(x,'g',9)+" 0 0 42\n");
    }
    file.close();
    const auto unfiltered = processCloud(path,0);
    DenoiseOptions noise{true,.5,3,true};
    const auto cleaned = processCloud(path,0,noise);
    if (!cleaned.error.isEmpty() || cleaned.originalCount!=10 || 
        cleaned.reducedCount!=6 || cleaned.after.size()!=6 || !cleaned.denoised ||
        cleaned.before.size()!=unfiltered.before.size() || cleaned.intensityP95!=unfiltered.intensityP95) return 1;
    for (const auto &p : cleaned.after) if (p.x>=1 || p.intensity!=42) return 1;
    noise.largestComponent = false;
    const auto densityOnly = processCloud(path,0,noise);
    if (!densityOnly.error.isEmpty() || densityOnly.reducedCount!=9) return 1;
    noise.minPoints = 100;
    const auto empty = processCloud(path,0,noise);
    if (!empty.error.isEmpty() || empty.reducedCount || !empty.after.empty() || empty.before.size()!=10) return 1;
    noise.cellSize = 0;
    if (processCloud(path,0,noise).error.isEmpty()) return 1;
    std::cout << "PASS: raw coordinates, negative coordinates, invalid coordinates, preview cap, truncated input\n";
    std::cout << "PASS: intensity in ASCII/both byte orders, missing intensity, original intensities, full-input color ranges\n";
    std::cout << "PASS: interpolated P5/P95, outliers, display independence, singleton and missing values\n";
    std::cout << "PASS: sparse cells, largest component, density-only mode, empty output, invalid settings\n";
    return 0;
}
