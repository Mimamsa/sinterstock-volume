/** @file ply_export.cpp
 * @brief 完整處理結果的二進位 PLY 匯出。 / Binary PLY export of full processed results.
 */
#include "cloud.h"
#include <QDir>
#include <QFileInfo>
#include <QTemporaryFile>
#include <QtEndian>
#include <cstring>
#include <stdexcept>

void exportProcessedPly(const QString &path,const std::vector<Point> &measured,const std::vector<Point> &estimated,bool calibrated) {
    const QFileInfo target(path);
    if(target.exists()) throw std::runtime_error("Output already exists; choose a different filename");
    if(!target.dir().exists()) throw std::runtime_error("Output directory does not exist");
    QTemporaryFile output(target.dir().filePath(".pointcloud-export-XXXXXX"));
    if(!output.open()) throw std::runtime_error("Cannot create temporary PLY output");
    /** @brief 步驟 1：宣告完整點數與來源標記；已校準座標寫入防止二次校準的標頭。
     *         Step 1: declare full counts and provenance; calibrated coordinates carry a marker preventing double registration.
     */
    QByteArray header="ply\nformat binary_little_endian 1.0\ncomment processed measured points followed by estimated fill; meters\n";
    if(calibrated) header+="comment calibrated building frame; meters; origin at reference roof corner\n";
    header+="comment is_estimated: 0 measured, 1 synthetic; missing intensity is NaN\n";
    header+="element vertex "+QByteArray::number(qulonglong(measured.size()+estimated.size()))+"\nproperty float x\nproperty float y\nproperty float z\nproperty float intensity\nproperty uchar is_estimated\nend_header\n";
    auto write=[&](const QByteArray &data) {if(output.write(data)!=data.size()) throw std::runtime_error("PLY write failed (check disk space)");};
    write(header);
    /** @brief 步驟 2：固定大小批次寫入小端 float32，不配置整份輸出陣列。
     *         Step 2: write little-endian float32 records in bounded batches, without allocating the full output.
     * 保留量測強度；補點強度為 NaN，另以 uchar 明確區分來源。
     * Retain measured intensity; synthetic intensity is NaN and an uchar explicitly identifies provenance.
     */
    QByteArray block;block.reserve(65536*17);
    auto append=[&](const std::vector<Point> &points,unsigned char synthetic) {
        for(const auto &p:points) {
            const float values[]{p.x,p.y,p.z,synthetic?std::numeric_limits<float>::quiet_NaN():p.intensity};
            for(float value:values) {quint32 bits;std::memcpy(&bits,&value,4);bits=qToLittleEndian(bits);block.append(reinterpret_cast<const char *>(&bits),4);}
            block.append(char(synthetic));
            if(block.size()>=65536*17) {write(block);block.clear();}
        }
    };
    append(measured,0);append(estimated,1);if(!block.isEmpty()) write(block);
    if(!output.flush()) throw std::runtime_error("Cannot flush PLY output");
    /** @brief 步驟 3：成功寫完才發布檔案；rename 不覆寫既有檔，失敗自動移除暫存。
     *         Step 3: publish only after successful writing; rename refuses existing files and failures remove the temporary output.
     */
    if(!output.rename(target.absoluteFilePath())) throw std::runtime_error("Cannot publish PLY (destination exists or is not writable)");
    output.setAutoRemove(false);
}
