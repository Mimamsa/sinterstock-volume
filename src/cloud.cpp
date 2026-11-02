/**
 * @file cloud.cpp
 * @brief 完整 PLY 量測清理、補面與體積處理。 / Full PLY measurement cleaning, reconstruction and volume processing.
 */
#include "cloud.h"
#include "config_paths.h"
#include "calibration_transform.h"
#include <QtEndian>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QCoreApplication>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace {
/** @brief PLY 標量欄位：名稱、型別、記錄內位移及位元組數。 / PLY scalar field: name, type, record offset, and byte size. */
struct Property { std::string name, type; size_t offset, size; };
/**
 * @brief 將 PLY 標量型別別名映射為位元組數。 / Map a PLY scalar type alias to its byte width.
 * @param t PLY 型別名稱。 / PLY type name.
 * @return 1、2、4 或 8 位元組。 / Width of 1, 2, 4, or 8 bytes.
 * @throws std::runtime_error 遇到未支援型別。 / If the type is unsupported.
 */
size_t typeSize(const std::string &t) {
    if (t == "char" || t == "uchar" || t == "int8" || t == "uint8") return 1;
    if (t == "short" || t == "ushort" || t == "int16" || t == "uint16") return 2;
    if (t == "int" || t == "uint" || t == "float" || t == "int32" || t == "uint32" || t == "float32") return 4;
    if (t == "double" || t == "float64") return 8;
    throw std::runtime_error("Unsupported PLY scalar type: " + t);
}
/**
 * @brief 以指定端序解碼一個二進位標量。 / Decode one binary scalar with the specified byte order.
 * @param p 指向欄位第一個位元組，至少有 v.size 個可讀位元組。 / Field bytes; at least v.size bytes must be readable.
 * @param v 由已驗證標頭取得的型別描述。 / Type descriptor from the validated header.
 * @param little true 為小端，false 為大端。 / True for little endian, false for big endian.
 * @return 轉成 double 的欄位值。 / Field value converted to double.
 * @details 先轉換整數位元順序，再用 memcpy 還原浮點位元，避免未對齊讀取與別名問題。
 *          Normalize integer byte order, then restore floating-point bits with memcpy to avoid unaligned access and aliasing issues.
 */
double scalar(const char *p, const Property &v, bool little) {
    quint64 bits = 0;
    if (v.size == 1) bits = static_cast<unsigned char>(*p);
    if (v.size == 2) bits = little ? qFromLittleEndian<quint16>(p) : qFromBigEndian<quint16>(p);
    if (v.size == 4) bits = little ? qFromLittleEndian<quint32>(p) : qFromBigEndian<quint32>(p);
    if (v.size == 8) bits = little ? qFromLittleEndian<quint64>(p) : qFromBigEndian<quint64>(p);
    if (v.type == "float" || v.type == "float32") { quint32 b = bits; float f; std::memcpy(&f, &b, 4); return f; }
    if (v.type == "double" || v.type == "float64") { double d; std::memcpy(&d, &bits, 8); return d; }
    if (v.type == "char" || v.type == "int8") return static_cast<qint8>(bits);
    if (v.type == "short" || v.type == "int16") return static_cast<qint16>(bits);
    if (v.type == "int" || v.type == "int32") return static_cast<qint32>(bits);
    return static_cast<double>(bits);
}
/** @brief 密度／連通性用三軸整數格索引。 / Integer 3D cell indices for density/connectivity. */
using Key = std::array<qint64, 3>;
/**
 * @brief 依輸入索引均勻挑選顯示點。 / Select display points at evenly spaced input indices.
 * @param p 不修改的來源陣列。 / Source array, left unchanged.
 * @param limit 最大點數，0 表示全部。 / Maximum point count; 0 means all.
 * @return 獨立的顯示樣本；空輸入回傳空陣列。 / Independent display sample; empty input yields an empty array.
 * @note 這不是體素降取樣或隨機取樣，也不保證空間分布均勻。
 *       This is neither voxel nor random sampling and does not guarantee uniform spatial coverage.
 */
std::vector<Point> sample(const std::vector<Point> &p, size_t limit) {
    size_t n = limit ? std::min(limit, p.size()) : p.size();
    std::vector<Point> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) out.push_back(p[i * p.size() / n]);
    return out;
}

/**
 * @brief 原地計算有限標量的 P5/P95，不配置另一份完整標量陣列。
 *        Compute finite scalar P5/P95 in place without allocating another full scalar array.
 * @param begin 有效範圍起點。 / Start of the valid range.
 * @param end 有效範圍終點（不包含）。 / Exclusive end of the valid range.
 * @param field 要比較的 Point 成員；範圍內該欄位必須皆有限。
 *              Point member to compare; this field must be finite throughout the range.
 * @return P5/P95；空範圍回傳 0/0，單點回傳同一數值。
 *         P5/P95; empty ranges yield 0/0 and singletons yield the same value twice.
 * @details 使用排序位置 r=(N-1)*p，在 floor(r) 與 ceil(r) 間線性內插。
 *          Use rank r=(N-1)*p and interpolate linearly between floor(r) and ceil(r).
 *          nth_element 找到下界，右側最小值為下一順序統計量；每個百分位平均 O(N)。
 *          nth_element finds the lower rank and the suffix minimum gives the next order statistic; each percentile is average O(N).
 * @note 會改變元素順序，必須在保留原始顯示樣本後呼叫。
 *       Reorders points; call only after preserving the original display sample.
 */
std::pair<float,float> percentileRange(std::vector<Point>::iterator begin,
                                     std::vector<Point>::iterator end, float Point::*field) {
    const auto count = end - begin;
    if (!count) return {0,0};
    auto compare = [field](const Point &a, const Point &b) { return a.*field < b.*field; };
    auto quantile = [&](double fraction) {
        const double rank = double(count-1) * fraction;
        const auto index = static_cast<std::ptrdiff_t>(std::floor(rank));
        auto selected = begin + index;
        std::nth_element(begin,selected,end,compare);
        const double lower = (*selected).*field;
        if (rank == double(index)) return float(lower);
        const double upper = (*std::min_element(selected+1,end,compare)).*field;
        return float(lower + (upper-lower)*(rank-double(index)));
    };
    const float lower = quantile(.05);
    const float upper = quantile(.95);
    return {lower,upper};
}

/** @brief 密度格整數索引雜湊。 / Hash for integer density-cell indices. */
struct KeyHash {
    size_t operator()(const Key &k) const {
        size_t h = 0;
        for (auto v : k) h ^= std::hash<qint64>{}(v) + size_t(0x9e3779b9) + (h<<6) + (h>>2);
        return h;
    }
};

/**
 * @brief 移除稀疏網格與分離的小點群，原地保留主要結構。
 *        Remove sparse cells and detached clusters, retaining the primary structure in place.
 * @param points 完整校準後原始量測；會刪除不符合條件的點。
 *               Full calibrated raw measurements; points failing the criteria are removed.
 * @param options 密度網格與連通區設定。 / Density-grid and component settings.
 * @note 這是幾何啟發式，不能辨認建築物；與主體相連的雜點可能保留，分離的真實構件可能移除。
 *       This geometric heuristic does not recognize buildings: attached noise may remain and detached real structures may be removed.
 */
void denoiseCloud(std::vector<Point> &points, const DenoiseOptions &options) {
    if (!std::isfinite(options.cellSize) || options.cellSize < .1 || !options.minPoints)
        throw std::runtime_error("Denoising requires cell size >= 0.1 m and minimum points >= 1");
    auto cellKey = [&](const Point &p) -> Key {
        return {qint64(std::floor(double(p.x)/options.cellSize)),
                qint64(std::floor(double(p.y)/options.cellSize)),
                qint64(std::floor(double(p.z)/options.cellSize))};
    };
    struct Cell { size_t count = 0, label = 0; };
    std::unordered_map<Key,Cell,KeyHash> cells;
    /** @brief 步驟 A：統計每個粗網格內的質心數，不配置逐點鄰居清單。
     *         Step A: count centroids in coarse cells without allocating per-point neighbor lists.
     */
    for (const auto &p : points) ++cells[cellKey(p)].count;
    size_t largestLabel = 0, largestSize = 0, nextLabel = 0;
    if (options.largestComponent) {
        std::vector<Key> pending;
        /** @brief 步驟 B：只讓密度足夠的格子參與 26 鄰接搜尋；累計點數選出主連通區。
         *         Step B: traverse 26-connected sufficiently dense cells and select the component with the most points.
         * 面、邊、角相接皆算鄰接；非遞迴搜尋避免大型建築造成呼叫堆疊溢位。
         * Face, edge, and corner contacts all connect; iterative traversal avoids call-stack overflow for large buildings.
         */
        for (auto &entry : cells) {
            if (entry.second.count < options.minPoints || entry.second.label) continue;
            const size_t label = ++nextLabel;
            entry.second.label = label; pending.push_back(entry.first);
            size_t size = 0;
            while (!pending.empty()) {
                const Key k = pending.back(); pending.pop_back();
                size += cells.find(k)->second.count;
                for (int dx=-1;dx<=1;++dx) for (int dy=-1;dy<=1;++dy) for (int dz=-1;dz<=1;++dz) {
                    if (!dx && !dy && !dz) continue;
                    auto found = cells.find({k[0]+dx,k[1]+dy,k[2]+dz});
                    if (found == cells.end() || found->second.count < options.minPoints || found->second.label) continue;
                    found->second.label = label; pending.push_back(found->first);
                }
            }
            if (size > largestSize) { largestSize = size; largestLabel = label; }
        }
    }
    /** @brief 步驟 C：原地壓縮保留點，XYZ 與原始強度均不修改。
     *         Step C: compact retained points in place without changing XYZ or original intensities.
     * 若沒有格子通過條件，結果為空，不靜默恢復原始雜點。
     * If no cells pass, return an empty result rather than silently restoring noise.
     */
    points.erase(std::remove_if(points.begin(),points.end(),[&](const Point &p) {
        const auto &cell = cells.find(cellKey(p))->second;
        return cell.count < options.minPoints || (options.largestComponent && cell.label != largestLabel);
    }),points.end());
}
}

size_t removeFloatingNoise(std::vector<Point> &points, const DenoiseOptions &options) {
    if(!options.enabled) return 0;
    const size_t before=points.size();
    denoiseCloud(points,options);
    return before-points.size();
}

Cloud processCloud(const QString &path, size_t displayLimit, const DenoiseOptions &denoise, const FillOptions &fill, bool applyCalibration, const StructureOptions &structures, const VolumeOptions &volume, const QString &exportPath) {
    /** @brief 步驟 1：建立結果並開啟串流。 / Step 1: initialize the result and open the stream.
     * 二進位模式可避免平台換行轉換；ASCII 標頭仍可用 getline 解析。
     * Binary mode avoids platform newline translation; getline can still parse the ASCII header.
     */
    Cloud result;
    try {
        const auto calibration=applyCalibration ? loadCalibrationTransform(path) : CalibrationTransform{};
        if (calibration.applied) result.coordinateFrame=calibration.frame;
        std::ifstream in(path.toStdString(), std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open PLY");
        /** @brief 步驟 2：確認 PLY 識別行。 / Step 2: verify the PLY signature line.
         * 接受 LF 與 CRLF；getline 移除 LF，但 CRLF 可能留下 CR。
         * Accept LF and CRLF; getline removes LF but may leave CR.
         */
        std::string line, format, element;
        std::getline(in, line);
        if (line != "ply" && line != "ply\r") throw std::runtime_error("Not a PLY file");
        /** @brief 步驟 3：準備標頭解析狀態。 / Step 3: prepare header parsing state.
         * props 保留 vertex 屬性順序；count 是宣告點數，stride 是單筆二進位記錄長度。
         * props preserves vertex property order; count is the declared vertex count and stride is the binary record size.
         * ended 確認標頭結束，vertexSeen 防止重複宣告 vertex。
         * ended tracks the header terminator; vertexSeen rejects duplicate vertex declarations.
         */
        std::vector<Property> props;
        size_t count = 0, stride = 0;
        bool ended = false, vertexSeen = false;
        /** @brief 步驟 4：逐行讀取標頭，直到 end_header。 / Step 4: read header lines through end_header.
         * format 決定資料解碼方式；未處理的標籤（例如 comment）會自然略過。
         * format selects the decoder; unhandled tags such as comment are naturally ignored.
         */
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string tag; row >> tag;
            if (tag == "end_header") { ended = true; break; }
            if (tag == "format") row >> format;
            /** @brief 步驟 4a：追蹤目前 element 與其數量。 / Step 4a: track the current element and its count.
             * 本讀取器不跳過 vertex 前的資料區，因此要求 vertex 為第一個非空 element。
             * This reader does not skip preceding data sections, so vertex must be the first non-empty element.
             * vertex 後的 face 等 element 不需讀取；本程式只使用點座標。
             * Later elements such as faces need not be read because only point coordinates are used.
             */
            if (tag == "element") {
                size_t n; if (!(row >> element >> n)) throw std::runtime_error("Invalid element header");
                if (element == "vertex") {
                    if (vertexSeen) throw std::runtime_error("Duplicate vertex element");
                    count = n; vertexSeen = true;
                } else if (!vertexSeen && n) throw std::runtime_error("Vertex must be the first non-empty element");
            }
            /** @brief 步驟 4b：建立 vertex 欄位配置。 / Step 4b: build the vertex field layout.
             * 記錄所有標量欄位，包含 intensity/RGB，才能正確計算 XYZ 的位移與 stride。
             * Record every scalar field, including intensity/RGB, to compute correct XYZ offsets and stride.
             * list 長度可變，無法使用固定 stride，因此明確拒絕。
             * Lists have variable length and cannot use a fixed stride, so they are explicitly rejected.
             */
            if (tag == "property" && element == "vertex") {
                std::string type, name; row >> type >> name;
                if (type == "list") throw std::runtime_error("List properties in vertices are unsupported");
                size_t size = typeSize(type);
                props.push_back({name, type, stride, size}); stride += size;
            }
        }
        /** @brief 步驟 5：驗證標頭完整性與編碼。 / Step 5: validate header completeness and encoding.
         * 必須有非空 vertex、end_header，以及支援的 ASCII／二進位格式。
         * Require a non-empty vertex element, end_header, and a supported ASCII/binary format.
         */
        if (!ended || !vertexSeen || !count) throw std::runtime_error("Missing or empty vertex element");
        if (format != "ascii" && format != "binary_little_endian" && format != "binary_big_endian")
            throw std::runtime_error("Unsupported PLY format");
        /** @brief 步驟 6：依名稱找出 XYZ 欄位索引。 / Step 6: locate XYZ property indices by name.
         * 不假設 XYZ 是前三欄或固定排列；-1 表示必要欄位尚未找到。
         * XYZ need not be the first three fields or in a fixed order; -1 marks a missing required field.
         */
        std::array<int, 3> xyz{-1, -1, -1};
        int intensityIndex = -1;
        for (size_t i = 0; i < props.size(); ++i) {
            if (props[i].name == "x") xyz[0] = i;
            if (props[i].name == "y") xyz[1] = i;
            if (props[i].name == "z") xyz[2] = i;
            if (props[i].name == "intensity") intensityIndex = int(i);
        }
        if (*std::min_element(xyz.begin(), xyz.end()) < 0) throw std::runtime_error("Missing x/y/z properties");
        /** @brief 步驟 7：預留完整 XYZ 工作陣列。 / Step 7: reserve the full XYZ working array.
         * 依標頭點數一次 reserve，減少讀取時擴容；無效點會使實際大小小於宣告數量。
         * Reserve once from the header count to reduce reallocations; invalid points may reduce the actual size.
         */
        std::vector<Point> points; points.reserve(count);
        /** @brief 步驟 8：共用的座標驗證與寫入流程。 / Step 8: shared coordinate validation and append path.
         * 任一座標為 NaN/Inf 就略過整個點；超出範圍則報錯，避免體素整數索引溢位。
         * Skip the entire point if any coordinate is NaN/Inf; reject excessive magnitudes to protect integer voxel indexing.
         * 驗證後轉為 float32 XYZ；更新從零開始的邊界，讓初始取景包含原點。
         * Convert accepted coordinates to float32 XYZ and update zero-initialized bounds so initial framing includes the origin.
         */
        /** @brief 額外讀取 intensity；缺值或超出 float32 範圍時保留 NaN，並以完整資料更新色階範圍。
         *         Also read intensity; missing/out-of-float32-range values remain NaN, and full input data determines color ranges.
         */
        auto append = [&](const std::array<double, 3> &v, double intensity) {
            for (double x : v) {
                if (!std::isfinite(x)) return; // discard NaN/Inf points
                if (std::abs(x) > 1e15) throw std::runtime_error("Coordinate exceeds supported range");
            }
            Point p{float(v[0]), float(v[1]), float(v[2])};
            /** @brief 在統計、清理及補點前套用同一建築座標系，避免重複轉換輸出檔。
             *         Apply the shared building frame before statistics, sampling, and filling; exported files have no source sidecar.
             */
            if (calibration.applied) {
                const auto &m=calibration.matrix;
                p.x=float(m[0]*v[0]+m[1]*v[1]+m[2]*v[2]+m[3]);
                p.y=float(m[4]*v[0]+m[5]*v[1]+m[6]*v[2]+m[7]);
                p.z=float(m[8]*v[0]+m[9]*v[1]+m[10]*v[2]+m[11]);
            }
            if (std::isfinite(intensity) && std::abs(intensity) <= std::numeric_limits<float>::max()) {
                p.intensity = float(intensity);
                if (!result.hasIntensity) result.intensityMin = result.intensityMax = p.intensity;
                result.intensityMin = std::min(result.intensityMin, p.intensity);
                result.intensityMax = std::max(result.intensityMax, p.intensity);
                result.hasIntensity = true;
            }
            if (points.empty()) result.zMin = result.zMax = p.z;
            result.zMin = std::min(result.zMin, p.z);
            result.zMax = std::max(result.zMax, p.z);
            points.push_back(p);
            result.low = {std::min(result.low.x, p.x), std::min(result.low.y, p.y), std::min(result.low.z, p.z)};
            result.high = {std::max(result.high.x, p.x), std::max(result.high.y, p.y), std::max(result.high.z, p.z)};
        };
        /** @brief 步驟 9：依資料編碼選擇讀取分支。 / Step 9: select the input path by encoding. */
        if (format == "ascii") {
            /** @brief 步驟 9a：ASCII 每行解析一筆 vertex。 / Step 9a: parse one ASCII vertex per line.
             * 依屬性數讀完所有欄位，確認 token 完整轉成數字，將 XYZ 放入 v，另外保存 intensity。
             * Consume all declared fields, require complete numeric token conversion, copy XYZ into v, and retain intensity separately.
             * 缺行、缺欄或數值格式錯誤會失敗；合法 NaN/Inf 交由 append 排除。
             * Missing lines/fields or malformed numbers fail; valid NaN/Inf tokens are filtered by append.
             */
            for (size_t i = 0; i < count; ++i) {
                if (!std::getline(in, line)) throw std::runtime_error("Truncated ASCII PLY");
                std::istringstream row(line); std::array<double, 3> v{};
                double intensity = std::numeric_limits<double>::quiet_NaN();
                for (size_t j = 0; j < props.size(); ++j) {
                    std::string token; if (!(row >> token)) throw std::runtime_error("Invalid ASCII vertex");
                    size_t used = 0; double value = std::stod(token, &used);
                    if (used != token.size()) throw std::runtime_error("Invalid scalar");
                    for (int a = 0; a < 3; ++a) if (int(j) == xyz[a]) v[a] = value;
                    if (int(j) == intensityIndex) intensity = value;
                }
                append(v, intensity);
            }
        } else {
            /** @brief 步驟 9b：二進位資料分批讀取。 / Step 9b: read binary records in batches.
             * 每批最多 65,536 筆，避免為整個原始檔配置位元組緩衝；尾批可較小。
             * Read at most 65,536 records per batch instead of buffering the entire raw file; the final batch may be smaller.
             */
            constexpr size_t batch = 65536;
            std::vector<char> buffer(batch * stride);
            for (size_t i = 0; i < count;) {
                size_t n = std::min(batch, count - i), bytes = n * stride;
                if (!in.read(buffer.data(), bytes)) throw std::runtime_error("Truncated binary PLY");
                /** @brief 步驟 9c：依記錄起點與欄位位移解碼 XYZ 與強度。 / Step 9c: decode XYZ and intensity using record starts and field offsets.
                 * j * stride 找到第 j 筆，p.offset 找到軸欄位，再由 scalar 處理型別與端序。
                 * j * stride locates record j, p.offset locates the axis field, and scalar handles its type and byte order.
                 */
                for (size_t j = 0; j < n; ++j) {
                    std::array<double, 3> v{};
                    for (int a = 0; a < 3; ++a) {
                        const auto &p = props[xyz[a]];
                        v[a] = scalar(buffer.data() + j * stride + p.offset, p, format == "binary_little_endian");
                    }
                    double intensity = std::numeric_limits<double>::quiet_NaN();
                    if (intensityIndex >= 0) {
                        const auto &p = props[intensityIndex];
                        intensity = scalar(buffer.data() + j * stride + p.offset, p, format == "binary_little_endian");
                    }
                    append(v, intensity);
                }
                i += n;
            }
        }
        /** @brief 步驟 10：確認有有效點，並在排序前保存原始顯示樣本。 / Step 10: require valid points and save the pre-sort display sample.
         * originalCount 是排除無效座標後的點數；displayLimit 不會改變完整工作陣列。
         * originalCount excludes invalid coordinates; displayLimit does not change the full working array.
         */
        if (points.empty()) throw std::runtime_error("No finite points");
        result.originalCount = points.size();
        result.before = sample(points, displayLimit);
        /** @brief 步驟 10a：從完整有效輸入計算自動色階的 P5/P95。
         *         Step 10a: compute automatic color P5/P95 from the full valid input.
         * 高度包含所有有效幾何點；強度先分隔有限值，只以有效強度計算。
         * Height includes all valid geometry; partition finite intensities first and use only that valid range.
         * partition 與百分位選擇只重排資料，不刪點、不改座標；與顯示上限無關。
         * Partition and quantile selection only reorder data, preserving all points/coordinates independently of the display cap.
         */
        std::tie(result.zP05,result.zP95) = percentileRange(points.begin(),points.end(),&Point::z);
        if (result.hasIntensity) {
            auto validEnd = std::partition(points.begin(),points.end(),[](const Point &p) { return std::isfinite(p.intensity); });
            std::tie(result.intensityP05,result.intensityP95) = percentileRange(points.begin(),validEnd,&Point::intensity);
        }
        /** @brief 校準後保留每一筆有效實測，不做三維體素質心化。
         *         Keep every calibrated finite measurement; never replace measurements with 3D voxel centroids.
         * 完整工作點雲依序去雜點、移柱、移除建築結構、清理懸浮點、補洞，再以 10 cm XY 格積分。
         * Clean noise/columns/shell/floating points, fill surfaces, then integrate on a 10 cm XY grid.
         * 顯示抽樣與著色色階不會改變計算使用的點數或座標。
         * Display sampling and color ranges never change computation coordinates or point counts.
         */
        // 柱體模式延後清理，避免密度／連通區篩選先刪掉細柱。
        // Defer cleanup in column mode so density/connectivity filtering cannot erase thin columns first.
        const bool detectBeforeDenoise=structures.enabled && structures.columns.enabled;
        if (denoise.enabled && !detectBeforeDenoise) { denoiseCloud(points,denoise); result.denoised = true; }
        /** @brief 步驟 13a-2：補洞前移除建築外殼、上層結構及懸浮殘點。
         *         Step 13a-2: remove the building shell, overhead structures, and floating remnants before filling.
         * 左側原始樣本不變；僅將清理後量測用作補洞支點，補點不回饋為量測。
         * Left samples stay unchanged; only cleaned measurements become fill supports, never synthetic points.
         */
        StructureOptions effectiveStructures=structures;
        if (structures.enabled) {
            if (!calibration.applied) throw std::runtime_error("Structure removal requires calibrated building coordinates; disable structure removal for raw scans");
            QFile manifest(QFileInfo(path).dir().filePath("calibration.json"));
            QJsonObject record;
            if (manifest.open(QIODevice::ReadOnly)) record=QJsonDocument::fromJson(manifest.readAll()).object();
            const QString model=resolveSiteModel(QFileInfo(path).dir(),record.value("site_model").toString());
            QFile file(model);
            if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error("Structure removal requires site_model.json");
            const auto object=QJsonDocument::fromJson(file.readAll()).object();
            const auto expected=record.value("site_id").toString();
            if (!expected.isEmpty() && expected!=object.value("id").toString()) throw std::runtime_error("Structure model differs from calibration site");
            std::vector<std::array<double,4>> planes;
            for (const auto &value : object.value("planes").toArray()) {
                const auto a=value.toArray(); if(a.size()!=4) throw std::runtime_error("Invalid structure plane");
                std::array<double,4> plane{};
                for(int i=0;i<4;++i) {if(!a[i].isDouble()) throw std::runtime_error("Invalid structure plane coefficient");plane[i]=a[i].toDouble();}
                planes.push_back(plane);
            }
            /** @brief 載入與校準模型獨立的真實料場輪廓，避免短牆平面截掉兩端料地。
             *         Load the actual site footprint independently of calibration, preserving end wings beyond fitted short walls.
             * 邊界以固定建築座標保存，未提供時仍支援既有四牆模型；模型識別碼須一致。
             * Persist the boundary in the fixed frame; legacy sites still use four walls. Site IDs must match.
             */
            const QString regionPath=QFileInfo(model).dir().filePath("processing_region.json");
            QFile regionFile(regionPath);
            if(regionFile.exists()) {
                if(!regionFile.open(QIODevice::ReadOnly)) throw std::runtime_error("Cannot read processing_region.json");
                const auto region=QJsonDocument::fromJson(regionFile.readAll()).object();
                if(region.value("version").toInt()!=1 || region.value("site_id").toString()!=object.value("id").toString())
                    throw std::runtime_error("Invalid processing region version/site ID");
                effectiveStructures.footprint.clear();
                for(const auto &vertex:region.value("footprint_xy_m").toArray()) {
                    const auto a=vertex.toArray();
                    if(a.size()!=2 || !a[0].isDouble() || !a[1].isDouble()) throw std::runtime_error("Invalid processing region vertex");
                    effectiveStructures.footprint.push_back({a[0].toDouble(),a[1].toDouble()});
                }
                if(effectiveStructures.footprint.empty()) throw std::runtime_error("Empty processing footprint");
                validateFootprint(effectiveStructures.footprint);
                result.footprintSource=regionPath;
            }
            if(volume.enabled && effectiveStructures.footprint.empty())
                throw std::runtime_error("Volume requires an explicit site footprint; load processing_region.json (no rectangular fallback)");
            /** @brief 步驟 13a-2a：在密度／連通區清理與一般結構裁切前提議柱體，保留可檢查的原量測證據。
             *         Step 13a-2a: propose columns before denoising or shell cropping, retaining measured evidence for review.
             * 自動提議保持待確認；只有固定場域審核檔中的 confirmed 才移除柱點。
             * Automatic proposals remain pending; only confirmed fixed-site semantics remove points.
             */
            if(structures.columns.enabled) {
                result.columns.modelPath=QFileInfo(model).dir().filePath("column_semantics.json");result.columns.siteId=object.value("id").toString();
                result.columns.items=loadColumnSemantics(result.columns.modelPath,result.columns.siteId);
                // 自動地面尚未估計時用語意參考地面；手動地面立即用於離地梁柱清理。
                // Before automatic ground estimation use the semantic reference; manual ground applies immediately.
                if(!volume.automaticGround) setColumnGround(result.columns.items,volume.ground);
                result.columns.finalized=columnSemanticsFinalized(result.columns.modelPath,result.columns.siteId);
                int nextId=0;for(const auto &c:result.columns.items) nextId=std::max(nextId,c.id);
                const int reviewedMaximumId=nextId;
                auto proposed=result.columns.finalized?std::vector<ColumnSemantic>{}:detectColumnCandidates(points,structures.columns);
                for(auto &c:proposed) {
                    const auto b=columnBounds(c);
                    if(!insideFootprint(effectiveStructures.footprint,(b[0]+b[1])*.5,(b[2]+b[3])*.5)) continue;
                    bool matched=false;
                    for(auto &existing:result.columns.items) {auto a=columnBounds(existing);
                        if(!existing.evidenceFootprint.empty()) {auto measured=existing;measured.footprint=existing.evidenceFootprint;a=columnBounds(measured);}
                        if(std::hypot((a[0]+a[1]-b[0]-b[1])*.5,(a[2]+a[3]-b[2]-b[3])*.5)<std::max(1.,structures.columns.xyCell*3)) {
                            if(existing.id>reviewedMaximumId) {
                                const double xmin=std::min(a[0],b[0]),xmax=std::max(a[1],b[1]),ymin=std::min(a[2],b[2]),ymax=std::max(a[3],b[3]);
                                if(xmax-xmin>structures.columns.maxWidth||ymax-ymin>structures.columns.maxWidth) continue;
                                existing.footprint=columnSection("rectangle",xmin,xmax,ymin,ymax);
                                existing.observedBottom=std::min(existing.observedBottom,c.observedBottom);existing.observedTop=std::max(existing.observedTop,c.observedTop);
                                existing.evidencePoints+=c.evidencePoints;existing.verticalCoverage=std::min(existing.verticalCoverage,c.verticalCoverage);existing.top=existing.observedTop;
                            } else {existing.observedBottom=c.observedBottom;existing.observedTop=c.observedTop;existing.evidencePoints=c.evidencePoints;existing.verticalCoverage=c.verticalCoverage;}
                            matched=true;break;
                        }
                    }
                    if(!matched) {c.id=++nextId;result.columns.items.push_back(c);}
                }
                validateColumns(result.columns.items);
                result.columnEvidence=sampleColumnEvidence(points,result.columns.items);
                /** @brief 柱體候選與證據已從保留完整結構的原始量測建立，現在才清理料堆工作點雲。
                 *         Candidates/evidence now come from structure-preserving raw measurements; clean the material working cloud only afterwards.
                 * 不複製整份點雲；證據採有限抽樣，且不使用補點或右側顯示抽樣。
                 * No full-cloud copy is retained; evidence is capped and never uses synthetic points or right-pane display samples.
                 */
                if(denoise.enabled) {denoiseCloud(points,denoise);result.denoised=true;}
                result.columns.removedMeasured=removeColumnPoints(points,result.columns.items);
            }
            result.structureRemoved=removeBuildingStructures(points,planes,effectiveStructures);
            /** @brief 保存建築座標屋頂平面供原始視圖顯示裁切，與計算點陣列獨立。
             *         Retain building-frame roof planes for original-view clipping independently of computation arrays.
             */
            if(planes.size()==6) {for(int roof=0;roof<2;++roof) {const auto &a=planes[roof+4];const double norm=std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
                for(int j=0;j<4;++j) result.roofPlanes[roof][j]=a[j]/(a[2]>0?norm:-norm);}
                result.hasRoofPlanes=true;
            }
            /** @brief 結構移除會切斷原本由牆頂連接的懸浮點群，故在此重新計算密度與 3D 連通區。
             *         Shell removal disconnects clusters previously linked by walls/roof; recompute density and 3D connectivity here.
             * 僅使用量測，先刪除懸浮殘點再建補洞 TIN，防止補出空中的假表面。
             * Use measurements only and remove detached remnants before TIN construction to avoid inventing airborne surfaces.
             */
            result.floatingRemoved=removeFloatingNoise(points,structures.floating);
            result.floatingFiltered=structures.floating.enabled;
            result.structuresFiltered=true;
            if(!effectiveStructures.footprint.empty() && planes.size()==6) {
                // 主矩形外的保留量測作為兩端料地復原的統計，不新增或融合量測。
                // Count retained measurements beyond the main rectangle; never invent or fuse measurements.
                for(const auto &p:points) {
                    const auto &a=planes[2],&b=planes[3];
                    const double da=a[0]*p.x+a[1]*p.y+a[2]*p.z+a[3],db=b[0]*p.x+b[1]*p.y+b[2]*p.z+b[3];
                    if(da<0) ++result.restoredEnds[0];
                    if(db>0) ++result.restoredEnds[1];
                    if(da<0 || db>0) ++result.restoredFootprintPoints;
                }
            }
        }
        result.reducedCount = points.size();
        if (!points.empty()) {
            result.cleanLow = result.cleanHigh = points.front();
            for (const auto &p : points) {
                result.cleanLow.x=std::min(result.cleanLow.x,p.x); result.cleanHigh.x=std::max(result.cleanHigh.x,p.x);
                result.cleanLow.y=std::min(result.cleanLow.y,p.y); result.cleanHigh.y=std::max(result.cleanHigh.y,p.y);
                result.cleanLow.z=std::min(result.cleanLow.z,p.z); result.cleanHigh.z=std::max(result.cleanHigh.z,p.z);
            }
        }
        /** @brief 步驟 13b：僅用完整清理後量測資料補洞，合成點獨立儲存並共用右側顯示預算。
         *         Step 13b: fill from full cleaned measurements only; keep synthetic points separate within the shared right-pane budget.
         */
        FillOptions effectiveFill=fill;
        if(structures.enabled && !effectiveStructures.footprint.empty()) {
            effectiveFill.footprint=effectiveStructures.footprint;
            effectiveFill.footprintMargin=structures.wallMargin;
        }
        validateFootprint(effectiveFill.footprint);
        if(!std::isfinite(effectiveFill.footprintMargin) || effectiveFill.footprintMargin<0) throw std::runtime_error("Invalid footprint setback");
        if (fill.enabled && fill.fullDomain && !points.empty()) {
            /** @brief 全場域以清理後量測邊界加一格邊距建立，不用原始離群點撐大網格。
             *         Derive full-domain bounds from cleaned measurements with a one-cell margin, excluding removed outliers.
             */
            effectiveFill.xMin=double(result.cleanLow.x)-fill.spacing;
            effectiveFill.xMax=double(result.cleanHigh.x)+fill.spacing;
            effectiveFill.yMin=double(result.cleanLow.y)-fill.spacing;
            effectiveFill.yMax=double(result.cleanHigh.y)+fill.spacing;
            effectiveFill.zMin=double(result.cleanLow.z)-fill.spacing;
            effectiveFill.zMax=double(result.cleanHigh.z)+fill.spacing;
            // 全場域必要時略放大網格，保持既有資源上限。 / Enlarge the grid spacing if needed to respect the resource cap.
            const double width=effectiveFill.xMax-effectiveFill.xMin,height=effectiveFill.yMax-effectiveFill.yMin;
            effectiveFill.spacing=std::max(fill.spacing,std::sqrt(width*height/1800000.));
            // 與世界原點對齊，避免最小量測點恰好變成格界而使規則掃描列因浮點誤差交錯併格。
            // Align to the world origin so the first measured point does not become a boundary that aliases regular scan rows.
            const double step=effectiveFill.spacing;
            effectiveFill.xMin=(std::floor(double(result.cleanLow.x)/step)-1)*step;
            effectiveFill.xMax=(std::ceil(double(result.cleanHigh.x)/step)+1)*step;
            effectiveFill.yMin=(std::floor(double(result.cleanLow.y)/step)-1)*step;
            effectiveFill.yMax=(std::ceil(double(result.cleanHigh.y)/step)+1)*step;
        }
        result.appliedFill=effectiveFill;
        std::vector<Point> generated;
        if (!points.empty() && effectiveFill.enabled) {
            if (effectiveFill.method==FillMethod::LayeredTin) {
                // 自動放大輸出格距時，支點間距也不得小於輸出格距。 / Adaptive support spacing must not be smaller than the output spacing.
                effectiveFill.supportSpacing=std::max(effectiveFill.supportSpacing,effectiveFill.spacing);
                result.appliedFill=effectiveFill;
                generated=reconstructStockpileSurface(points,effectiveFill,result.fillReport);
                result.filledHoles=result.fillReport.patches;
            } else generated=fillStockpileHoles(points,effectiveFill,result.filledHoles);
        }
        // 所有補洞方法都須遵守凹邊界；禁止把料場凹角或外部當作缺口補平。
        // Every fill method respects the concave perimeter; exterior recesses are not holes to fill.
        generated.erase(std::remove_if(generated.begin(),generated.end(),[&](const Point &p) {
            return !insideFootprint(effectiveFill.footprint,p.x,p.y,effectiveFill.footprintMargin);
        }),generated.end());
        result.reducedCount = points.size();
        result.filledCount = generated.size();
        /** @brief 體積前末段清理：以實測鄰格穩健平面刪除高／低離群點，補點不作支點。
         *         Final pre-volume cleanup removes high/low residual outliers using measured-only robust neighborhoods.
         * 保留刪除點供固定同一地面的前後體積比較；正式地面只從過濾後實測估計。
         * Retain rejects for a same-ground before/after comparison; fit the official ground only from filtered measurements.
         */
        std::vector<Point> surfaceRejectedMeasured,surfaceRejectedEstimated;
        if(structures.enabled) result.surfaceFilter=filterSurfaceOutliers(points,generated,structures.surfaceFilter,surfaceRejectedMeasured,surfaceRejectedEstimated);
        result.reducedCount=points.size();result.filledCount=generated.size();
        /** @brief 步驟 13c：在選取顯示樣本前估地面、分割邊界與計算 10 cm 柵格體積。
         *         Step 13c: estimate ground, segment boundaries, and integrate 10 cm raster volumes before display sampling.
         * 此處使用完整量測與補點，顯示點數上限不影響體積；體積錯誤獨立顯示，不阻擋點雲。
         * Full measurements and estimates are used; display limits never affect volume. Volume errors do not block rendering.
         */
        VolumeOptions effectiveVolume=volume;
        if(structures.columns.enabled) effectiveVolume.columns=result.columns.items;
        if(volume.enabled && !structures.enabled) {
            result.volume.enabled=true;result.volume.error="Enable structure removal before volume calculation";
        } else result.volume=calculateStockpileVolumes(points,generated,effectiveFill.footprint,effectiveFill.footprintMargin,effectiveVolume);

        if(result.volume.groundValid && result.volume.error.isEmpty() && structures.enabled && structures.surfaceFilter.enabled) {
            result.surfaceFilter.volumeAfter=result.volume.totalVolume;
            if(surfaceRejectedMeasured.empty()&&surfaceRejectedEstimated.empty()) {
                result.surfaceFilter.volumeBefore=result.volume.totalVolume;result.surfaceFilter.comparisonValid=true;
            } else {
                // 比較同一地面與同一柱體模型，避免把地面重估差異誤當作移除點的影響。
                // Hold ground and column models fixed, isolating point-removal effects from ground refitting.
                const auto n=points.size(),m=generated.size();
                points.insert(points.end(),surfaceRejectedMeasured.begin(),surfaceRejectedMeasured.end());generated.insert(generated.end(),surfaceRejectedEstimated.begin(),surfaceRejectedEstimated.end());
                auto comparison=effectiveVolume;comparison.automaticGround=false;comparison.ground=result.volume.ground;
                const auto before=calculateStockpileVolumes(points,generated,effectiveFill.footprint,effectiveFill.footprintMargin,comparison);
                points.resize(n);generated.resize(m);
                if(before.error.isEmpty()) {result.surfaceFilter.volumeBefore=before.totalVolume;result.surfaceFilter.comparisonValid=true;}
            }
        }
        /** @brief 以同一份末段清理資料／地面／場域／柱體另算網格積分，不改寫既有柵格結果。
         *         Independently integrate a TIN over the same filtered data, ground, site and columns, preserving raster results.
         * 柱內虛擬料面仍保留供相交積分，之後才移除顯示／匯出點雲內的柱點。
         * Retain virtual surfaces inside columns for intersection integration before removing display/export points.
         */
        if(volume.enabled&&volume.mesh.enabled) {
            if(result.volume.groundValid&&result.volume.error.isEmpty()) result.meshVolume=calculateMeshVolume(points,generated,effectiveFill.footprint,effectiveFill.footprintMargin,result.volume.ground,effectiveVolume);
            else {result.meshVolume.enabled=true;result.meshVolume.error="Mesh comparison requires a valid shared ground plane";}
        }
        // 最終補點清理與繪圖共用正式地面；前置實測清理仍維持既定流程順序。
        // Final fill cleanup and rendering use the official ground, preserving early measurement-cleanup order.
        if(result.volume.groundValid) setColumnGround(result.columns.items,result.volume.ground);
        for(auto &column:result.columns.items) for(const auto &part:result.volume.columnVolumes) if(column.id==part.first) column.excludedVolume=part.second;
        // 虛擬料面先供柱料相交積分，再從視覺化與匯出移除柱體內補點，避免把柱子補回來。
        // Virtual material supports first inform displacement integration; remove estimates inside columns before display/export.
        if(structures.columns.enabled) result.columns.removedEstimated=removeColumnPoints(generated,result.columns.items);
        result.filledCount=generated.size();
        /** @brief 步驟 13d：在顯示抽樣前匯出完整量測及補點；顯示上限不減少輸出點數。
         *         Step 13d: export full measurements and estimates before display sampling; display caps never reduce exported counts.
         * 匯出錯誤獨立記錄，仍保留成功的處理與視覺化結果。
         * Export failures are reported separately, retaining successfully processed visualization data.
         */
        if(!exportPath.isEmpty()) {
            try {exportProcessedPly(exportPath,points,generated,calibration.applied);result.exportedPath=QFileInfo(exportPath).absoluteFilePath();}
            catch(const std::exception &e) {result.exportError=QString::fromUtf8(e.what());}
        }
        if (displayLimit && !generated.empty() && points.size()+generated.size()>displayLimit) {
            const size_t total = points.size()+generated.size();
            const size_t budget = std::min(displayLimit,total);
            // 少量補點優先完整顯示；最多預留一成預算，避免淹沒在大量量測點的抽樣中。
            // Show small fills in full, reserving up to 10% of the budget so they are not lost among measured samples.
            const size_t fillBudget = std::min(generated.size(),std::max(size_t(1),budget/10));
            result.filled = sample(generated,fillBudget);
            // sample 的 0 表示全部，零預算必須另行處理。 / sample(0) means all, so handle zero budgets explicitly.
            if (budget>fillBudget) result.after=sample(points,budget-fillBudget);
        } else {
            result.after=sample(points,displayLimit);
            result.filled=std::move(generated);
        }
    /** @brief 步驟 14：將標準例外轉為結果錯誤，供 GUI 或 CLI 顯示。 / Step 14: convert standard exceptions into a result error for GUI or CLI reporting.
     * 失敗時欄位可能只有部分資料；呼叫端必須先檢查 error，再使用結果。
     * Failure may leave partially populated fields; callers must check error before consuming the result.
     */
    } catch (const std::exception &e) { result.error = QString::fromUtf8(e.what()); }
    return result;
}
