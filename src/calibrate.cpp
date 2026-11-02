/** @file calibrate.cpp
 * @brief 由固定屋頂與牆面建立多次掃描共用的建築物座標系。
 *        Establish a shared building frame from static roof and wall geometry across scans.
 */
#include "config_paths.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QDateTime>
#include <QSaveFile>
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonDocument>
#include <QCryptographicHash>
#include <QLockFile>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

using V=Eigen::Vector3d;
using M=Eigen::Matrix4d;
/** @brief 本專案 PLY 的原始記錄；強度原位保留。 / Native project PLY record; intensity is preserved verbatim. */
struct Record {float x,y,z,intensity;};
/** @brief 平面 n·p+d=0 與擬合品質。 / Plane n·p+d=0 with fitting quality. */
struct Plane {V n; double d=0; size_t count=0; double rmse=0;};
/** @brief 來源串流標頭，限定實際三份掃描的無損記錄配置。 / Stream header restricted to the lossless layout of these three scans. */
struct Header {std::string text; size_t count=0;};

/** @brief 驗證 binary little-endian float XYZ/intensity，不支援的配置明確拒絕而不遺失欄位。
 *         Validate binary little-endian float XYZ/intensity, rejecting unsupported layouts instead of dropping fields.
 * @param in 原始二進位串流。 / Source binary stream.
 * @return 原始標頭與記錄數。 / Original header and record count.
 */
Header header(std::istream &in) {
    Header h; std::string line; std::vector<std::string> properties; bool vertex=false,ended=false;
    while (std::getline(in,line)) {
        h.text+=line+"\n";
        std::istringstream row(line); std::string tag; row>>tag;
        if (tag=="format" && line.find("binary_little_endian 1.0")==std::string::npos) throw std::runtime_error("Calibration export requires binary little-endian PLY");
        if (tag=="element") {std::string name; size_t n; row>>name>>n; vertex=name=="vertex"; if (vertex) h.count=n;}
        if (tag=="property" && vertex) properties.push_back(line);
        if (tag=="end_header") {ended=true; break;}
    }
    if (!ended || !h.count || properties!=std::vector<std::string>{"property float x","property float y","property float z","property float intensity"})
        throw std::runtime_error("Calibration requires float x/y/z/intensity vertex layout");
    return h;
}
/** @brief 串流抽樣，避免同時持有多份完整點雲。 / Stream uniform samples without retaining full scans.
 * @param path 原始 PLY。 / Source PLY.
 * @param signature 可選的完整來源 SHA-256。 / Optional SHA-256 of the complete source.
 * @return 有限座標樣本。 / Finite-coordinate samples.
 */
std::vector<V> sample(const QString &path,QByteArray *signature=nullptr) {
    std::ifstream in(path.toStdString(),std::ios::binary); if (!in) throw std::runtime_error("Cannot read source PLY");
    const auto h=header(in); const size_t step=std::max(size_t(1),h.count/80000);
    QCryptographicHash hash(QCryptographicHash::Sha256); hash.addData(h.text.data(),qint64(h.text.size()));
    std::vector<Record> buffer(65536); std::vector<V> points; points.reserve(h.count/step+1);
    for (size_t start=0;start<h.count;) {
        const size_t n=std::min(buffer.size(),h.count-start);
        if (!in.read(reinterpret_cast<char*>(buffer.data()),n*sizeof(Record))) throw std::runtime_error("Truncated PLY");
        hash.addData(reinterpret_cast<const char*>(buffer.data()),qint64(n*sizeof(Record)));
        for (size_t i=0;i<n;++i) if ((start+i)%step==0) {
            const auto &p=buffer[i]; if (std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)) points.emplace_back(p.x,p.y,p.z);
        }
        start+=n;
    }
    std::array<char,65536> tail{};
    while (in.read(tail.data(),tail.size()) || in.gcount()) hash.addData(tail.data(),in.gcount());
    if (signature) *signature=hash.result().toHex();
    return points;
}
/** @brief 對平面內點以 PCA 最小特徵向量重估法向。 / Refine a plane with the smallest PCA eigenvector.
 * @param points 平面內點。 / Plane inliers.
 * @return 朝上的單位法向與平面誤差。 / Upward unit normal and residual statistics.
 */
Plane fit(const std::vector<V> &points) {
    if (points.size()<30) throw std::runtime_error("Insufficient structural plane points");
    V mean=V::Zero(); for (const auto &p:points) mean+=p; mean/=points.size();
    Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero(); for (const auto &p:points) {const V q=p-mean; covariance+=q*q.transpose();}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
    Plane plane{solver.eigenvectors().col(0),0,points.size(),0};
    if (plane.n.z()<0) plane.n=-plane.n;
    plane.d=-plane.n.dot(mean);
    for (const auto &p:points) plane.rmse+=std::pow(plane.n.dot(p)+plane.d,2);
    plane.rmse=std::sqrt(plane.rmse/points.size()); return plane;
}
/** @brief RANSAC 擷取大型向上斜面，只使用固定屋頂而不配準變動料堆。
 *         Extract a dominant upward roof plane with RANSAC instead of registering variable piles.
 * @param points 輸入樣本。 / Input samples.
 * @param seed 可重現的亂數種子。 / Reproducible random seed.
 * @return 平面模型。 / Plane model.
 */
Plane roofPlane(const std::vector<V> &points,unsigned seed) {
    if (points.size()<300) throw std::runtime_error("Too few finite points to identify static roof geometry");
    std::vector<V> test; for (size_t i=0;i<points.size();i+=std::max(size_t(1),points.size()/18000)) test.push_back(points[i]);
    std::mt19937 rng(seed); std::uniform_int_distribution<size_t> select(0,test.size()-1);
    Plane best; size_t bestCount=0;
    for (int iteration=0;iteration<5000;++iteration) {
        const V a=test[select(rng)],b=test[select(rng)],c=test[select(rng)];
        V n=(b-a).cross(c-a); if (n.norm()<1e-8) continue; n.normalize(); if (n.z()<0)n=-n;
        if (n.z()<.85) continue;
        const double d=-n.dot(a); size_t count=0;
        for (const auto &p:test) if (std::abs(n.dot(p)+d)<.15) ++count;
        if (count>bestCount) {best={n,d,count,0}; bestCount=count;}
    }
    if (bestCount<200) throw std::runtime_error("No dominant roof plane found");
    for (int iteration=0;iteration<4;++iteration) {
        std::vector<V> inliers; for (const auto &p:points) if (std::abs(best.n.dot(p)+best.d)<.15) inliers.push_back(p);
        best=fit(inliers);
    }
    return best;
}
/** @brief 排序值的線性百分位數。 / Linearly interpolated quantile of a sorted value array.
 * @param sorted 已排序數值。 / Sorted values.
 * @param fraction 0 到 1。 / Fraction from zero to one.
 * @return 百分位值。 / Quantile value.
 */
double quantile(const std::vector<double> &sorted,double fraction) {
    const double pos=(sorted.size()-1)*fraction; const size_t i=size_t(pos);
    return sorted[i]+(sorted[std::min(i+1,sorted.size()-1)]-sorted[i])*(pos-i);
}
/** @brief 在投影外側區域找穩定牆面密度峰，再對平面內點擬合。
 *         Find a wall-density peak in an outer projection range, then fit the structural inliers.
 * @param points 原始樣本。 / Raw samples.
 * @param direction 預期牆面法向。 / Expected wall normal.
 * @param upper 搜尋較大座標一側。 / Search the higher-coordinate side.
 * @return 牆面模型。 / Wall plane.
 */
Plane wall(const std::vector<V> &points,const V &direction,bool upper) {
    std::vector<double> values; values.reserve(points.size()); for (const auto &p:points) values.push_back(direction.dot(p));
    std::sort(values.begin(),values.end());
    const double low=quantile(values,upper ? .75:.005),high=quantile(values,upper ? .995:.25);
    if (high<=low || high-low>2000) throw std::runtime_error("Unsupported wall search extent");
    const double width=.1; std::vector<size_t> counts(size_t(std::ceil((high-low)/width))+1,0);
    for (double v:values) if (v>=low&&v<=high) ++counts[size_t((v-low)/width)];
    const auto max=std::max_element(counts.begin(),counts.end());
    const double center=low+(std::distance(counts.begin(),max)+.5)*width;
    std::vector<V> inliers; for (const auto &p:points) if (std::abs(direction.dot(p)-center)<.13) inliers.push_back(p);
    auto result=fit(inliers);
    if (result.n.dot(direction)<0) {result.n=-result.n; result.d=-result.d;}
    if (result.n.dot(direction)<.995) throw std::runtime_error("Wall peak is not a stable plane perpendicular to the building axis");
    return result;
}
/** @brief 尺寸不符時，從多個外側牆峰尋找符合永久建築跨度的平面對。
 *         When dominant peaks give a wrong span, select an outer-plane pair matching the frozen building span.
 * @param points 擬合樣本，不含保留驗證點。 / Fitting samples, excluding held-out validation points.
 * @param direction 建築軸向單位向量。 / Unit building-axis direction.
 * @param expected 永久模型的牆距（公尺）。 / Frozen wall separation in meters.
 * @param tolerance 原有尺寸允差（公尺），不放寬。 / Existing dimension tolerance in meters, not relaxed.
 * @return 低側／高側牆面；無合格對時明確失敗。 / Lower/upper planes; fail explicitly if no pair qualifies.
 * @details 延伸料地端牆可能比主建築牆更密，不能只挑單一最大峰。
 *          Extension walls may dominate density; the highest peak alone need not bound the main building.
 */
std::array<Plane,2> fixedSpanWalls(const std::vector<V> &points,const V &direction,double expected,double tolerance) {
    std::vector<double> values;for(const auto &p:points) values.push_back(direction.dot(p));
    std::sort(values.begin(),values.end());
    auto candidates=[&](bool upper) {
        const double low=quantile(values,upper?.5:.005),high=quantile(values,upper?.995:.5),step=.1;
        if(high<=low||high-low>2000) throw std::runtime_error("Unsupported fixed-span wall search extent");
        std::vector<size_t> counts(size_t(std::ceil((high-low)/step))+1,0);
        for(double v:values) if(v>=low&&v<=high) ++counts[size_t((v-low)/step)];
        std::vector<Plane> result;
        for(int attempt=0;attempt<24;++attempt) {
            const auto peak=std::max_element(counts.begin(),counts.end());if(*peak<30) break;
            const size_t index=size_t(peak-counts.begin());const double center=low+(index+.5)*step;
            for(size_t j=0;j<counts.size();++j) if(std::abs(double(j)-double(index))*step<=.4) counts[j]=0;
            std::vector<V> inliers;for(const auto &p:points) if(std::abs(direction.dot(p)-center)<.13) inliers.push_back(p);
            if(inliers.size()<50) continue;
            auto plane=fit(inliers);if(plane.n.dot(direction)<0) {plane.n=-plane.n;plane.d=-plane.d;}
            if(plane.n.dot(direction)<.995||plane.rmse>.12) continue;
            result.push_back(plane);
        }
        return result;
    };
    const auto low=candidates(false),high=candidates(true);
    double best=-1;std::array<Plane,2> pair{};
    for(const auto &a:low) for(const auto &b:high) {
        const double span=-b.d/b.n.dot(direction)+a.d/a.n.dot(direction),error=std::abs(span-expected);
        if(error>tolerance) continue;
        const double score=double(std::min(a.count,b.count))/(1+error/tolerance);
        if(score>best) {best=score;pair={a,b};}
    }
    if(best<0) throw std::runtime_error("No supported wall pair matches frozen site dimensions; source remains unregistered");
    return pair;
}
/** @brief 每次掃描的結構模型。 / Structural model of one scan. */
struct Building {Plane roof[2],walls[4]; Eigen::Matrix3d axes; V origin; M transform; double width,length,pitch;};

/** @brief 保留樣本對固定結構的殘差，與擬合資料分開。 / Held-out residuals against fixed structure, separate from fitting data. */
struct Quality {size_t eligible=0,matched=0; double rmse=0,p95=0;};

/** @brief 由等斜率雙屋面建立右手座標系，並以三平面交點選定同一側屋頂角點。
 *         Build a right-handed frame from equal-pitch roof faces and locate one roof corner by three-plane intersection.
 * @param points 原始掃描樣本。 / Raw scan samples.
 * @param pitch 已知屋頂斜率。 / Known roof pitch.
 * @param reference 可選的永久模型，以既有跨度辨別延伸區端牆。 / Optional frozen model to distinguish extension walls using existing spans.
 * @return 建築座標系與靜態平面。 / Building frame and static planes.
 */
Building model(const std::vector<V> &points,double pitch,const Building *reference=nullptr) {
    Building b;
    b.roof[0]=roofPlane(points,42);
    std::vector<V> remaining; for (const auto &p:points) if (std::abs(b.roof[0].n.dot(p)+b.roof[0].d)>.25) remaining.push_back(p);
    b.roof[1]=roofPlane(remaining,71);
    b.pitch=std::tan(.5*std::acos(std::clamp(b.roof[0].n.dot(b.roof[1].n),-1.,1.)));
    if (std::abs(b.pitch-pitch)>.035) throw std::runtime_error("Detected roof pair disagrees with the known pitch");
    V z=(b.roof[0].n+b.roof[1].n).normalized(),x=(b.roof[1].n-b.roof[0].n).normalized();
    if (x.x()<0)x=-x;
    const V y=z.cross(x).normalized(); x=y.cross(z).normalized();
    b.axes.row(0)=x; b.axes.row(1)=y; b.axes.row(2)=z;
    b.walls[0]=wall(points,x,false); b.walls[1]=wall(points,x,true);
    b.walls[2]=wall(points,y,false); b.walls[3]=wall(points,y,true);
    /** @brief 使用固定尺寸消除主牆與延伸區端牆的峰值混淆；初始化場域仍保留原行為。
     *         Resolve main-wall versus extension-wall peaks using frozen dimensions; site initialization is unchanged.
     */
    if(reference) for(int axis=0;axis<2;++axis) {
        const V direction=axis?y:x;const int i=axis*2;
        const double span=-b.walls[i+1].d/b.walls[i+1].n.dot(direction)+b.walls[i].d/b.walls[i].n.dot(direction);
        const double expected=axis?reference->length:reference->width,tolerance=axis?1.:.7;
        if(std::abs(span-expected)>tolerance) {
            // 屋面切片也能形成假的窄投影平面，搜尋牆峰前先排除兩個已知屋面。
            // Roof slices can mimic projected planes; exclude both fitted roofs before searching wall peaks.
            std::vector<V> wallSamples;
            for(const auto &p:points) if(std::abs(b.roof[0].n.dot(p)+b.roof[0].d)>.25&&std::abs(b.roof[1].n.dot(p)+b.roof[1].d)>.25) wallSamples.push_back(p);
            const auto pair=fixedSpanWalls(wallSamples,direction,expected,tolerance);b.walls[i]=pair[0];b.walls[i+1]=pair[1];
            std::cerr<<"Re-selected "<<(axis?"end":"side")<<" walls against frozen span "<<expected<<" m (dominant peaks gave "<<span<<" m)\n";
        }
    }
    // 僅用固定屋頂與牆面校準；不將料堆形狀用於跨次掃描對齊。
    // Use fixed roofs and walls only; do not align changing stockpile geometry across scans.
    const double left=-b.walls[0].d/b.walls[0].n.dot(x),right=-b.walls[1].d/b.walls[1].n.dot(x);
    const double far=-b.walls[2].d/b.walls[2].n.dot(y),near=-b.walls[3].d/b.walls[3].n.dot(y);
    b.width=right-left; b.length=near-far;
    if (b.width<=0 || b.length<b.width*1.5) throw std::runtime_error("Detected long/short walls are inconsistent");
    // 選左側長牆所接的屋面：其法向 X 分量應為負。
    // Select the roof face adjoining the lower-X long wall: its normal has negative X component.
    if (b.roof[0].n.dot(x)>b.roof[1].n.dot(x)) std::swap(b.roof[0],b.roof[1]);
    Eigen::Matrix3d equation; equation.row(0)=b.walls[0].n; equation.row(1)=b.walls[2].n; equation.row(2)=b.roof[0].n;
    b.origin=equation.fullPivLu().solve(V(-b.walls[0].d,-b.walls[2].d,-b.roof[0].d));
    b.transform=M::Identity(); b.transform.block<3,3>(0,0)=b.axes; b.transform.block<3,1>(0,3)=-b.axes*b.origin;
    return b;
}

/** @brief 將平面轉入共同建築座標。 / Transform a plane into the shared building frame.
 * @param plane 原始平面。 / Source plane.
 * @param transform 來源到建築的剛體矩陣。 / Source-to-building rigid matrix.
 * @return 共同座標中的平面。 / Plane in shared coordinates.
 */
Plane transformedPlane(const Plane &plane,const M &transform) {
    const V n=transform.block<3,3>(0,0)*plane.n;
    return {n,plane.d-n.dot(transform.block<3,1>(0,3)),plane.count,plane.rmse};
}
/** @brief 只依固定結構平面做平衡權重 Huber 點對平面剛體精修。
 *         Refine rigid alignment with balanced Huber point-to-plane residuals on static structure only.
 * @param points 待配準樣本。 / Source samples.
 * @param reference 第一份掃描的固定結構。 / Static structure from the first scan.
 * @param initial 幾何初始對齊。 / Geometric initial alignment.
 * @return 精修後剛體矩陣。 / Refined rigid matrix.
 */
M refine(const std::vector<V> &points,const Building &reference,M initial) {
    std::vector<Plane> targets;
    for (const auto &p:reference.walls) targets.push_back(transformedPlane(p,reference.transform));
    for (const auto &p:reference.roof) targets.push_back(transformedPlane(p,reference.transform));
    for (int iteration=0;iteration<20;++iteration) {
        struct Match {V p; int target; double residual;}; std::vector<Match> matches;
        std::array<size_t,6> counts{};
        for (const auto &raw:points) {
            const V p=initial.block<3,3>(0,0)*raw+initial.block<3,1>(0,3);
            if (p.x()<-.5 || p.x()>reference.width+.5 || p.y()<-.5 || p.y()>reference.length+.5 || p.z()>15 || p.z()< -60) continue;
            int best=-1; double distance=.3;
            for (int i=0;i<6;++i) {
                if (i>=4 && p.z()< -1) continue;
                const double d=std::abs(targets[i].n.dot(p)+targets[i].d);
                if (d<distance) {best=i; distance=d;}
            }
            if (best>=0) {matches.push_back({p,best,targets[best].n.dot(p)+targets[best].d}); ++counts[best];}
        }
        for (auto count:counts) if (count<30) throw std::runtime_error("Insufficient static plane overlap for rigid refinement");
        Eigen::Matrix<double,6,6> h=Eigen::Matrix<double,6,6>::Zero(); Eigen::Matrix<double,6,1> g=Eigen::Matrix<double,6,1>::Zero();
        for (const auto &m:matches) {
            const auto &n=targets[m.target].n;
            Eigen::Matrix<double,6,1> j; j.head<3>()=m.p.cross(n); j.tail<3>()=n;
            const double weight=std::min(1.,.05/std::max(1e-12,std::abs(m.residual)))/counts[m.target];
            h+=weight*j*j.transpose(); g+=weight*j*m.residual;
        }
        const Eigen::Matrix<double,6,1> delta=-(h+Eigen::Matrix<double,6,6>::Identity()*1e-9).ldlt().solve(g);
        if (!delta.allFinite() || delta.head<3>().norm()>.05 || delta.tail<3>().norm()>2) throw std::runtime_error("Unstable structural alignment");
        M correction=M::Identity(); const double angle=delta.head<3>().norm();
        if (angle>1e-12) correction.block<3,3>(0,0)=Eigen::AngleAxisd(angle,delta.head<3>()/angle).toRotationMatrix();
        correction.block<3,1>(0,3)=delta.tail<3>(); initial=correction*initial;
        if (delta.norm()<1e-7) break;
    }
    return initial;
}

/** @brief 對未參與擬合的樣本檢查建築結構平面殘差，不把料堆當對應。
 *         Evaluate structural plane residuals on held-out samples, without treating stockpiles as correspondences.
 * @param points 保留的原始樣本。 / Held-out source samples.
 * @param reference 共同座標的參考建築。 / Reference building defining the shared frame.
 * @param transform 來源到共同座標。 / Source-to-shared transform.
 * @return 0.3 m 對應門檻內的殘差統計與樣本數；不是測量精度保證。
 *         Residuals/counts within the 0.3 m association gate, not a surveying accuracy guarantee.
 */
Quality quality(const std::vector<V> &points,const Building &reference,const M &transform) {
    Quality q; std::vector<Plane> targets;
    for (const auto &p:reference.walls) targets.push_back(transformedPlane(p,reference.transform));
    for (const auto &p:reference.roof) targets.push_back(transformedPlane(p,reference.transform));
    std::vector<double> residuals;
    for (const auto &raw:points) {
        const V p=transform.block<3,3>(0,0)*raw+transform.block<3,1>(0,3);
        if (p.x()<-.5 || p.x()>reference.width+.5 || p.y()<-.5 || p.y()>reference.length+.5 || p.z()>15 || p.z()< -60) continue;
        ++q.eligible; double distance=.3;
        for (int i=0;i<6;++i) {
            if (i>=4 && p.z()< -1) continue;
            distance=std::min(distance,std::abs(targets[i].n.dot(p)+targets[i].d));
        }
        if (distance<.3) {residuals.push_back(distance); q.rmse+=distance*distance;}
    }
    q.matched=residuals.size();
    if (q.matched) {q.rmse=std::sqrt(q.rmse/q.matched);std::sort(residuals.begin(),residuals.end());q.p95=quantile(residuals,.95);}
    return q;
}

/** @brief 串流輸出所有原始記錄，只轉換 XYZ；強度、點數與尾端資料不變。
 *         Stream every original record, transforming XYZ only while preserving intensity, counts, and trailing data.
 * @param source 原始檔。 / Source file.
 * @param target 新的校準檔，不覆寫既有資料。 / New calibrated file; existing files are not overwritten.
 * @param transform 來源到建築的矩陣。 / Source-to-building matrix.
 * @param identifier 永久場域識別碼，寫入輸出標頭。 / Persistent site ID stored in the output header.
 */
void exportCloud(const QString &source,const QString &target,const M &transform,const QString &identifier={}) {
    if (QFile::exists(target)) throw std::runtime_error("Calibrated output exists; refusing overwrite");
    std::ifstream in(source.toStdString(),std::ios::binary); const auto h=header(in);
    QSaveFile out(target); if (!out.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot create calibrated PLY");
    auto text=h.text;
    text.insert(text.find('\n')+1,"comment calibrated building frame; meters; origin at reference roof corner\n");
    if(!identifier.isEmpty()) text.insert(text.find('\n')+1,"comment calibration_site_id "+identifier.toStdString()+"\n");
    out.write(text.data(),qint64(text.size())); std::vector<Record> buffer(65536);
    for (size_t start=0;start<h.count;) {
        const size_t n=std::min(buffer.size(),h.count-start);
        if (!in.read(reinterpret_cast<char*>(buffer.data()),n*sizeof(Record))) throw std::runtime_error("Truncated PLY during export");
        for (size_t i=0;i<n;++i) {
            auto &p=buffer[i]; if (!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)) continue;
            const V q=transform.block<3,3>(0,0)*V(p.x,p.y,p.z)+transform.block<3,1>(0,3); p.x=float(q.x());p.y=float(q.y());p.z=float(q.z());
        }
        const qint64 bytes=qint64(n*sizeof(Record)); if (out.write(reinterpret_cast<const char*>(buffer.data()),bytes)!=bytes) throw std::runtime_error("Cannot write calibrated PLY");
        start+=n;
    }
    std::array<char,65536> tail{};
    while (in.read(tail.data(),tail.size()) || in.gcount()) if (out.write(tail.data(),in.gcount())!=in.gcount()) throw std::runtime_error("Cannot preserve PLY trailing data");
    if (!out.commit()) throw std::runtime_error("Cannot commit calibrated PLY");
}

/** @brief 讀取 JSON 物件，無效格式立即拒絕。 / Read a JSON object, rejecting invalid data.
 * @param path JSON 路徑。 / JSON path.
 * @return JSON 物件。 / JSON object.
 */
QJsonObject readObject(const QString &path) {
    QFile file(path); if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error(("Cannot read "+path).toStdString());
    QJsonParseError error; const auto doc=QJsonDocument::fromJson(file.readAll(),&error);
    if (error.error!=QJsonParseError::NoError || !doc.isObject()) throw std::runtime_error("Invalid JSON object");
    return doc.object();
}
/** @brief 原子提交 JSON，不留下半份矩陣。 / Atomically commit JSON without leaving partial matrices.
 * @param path 目標路徑。 / Destination.
 * @param object JSON 物件。 / JSON object.
 */
void saveObject(const QString &path,const QJsonObject &object) {
    QSaveFile out(path); if (!out.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot write JSON output");
    if (out.write(QJsonDocument(object).toJson())<0 || !out.commit()) throw std::runtime_error("Cannot commit JSON output");
}
/** @brief 將固定矩陣寫成列優先陣列。 / Serialize a matrix in row-major order.
 * @param matrix 剛體矩陣。 / Rigid matrix.
 * @return 16 個係數。 / Sixteen coefficients.
 */
QJsonArray matrixJson(const M &matrix) {
    QJsonArray values; for (int i=0;i<4;++i) for (int j=0;j<4;++j) values.append(matrix(i,j)); return values;
}
/** @brief 驗證可重複使用的剛體矩陣。 / Validate a reusable rigid transform.
 * @param values 列優先 JSON 陣列。 / Row-major JSON array.
 * @return 正交、右手、無縮放的矩陣。 / Orthonormal, right-handed matrix without scaling.
 */
M parseMatrix(const QJsonArray &values) {
    if (values.size()!=16) throw std::runtime_error("Invalid persistent transform");
    M m; for (int i=0;i<4;++i) for (int j=0;j<4;++j) {
        const auto v=values[4*i+j]; if (!v.isDouble()) throw std::runtime_error("Invalid transform coefficient"); m(i,j)=v.toDouble();
    }
    const Eigen::Matrix3d r=m.block<3,3>(0,0);
    if (!m.allFinite() || (m.row(3)-Eigen::RowVector4d(0,0,0,1)).norm()>1e-8 ||
        (r*r.transpose()-Eigen::Matrix3d::Identity()).norm()>1e-5 || std::abs(r.determinant()-1)>1e-5)
        throw std::runtime_error("Persistent matrix must be a proper rigid transform");
    return m;
}
/** @brief 對場域模型內容產生識別碼，防止誤用不同版本的座標系。 / Identify site geometry to prevent mixing frame versions.
 * @param object 不含 id 的場域內容。 / Site contents excluding id.
 * @return SHA-256 識別碼。 / SHA-256 identifier.
 */
QString siteId(QJsonObject object) {
    object.remove("id"); return QString::fromLatin1(QCryptographicHash::hash(QJsonDocument(object).toJson(QJsonDocument::Compact),QCryptographicHash::Sha256).toHex());
}
/** @brief 固定場域參考：平面已在共同座標中，不因新檔案重新取原點。 / Frozen site reference: planes already use the shared frame.
 * @param object 永久場域模型。 / Persistent site model.
 * @return 可供每次新掃描配準的參考。 / Reference for every future scan.
 */
Building loadSite(const QJsonObject &object) {
    if (object.value("kind").toString()!="stockyard_site" || object.value("version").toInt()!=1 || object.value("id").toString()!=siteId(object))
        throw std::runtime_error("Invalid or modified site model; initialize a new version deliberately");
    Building b; b.transform=M::Identity();b.axes=Eigen::Matrix3d::Identity();b.origin=V::Zero();
    b.width=object.value("width_m").toDouble();b.length=object.value("length_m").toDouble();b.pitch=object.value("roof_pitch").toDouble();
    if (b.width<=0 || b.length<b.width*1.5 || b.pitch<=0) throw std::runtime_error("Invalid site dimensions");
    const auto planes=object.value("planes").toArray(); if (planes.size()!=6) throw std::runtime_error("Site needs four walls and two roof planes");
    for (int i=0;i<6;++i) {
        const auto a=planes[i].toArray(); if (a.size()!=4) throw std::runtime_error("Invalid site plane");
        Plane p{V(a[0].toDouble(),a[1].toDouble(),a[2].toDouble()),a[3].toDouble(),0,0};
        if (!p.n.allFinite() || !std::isfinite(p.d) || std::abs(p.n.norm()-1)>1e-5) throw std::runtime_error("Invalid site normal");
        if (i<4) b.walls[i]=p; else b.roof[i-4]=p;
    }
    return b;
}
/** @brief 凍結一次參考模型，可承接既有三次校準的座標而不搬移原點。
 *         Freeze the reference once, preserving the frame from an existing calibration when available.
 * @param source 指定的參考掃描。 / Explicit reference scan.
 * @param output 不可覆寫的場域模型檔。 / Non-overwritable site model file.
 * @param pitch 已知等斜率。 / Known equal roof pitch.
 */
void initializeSite(const QString &source,const QString &output,double pitch) {
    if (QFile::exists(output)) throw std::runtime_error("Site model already exists; it is immutable");
    QByteArray signature; const auto points=sample(source,&signature); std::vector<V> fitting;
    for (size_t i=0;i<points.size();++i) if (i%3) fitting.push_back(points[i]);
    auto reference=model(fitting,pitch);
    const QFileInfo info(source);const auto manifestPath=info.dir().filePath("calibration.json");
    QJsonObject existing;
    if (QFile::exists(manifestPath)) {
        existing=readObject(manifestPath);
        if (existing.value("reference_file").toString()!=info.fileName()) throw std::runtime_error("Choose the existing reference_file to preserve its frame");
        const auto record=existing.value("files").toObject().value(info.fileName()).toObject();
        reference.transform=parseMatrix(record.value("source_to_building").toArray());
    }
    QJsonArray planes;
    for (const auto &p:reference.walls) {const auto q=transformedPlane(p,reference.transform);planes.append(QJsonArray{q.n.x(),q.n.y(),q.n.z(),q.d});}
    for (const auto &p:reference.roof) {const auto q=transformedPlane(p,reference.transform);planes.append(QJsonArray{q.n.x(),q.n.y(),q.n.z(),q.d});}
    QJsonObject site{{"version",1},{"kind","stockyard_site"},{"frame","Building: fixed roof corner / X across width / Y along length / Z up"},
        {"reference_file",info.fileName()},{"reference_sha256",QString::fromLatin1(signature)},{"reference_source_to_building",matrixJson(reference.transform)},
        {"width_m",reference.width},{"length_m",reference.length},{"roof_pitch",pitch},{"planes",planes},
        {"assumption","Equal-pitch roof bisector represents ground normal; approximately upright input, both roofs and end walls visible"},
        {"quality",QJsonObject{{"max_rmse_m",.18},{"max_p95_m",.28},{"min_match_fraction",.25},{"orientation_margin",.01}}}};
    site["id"]=siteId(site); loadSite(site);
    // 先驗證並雜湊既有註冊檔案，再綁定同一永久座標；不修改任一舊矩陣。
    // Validate/hash existing scans before binding their unchanged matrices to the persistent frame.
    if (!existing.isEmpty()) {
        auto records=existing.value("files").toObject();
        for (auto it=records.begin();it!=records.end();++it) {
            auto record=it.value().toObject(); parseMatrix(record.value("source_to_building").toArray());
            const auto path=info.dir().filePath(it.key());
            if (!QFile::exists(path) || qint64(record.value("source_bytes").toDouble())!=QFileInfo(path).size()) throw std::runtime_error("Existing scan changed or is missing");
            QByteArray hash; if (it.key()==info.fileName()) hash=signature;else sample(path,&hash);
            record["source_sha256"]=QString::fromLatin1(hash);record["site_id"]=site.value("id");
            const auto exported=info.dir().filePath("calibrated/"+it.key());
            if (QFile::exists(exported)) {record["exported_source_sha256"]=QString::fromLatin1(hash);record["exported_path"]=QFileInfo(exported).absoluteFilePath();}
            record["source_modified_ms"]=double(QFileInfo(path).lastModified().toMSecsSinceEpoch());it.value()=record;
        }
        existing["files"]=records;existing["site_id"]=site.value("id");existing["site_model"]=QFileInfo(output).absoluteFilePath();
    }
    saveObject(output,site);
    if (!existing.isEmpty()) saveObject(manifestPath,existing);
    std::cout<<"Frozen site "<<site.value("id").toString().toStdString()<<" at "<<output.toStdString()<<"\n";
}
/** @brief 新掃描的兩個可能長軸方向，利用固定屋頂不對稱與結構重疊選擇；不確定則拒絕。
 *         Evaluate both long-axis directions using fixed roof asymmetry/overlap, rejecting ambiguity.
 * @param points 擬合樣本。 / Fitting samples.
 * @param heldout 獨立檢核樣本。 / Held-out samples.
 * @param reference 固定場域模型。 / Frozen site reference.
 * @param source 此次掃描的初始模型。 / Model of this scan.
 * @param site 固定品質門檻。 / Frozen quality limits.
 * @param orientation auto/forward/reverse，對稱情況可由使用者提供明確方向。
 *                    auto/forward/reverse; explicit hints resolve symmetric cases.
 * @return 通過品質檢核的矩陣。 / Transform passing quality checks.
 */
M registerScan(const std::vector<V> &points,const std::vector<V> &heldout,const Building &reference,const Building &source,
               const QJsonObject &site,const QString &orientation) {
    if (std::abs(source.width-reference.width)>.7 || std::abs(source.length-reference.length)>1.)
        throw std::runtime_error("Site dimensions mismatch or partial architecture; not enough evidence to register");
    struct Candidate {M t;double score;}; std::vector<Candidate> candidates;
    for (int direction=0;direction<2;++direction) {
        if ((orientation=="forward" && direction) || (orientation=="reverse" && !direction)) continue;
        M initial=source.transform;
        if (direction) {
            const auto right=transformedPlane(source.roof[1],source.transform);
            const double height=-(right.n.x()*source.width+right.n.y()*source.length+right.d)/right.n.z();
            M flip=M::Identity();flip(0,0)=-1;flip(1,1)=-1;flip(0,3)=source.width;flip(1,3)=source.length;flip(2,3)=-height;
            initial=flip*initial;
        }
        try {
            const M aligned=refine(points,reference,initial); const auto q=quality(heldout,reference,aligned);
            const auto limits=site.value("quality").toObject();
            const double fraction=q.eligible ? double(q.matched)/q.eligible:0;
            if (q.matched<200 || fraction<limits.value("min_match_fraction").toDouble() || q.rmse>limits.value("max_rmse_m").toDouble() || q.p95>limits.value("max_p95_m").toDouble()) continue;
            candidates.push_back({aligned,q.rmse+.3*(1-fraction)});
        } catch (const std::exception &) { /* 此方向失敗時仍檢查另一方向。 / Still evaluate the other direction. */ }
    }
    if (candidates.empty()) throw std::runtime_error("Registration failed static-plane quality checks; source remains unregistered");
    std::sort(candidates.begin(),candidates.end(),[](const Candidate &a,const Candidate &b){return a.score<b.score;});
    if (candidates.size()>1 && candidates[1].score-candidates[0].score<site.value("quality").toObject().value("orientation_margin").toDouble())
        throw std::runtime_error("180-degree orientation is ambiguous; supply --orientation forward/reverse or a surveyed landmark");
    return candidates.front().t;
}
/** @brief 永久場域初始化與未來掃描註冊入口。 / Persistent site initialization and future-scan registration entry point.
 * @param argc 參數數量。 / Argument count.
 * @param argv 檔案／目錄與固定場域選項。 / File/directory and fixed-site options.
 * @return 0 成功，1 失敗。 / Zero on success, one on failure.
 */
int main(int argc,char **argv) {
    QCoreApplication app(argc,argv); QCommandLineParser parser;parser.addHelpOption();
    parser.addPositionalArgument("input","Source PLY file or directory","[input]");
    parser.addOption({"site","Immutable site model","json"});parser.addOption({"init-site","Initialize an immutable site model once","json"});
    parser.addOption({"reference","Explicit source PLY for site initialization","ply"});
    parser.addOption({"roof-pitch","Known equal roof pitch for initialization","slope","0.2"});
    parser.addOption({"matrix-only","Save calibration matrices without exporting full PLY copies"});
    parser.addOption({"output-dir","Directory for new calibrated PLY copies","directory"});
    parser.addOption({"orientation","Resolve geometric direction: auto, forward, reverse","direction","auto"});
    parser.addOption({"analyze-only","Check registration without writing any results"});parser.process(app);
    try {
        bool ok;const double pitch=parser.value("roof-pitch").toDouble(&ok);
        if (!ok || pitch<=0 || pitch>.8) throw std::runtime_error("Invalid roof pitch");
        if (parser.isSet("init-site")) {
            if(parser.isSet("analyze-only"))throw std::runtime_error("--init-site writes a model and cannot be combined with --analyze-only");
            if (!parser.isSet("reference")) throw std::runtime_error("Initialization requires --reference source.ply");
            initializeSite(QFileInfo(parser.value("reference")).absoluteFilePath(),QFileInfo(parser.value("init-site")).absoluteFilePath(),pitch);return 0;
        }
        const auto args=parser.positionalArguments();if(args.size()>1) throw std::runtime_error("Provide one file or directory");
        const QFileInfo input(args.isEmpty()?QDir::currentPath():args[0]);
        if (!input.exists()) throw std::runtime_error("Input does not exist");
        QDir dir=input.isDir()?QDir(input.absoluteFilePath()):input.dir();QStringList files;
        if (input.isFile()) files<<input.fileName();else for(const auto &file:dir.entryList(QDir::Files,QDir::Name)) if(file.endsWith(".ply",Qt::CaseInsensitive)) files<<file;
        if(files.isEmpty()) throw std::runtime_error("No source PLY files");
        QString sitePath=parser.value("site");
        if(sitePath.isEmpty()) sitePath=resolveSiteModel(dir);
        const auto site=readObject(sitePath);const auto reference=loadSite(site);const auto identifier=site.value("id").toString();
        const QString orientation=parser.value("orientation");
        if(orientation!="auto"&&orientation!="forward"&&orientation!="reverse") throw std::runtime_error("Invalid orientation hint");
        const bool dry=parser.isSet("analyze-only");const auto manifestPath=dir.filePath("calibration.json");
        QLockFile lock(manifestPath+".lock");if(!dry && !lock.tryLock(5000)) throw std::runtime_error("Another registration is writing this directory");
        QJsonObject manifest=QFile::exists(manifestPath)?readObject(manifestPath):QJsonObject{};
        if(!manifest.isEmpty() && manifest.value("site_id").toString()!=identifier) throw std::runtime_error("Directory belongs to another/unversioned site frame; do not mix coordinate systems");
        auto records=manifest.value("files").toObject();
        const QDir output(parser.isSet("output-dir")?parser.value("output-dir"):dir.filePath("calibrated"));
        for(const auto &name:files) {
            const auto path=dir.filePath(name);QByteArray hash;const auto points=sample(path,&hash);
            const auto cached=records.value(name).toObject();M transform;
            const bool reuse=cached.value("source_sha256").toString()==QString::fromLatin1(hash)&&cached.value("site_id").toString()==identifier;
            if(reuse) {transform=parseMatrix(cached.value("source_to_building").toArray());std::cout<<"Cached fixed-frame registration: "<<name.toStdString()<<"\n";}
            else {
                std::vector<V> fitting,heldout;for(size_t i=0;i<points.size();++i)(i%3?fitting:heldout).push_back(points[i]);
                if(QString::fromLatin1(hash)==site.value("reference_sha256").toString()) transform=parseMatrix(site.value("reference_source_to_building").toArray());
                else {const auto source=model(fitting,site.value("roof_pitch").toDouble(),&reference);transform=registerScan(fitting,heldout,reference,source,site,orientation);}
                const auto q=quality(heldout,reference,transform);
                records[name]=QJsonObject{{"source_to_building",matrixJson(transform)},{"source_bytes",double(QFileInfo(path).size())},
                    {"source_modified_ms",double(QFileInfo(path).lastModified().toMSecsSinceEpoch())},
                    {"source_sha256",QString::fromLatin1(hash)},{"site_id",identifier},{"heldout_roi_samples",double(q.eligible)},
                    {"heldout_structure_matches",double(q.matched)},{"heldout_structure_rmse_m",q.rmse},{"heldout_structure_p95_m",q.p95}};
                std::cout<<"Registered "<<name.toStdString()<<" to frozen site; structure RMSE "<<q.rmse<<" m; P95 "<<q.p95<<" m\n";
            }
            if(!dry && !parser.isSet("matrix-only")) {
                const auto target=output.filePath(name);
                if(QFile::exists(target)) {
                    if(!reuse || cached.value("exported_source_sha256").toString()!=QString::fromLatin1(hash) || cached.value("exported_path").toString()!=QFileInfo(target).absoluteFilePath())
                        throw std::runtime_error("Output exists for changed/untracked scan; choose a new output directory");
                } else {
                    if(!QDir().mkpath(output.absolutePath())) throw std::runtime_error("Cannot create output directory");
                    exportCloud(path,target,transform,identifier);
                    auto record=records.value(name).toObject();record["exported_source_sha256"]=QString::fromLatin1(hash);record["exported_path"]=QFileInfo(target).absoluteFilePath();records[name]=record;
                }
            }
            if(!dry) {
                manifest["version"]=1;manifest["frame"]=site.value("frame");manifest["site_id"]=identifier;
                manifest["site_model"]=QFileInfo(sitePath).absoluteFilePath();manifest["reference_file"]=site.value("reference_file");manifest["files"]=records;saveObject(manifestPath,manifest);
            }
        }
    }catch(const std::exception &e){std::cerr<<e.what()<<"\n";return 1;}
    return 0;
}
