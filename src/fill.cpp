/**
 * @file fill.cpp
 * @brief 對 Z 朝上的單層堆體補償小型封閉遮擋缺口。 / Compensate small enclosed occlusion holes in Z-up single-surface piles.
 */
#include "cloud.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

std::vector<Point> fillStockpileHoles(const std::vector<Point> &points, const FillOptions &o, size_t &holes) {
    holes = 0;
    if (!o.enabled) return {};
    /** @brief 步驟 1：檢查 ROI 與資源限制，不在建築物全域自動推測堆體。
     *         Step 1: validate the ROI and resource limits without guessing pile locations throughout the building.
     */
    for (double v : {o.xMin,o.xMax,o.yMin,o.yMax,o.zMin,o.zMax,o.spacing,o.maxHole,o.maxError,o.maxSlope})
        if (!std::isfinite(v)) throw std::runtime_error("Non-finite hole-filling parameter");
    if (o.xMin>=o.xMax || o.yMin>=o.yMax || o.zMin>=o.zMax || o.spacing<.1 ||
        o.maxHole<o.spacing || o.maxError<=0 || o.maxSlope<=0)
        throw std::runtime_error("Select stockpile XY/Z bounds; spacing >= 0.1 m, max hole >= spacing, error/slope > 0");
    const double width = std::ceil((o.xMax-o.xMin)/o.spacing), height = std::ceil((o.yMax-o.yMin)/o.spacing);
    if (width<3 || height<3 || width*height>2000000)
        throw std::runtime_error("Fill ROI needs 3x3 .. 2,000,000 cells; reduce ROI or increase spacing");
    const int nx = int(width), ny = int(height);
    struct Cell {
        double sum = 0, sumX = 0, sumY = 0;
        float low = std::numeric_limits<float>::infinity(), high = -std::numeric_limits<float>::infinity();
        unsigned count = 0;
        bool visited = false;
    };
    std::vector<Cell> grid(size_t(nx)*ny);
    /** @brief 步驟 2：只有 ROI 內的量測點形成高度網格；單格高度跨度太大代表可能有多層表面。
     *         Step 2: build the height grid from measured ROI points only; excessive within-cell spread indicates possible multiple surfaces.
     * 單格取平均高度；高度跨度大於 2 倍允許誤差的格子視為不可信，不能當作補洞邊界。
     * Use mean cell height; cells spanning more than twice the allowed error cannot support filling.
     */
    for (const auto &p : points) {
        if (p.x<o.xMin || p.x>=o.xMax || p.y<o.yMin || p.y>=o.yMax || p.z<o.zMin || p.z>o.zMax) continue;
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) continue;
        int x = int(std::floor((double(p.x)-o.xMin)/o.spacing));
        int y = int(std::floor((double(p.y)-o.yMin)/o.spacing));
        if (x<0 || x>=nx || y<0 || y>=ny) continue;
        auto &c = grid[size_t(y)*nx+x]; c.sum+=p.z; c.sumX+=p.x; c.sumY+=p.y; ++c.count;
        c.low=std::min(c.low,p.z); c.high=std::max(c.high,p.z);
    }
    std::vector<Point> generated;
    std::vector<int> pending, region, boundary;
    /** @brief 步驟 3：以四鄰接找空格連通區；接觸 ROI 邊緣的區域不補，避免外插。
     *         Step 3: find four-connected empty regions; reject regions touching the ROI edge to avoid extrapolation.
     * 空格連通區不因新增補點而改變，補點不會成為其他缺口的量測依據。
     * Generated points never become observations for other holes.
     */
    for (int start=0; start<nx*ny; ++start) {
        if (grid[start].count || grid[start].visited) continue;
        pending.clear(); region.clear(); boundary.clear();
        pending.push_back(start); grid[start].visited=true;
        bool exterior=false, ambiguous=false, tooLarge=false;
        int minX=nx,maxX=0,minY=ny,maxY=0;
        while (!pending.empty()) {
            int id=pending.back(); pending.pop_back(); const int x=id%nx,y=id/nx;
            minX=std::min(minX,x); maxX=std::max(maxX,x); minY=std::min(minY,y); maxY=std::max(maxY,y);
            if (!x || !y || x==nx-1 || y==ny-1) exterior=true;
            // 持續標記整個空區，但限制待補點清單大小。 / Visit the entire region but bound the candidate list.
            if (region.size()<4096) region.push_back(id); else tooLarge=true;
            for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
                if (!dx && !dy) continue;
                int xx=x+dx,yy=y+dy;
                if (xx<0 || yy<0 || xx>=nx || yy>=ny) continue;
                int neighbor=yy*nx+xx; auto &c=grid[neighbor];
                if (c.count) {
                    if (double(c.high)-c.low>2*o.maxError) ambiguous=true;
                    if (boundary.size()<32768) boundary.push_back(neighbor); else tooLarge=true;
                } else if ((!dx || !dy) && !c.visited) { c.visited=true; pending.push_back(neighbor); }
            }
        }
        const double diameter=std::hypot(maxX-minX+1,maxY-minY+1)*o.spacing;
        if (exterior || ambiguous || tooLarge || diameter>o.maxHole) continue;
        std::sort(boundary.begin(),boundary.end()); boundary.erase(std::unique(boundary.begin(),boundary.end()),boundary.end());
        if (boundary.size()<4) continue;
        /** @brief 步驟 4：用缺口周圍量測格擬合 z=a*u+b*v+c，u/v 以缺口中心為基準。
         *         Step 4: fit z=a*u+b*v+c to surrounding measured cells, centering u/v on the hole.
         * 高斯消去法解 3x3 最小平方法；共線或病態邊界不補。
         * Solve the 3x3 least-squares system with Gaussian elimination; reject collinear/ill-conditioned boundaries.
         */
        const double cx=(minX+maxX)*.5, cy=(minY+maxY)*.5;
        const double centerX=o.xMin+(cx+.5)*o.spacing,centerY=o.yMin+(cy+.5)*o.spacing;
        double matrix[3][4]{};
        for (int id : boundary) {
            // 邊界使用量測點的平均 XY，不假定它們恰好位於格心。 / Use actual mean measured XY, not assumed cell centers.
            const double row[]{grid[id].sumX/grid[id].count-centerX,grid[id].sumY/grid[id].count-centerY,1};
            const double z=grid[id].sum/grid[id].count;
            for (int i=0;i<3;++i) { for (int j=0;j<3;++j) matrix[i][j]+=row[i]*row[j]; matrix[i][3]+=row[i]*z; }
        }
        bool singular=false;
        for (int i=0;i<3;++i) {
            int pivot=i;
            for (int j=i+1;j<3;++j) if (std::abs(matrix[j][i])>std::abs(matrix[pivot][i])) pivot=j;
            if (std::abs(matrix[pivot][i])<1e-10) { singular=true; break; }
            for (int j=0;j<4;++j) std::swap(matrix[i][j],matrix[pivot][j]);
            const double divisor=matrix[i][i]; for (int j=0;j<4;++j) matrix[i][j]/=divisor;
            for (int r=0;r<3;++r) if (r!=i) { const double factor=matrix[r][i]; for (int j=0;j<4;++j) matrix[r][j]-=factor*matrix[i][j]; }
        }
        if (singular) continue;
        const double a=matrix[0][3],b=matrix[1][3],c=matrix[2][3];
        auto predict=[&](int id) { return a*(id%nx-cx)*o.spacing+b*(id/nx-cy)*o.spacing+c; };
        /** @brief 步驟 5：限制坡度與邊界誤差，避免跨牆、跨堆或跨高度斷層補平面。
         *         Step 5: bound slope and boundary residuals to avoid bridging walls, separate piles, or height discontinuities.
         */
        if (std::hypot(a,b)>o.maxSlope) continue;
        bool acceptable=true;
        for (int id : boundary) {
            const auto &cell=grid[id];
            const double fitted=a*(cell.sumX/cell.count-centerX)+b*(cell.sumY/cell.count-centerY)+c;
            if (std::abs(fitted-cell.sum/cell.count)>o.maxError) acceptable=false;
        }
        for (int id : region) if (!std::isfinite(predict(id)) || predict(id)<o.zMin || predict(id)>o.zMax) acceptable=false;
        if (!acceptable) continue;
        /** @brief 步驟 6：只在空格中心新增估計點，不修改量測資料，也不捏造強度。
         *         Step 6: add estimates only at empty cell centers, without modifying measurements or inventing intensity.
         */
        for (int id : region) {
            const double x=o.xMin+(id%nx+.5)*o.spacing, y=o.yMin+(id/nx+.5)*o.spacing;
            if (x<o.xMax && y<o.yMax) generated.push_back({float(x),float(y),float(predict(id))});
        }
        ++holes;
    }
    return generated;
}
