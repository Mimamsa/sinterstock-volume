/**
 * @file main.cpp
 * @brief Qt5 雙視區點雲檢視器與命令列入口。 / Qt5 dual-pane point-cloud viewer and command-line entry point.
 */
#include "cloud.h"
#include <QApplication>
#include <QComboBox>
#include <QCheckBox>
#include <QCommandLineParser>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFutureWatcher>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QLinearGradient>
#include <QMatrix4x4>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPainter>
#include <QPushButton>
#include <QSurfaceFormat>
#include <QSignalBlocker>
#include <QTimer>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <QtConcurrent>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <QTableWidget>
#include <QHeaderView>
#include <memory>
#include <functional>

/**
 * @brief 以共享相機繪製處理前後點雲。 / Render before/after clouds with a shared camera.
 * @details 各側使用獨立 XYZ VBO；只有 dirty 時上傳資料，操作或視窗事件才要求重繪。
 *          Each pane has an XYZ VBO; data uploads occur only when dirty, and interaction/window events request repainting.
 * @note 所有 OpenGL 操作在 GUI 執行緒及有效 context 中進行。
 *       All OpenGL operations run on the GUI thread with a valid context.
 */
class Viewer : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    /// 點座標轉換與依強度／高度著色的著色器。 / Position-transform and intensity/height color shaders.
    QOpenGLShaderProgram shader;
    /// 頂點格式狀態。 / Vertex format state.
    QOpenGLVertexArrayObject vao;
    /// [0] 原始樣本，[1] 清理後量測樣本，[2] 合成補點樣本。 / [0] input, [1] cleaned measurements, [2] synthetic samples.
    QOpenGLBuffer buffers[4]; ///< 第四組為審核證據，不是處理結果。 / Fourth stream is review evidence, never material output.
    /// 保留 CPU 樣本，供資料上傳或 context 重建使用。 / Retain CPU samples for uploads and context recreation.
    std::shared_ptr<Cloud> cloud;
    /// dirty：待上傳資料；ready：GL 資源可用。 / dirty: upload pending; ready: GL resources usable.
    bool dirty = false, ready = false;
    /// 初始化失敗時顯示的診斷訊息。 / Diagnostic text displayed on initialization failure.
    QString glError;
    /// 共享相機的取景中心。 / Shared camera framing center.
    QVector3D center;
    /// 取景尺度、縮放比例與旋轉角度（度）。 / Framing scale, zoom factor, and rotation angles in degrees.
    float radius = 1, zoom = 1, yaw = -35, pitch = 25;
    /// 視平面平移與座標軸公尺長度。 / View-plane pan and coordinate-axis length in meters.
    float panX = 0, panY = 0, axisSize = 1;
    /// 上次滑鼠位置，用於計算拖曳位移。 / Last mouse position for drag deltas.
    QPoint last;
    /// 兩側實際上傳的點數。 / Actual uploaded point counts for both panes.
    int counts[4]{0, 0, 0, 0};
    int paneModes[2]{0,1};
    bool paneColumns[2]{false,false},paneGround[2]{true,true},showAxes=true;
    static bool visibleLayer(int mode,int layer) {return (mode==0&&layer==0)||(mode==1&&(layer==1||layer==2))||(mode==2&&layer==1)||(mode==3&&layer==2)||(mode==4&&layer==3);}

    bool highlightFill = true, selectingROI = false, draggingROI = false;
    float roofHideDepth=3.f; ///< 屋頂向下的垂直顯示裁切深度。 / Vertical display clipping depth below the roof.
    bool hideLeftRoof=false; ///< 左側 GPU 屋頂裁切，不改變點雲。 / Left-pane GPU clipping, never mutating clouds.
    bool zSlice = false;
    float sliceMin = 0, sliceMax = 0;
    int roiPane = 0;
    QPoint roiStart, roiEnd;
    std::function<void(double,double,double,double)> roiCallback;
    /// 0 單色、1 強度、2 高度；兩側使用同一色階。 / 0 uniform, 1 intensity, 2 height; both panes share the same scale.
    int colorMode = 2;
    float colorMin = 0, colorMax = 1;

    /** @brief 將 context 設為 current 後釋放 GL 資源。 / Make the context current and release GL resources. */
    void cleanup() {
        if (!context()) return;
        makeCurrent();
        for (auto &buffer : buffers) buffer.destroy();
        vao.destroy(); shader.removeAllShaders(); ready = false;
        doneCurrent();
    }
    /**
     * @brief 建立兩側共用姿態的正交投影矩陣。 / Build an orthographic matrix with shared camera pose.
     * @param paneWidth 單側寬度，單位為 Qt 邏輯像素。 / Pane width in Qt logical pixels.
     * @return projection * view，將世界座標轉成裁切座標。 / projection * view, mapping world to clip coordinates.
     * @note 合成矩陣實際先將取景中心移到原點，再旋轉，最後在視平面平移。
     *       The composed view applies centering first, then rotation, then view-plane translation.
     */
    QMatrix4x4 camera(int paneWidth) const {
        QMatrix4x4 projection, view;
        float extent = radius * zoom;
        float aspect = float(paneWidth) / std::max(1, height());
        projection.ortho(-extent * aspect, extent * aspect, -extent, extent, -radius * 10, radius * 10);
        view.translate(panX, panY, 0);
        view.rotate(pitch, 1, 0, 0); view.rotate(yaw, 0, 0, 1);
        view.translate(-center);
        return projection * view;
    }
    /**
     * @brief 由 Qt 在 context 可用時初始化 GL 3.3 資源。 / Initialize GL 3.3 resources when Qt supplies a valid context.
     * @details 編譯著色器、建立 VAO/VBO，並註冊 context 銷毀時的清理回呼。
     *          Compile shaders, create VAO/VBOs, and register cleanup for context destruction.
     */
    void initializeGL() override {
        if (!initializeOpenGLFunctions()) { glError = "OpenGL 3.3 is required"; return; }
        connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, [this] { cleanup(); }, Qt::DirectConnection);
        /** @brief 以 uniform 切換單色／強度／Z；正規化後在藍紅間插值，無效值使用灰色。
         *         Switch uniform/intensity/Z modes via uniforms; normalize and interpolate blue/red, using gray for invalid values.
         * 先縮放上下限再相減，避免極端 float32 強度範圍相減溢位；等值範圍取中間色。
         * Scale limits before subtraction to avoid overflow for extreme float32 ranges; equal limits use the midpoint color.
         */
        bool ok = shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
            "#version 330 core\nlayout(location=0) in vec3 position; layout(location=1) in float intensity;"
            "uniform mat4 mvp; uniform int colorMode; uniform float colorMin; uniform float colorMax; uniform bool fillHighlight; uniform bool semanticHighlight;"
            "uniform bool hideRoof; uniform vec4 roof0; uniform vec4 roof1; uniform float roofDepth; uniform bool zSlice; uniform float sliceMin; uniform float sliceMax; uniform float pointSize; out vec3 pointColor;"
            "void main(){gl_Position=mvp*vec4(position,1);gl_PointSize=pointSize;"
            "if(zSlice&&(position.z<sliceMin||position.z>sliceMax))gl_Position=vec4(2,2,2,1);"
            "if(hideRoof&&(dot(roof0,vec4(position,1))>=-roofDepth*roof0.z||dot(roof1,vec4(position,1))>=-roofDepth*roof1.z))gl_Position=vec4(2,2,2,1);"
            "float value=colorMode==1?intensity:position.z;"
            "float scale=max(abs(colorMin),max(abs(colorMax),1.0));"
            "float span=colorMax/scale-colorMin/scale;"
            "float t=span>0.0?clamp((value/scale-colorMin/scale)/span,0.0,1.0):0.5;"
            "pointColor=semanticHighlight?vec3(1.0,.1,.8):fillHighlight?vec3(1.0,.65,.05):(colorMode==0?vec3(.6,.78,.95):"
            "(isnan(value)||isinf(value)?vec3(.55):mix(vec3(0.1,0.25,1.0),vec3(1.0,0.15,0.05),t)));}") &&
            shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
            "#version 330 core\nin vec3 pointColor; out vec4 frag; void main(){frag=vec4(pointColor,1);}") && shader.link();
        if (!ok) { glError = shader.log(); return; }
        vao.create(); for (auto &buffer : buffers) buffer.create();
        glEnable(GL_DEPTH_TEST); glEnable(GL_PROGRAM_POINT_SIZE);
        ready = true; dirty = true;
    }
    /**
     * @brief 上傳必要資料、繪製雙視區，再疊加座標與文字。 / Upload pending data, draw both panes, then overlay axes and text.
     * @details 資料上傳只在 dirty 時發生；平移旋轉只更新相機矩陣。
     *          Uploads occur only when dirty; pan/rotation only update camera matrices.
     */
    void paintGL() override {
        if (!ready) {
            QPainter p(this); p.fillRect(rect(), Qt::black); p.setPen(Qt::white);
            p.drawText(rect(), Qt::AlignCenter, glError); return;
        }
        /** @brief 恢復可能被 QPainter 改動的 GL 狀態，並清除色彩與深度。
         *         Restore GL state potentially changed by QPainter and clear color/depth buffers.
         */
        glEnable(GL_DEPTH_TEST); glEnable(GL_PROGRAM_POINT_SIZE); glDisable(GL_SCISSOR_TEST);
        glClearColor(.045f, .055f, .07f, 1); glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        vao.bind(); shader.bind();
        /** @brief 各側樣本只在資料變動或 context 重建後重新上傳。
         *         Re-upload each pane's sample only after data changes or context recreation.
         */
        if (cloud && dirty) {
            const std::vector<Point> *samples[]{&cloud->before, &cloud->after, &cloud->filled, &cloud->columnEvidence};
            for (int i = 0; i < 4; ++i) {
                buffers[i].bind(); buffers[i].setUsagePattern(QOpenGLBuffer::StaticDraw);
                buffers[i].allocate(samples[i]->data(), int(samples[i]->size() * sizeof(Point)));
                counts[i] = int(samples[i]->size());
            }
            dirty = false;
        }
        /** @brief Qt 使用邏輯像素，GL viewport 使用實體像素；以 DPR 轉換後繪製兩側。
         *         Qt uses logical pixels and GL viewports use physical pixels; convert by DPR before drawing both panes.
         */
        const float dpr = devicePixelRatioF();
        // 兩側共用既有 VBO，切換顯示只改繪製層，不重算、不複製 GPU 點資料。
        // Both panes reuse existing VBOs; display changes never recompute or duplicate GPU point data.
        for(int pane=0;pane<2;++pane) for (int i = 0; i < 4; ++i) {
            if(!visibleLayer(paneModes[pane],i)) continue;
            int left = pane ? width() / 2 : 0;
            int paneWidth = pane ? width() - left : width() / 2;
            glViewport(int(left * dpr), 0, int(paneWidth * dpr), int(height() * dpr));
            const auto mvp = camera(paneWidth);
            shader.setUniformValue("mvp", mvp);
            // 僅左側裁切屋頂與其上方；使用雙斜面而非固定 Z 高度。
            // Clip only the left pane using both sloped roof planes rather than a fixed Z cutoff.
            shader.setUniformValue("hideRoof",pane==0&&hideLeftRoof&&cloud&&cloud->hasRoofPlanes);
            shader.setUniformValue("roofDepth",roofHideDepth);
            for(int roof=0;roof<2;++roof) {const auto a=cloud?cloud->roofPlanes[roof]:std::array<double,4>{};
                shader.setUniformValue(roof==0?"roof0":"roof1",QVector4D(float(a[0]),float(a[1]),float(a[2]),float(a[3])));}

            shader.setUniformValue("colorMode", colorMode);
            shader.setUniformValue("colorMin", colorMin);
            shader.setUniformValue("colorMax", colorMax);
            shader.setUniformValue("fillHighlight", i==2 && highlightFill);
            shader.setUniformValue("semanticHighlight",i==3);
            shader.setUniformValue("pointSize",i==2 && highlightFill ? 3.f : 1.5f);
            shader.setUniformValue("zSlice",zSlice);
            shader.setUniformValue("sliceMin",sliceMin); shader.setUniformValue("sliceMax",sliceMax);
            buffers[i].bind();
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point), nullptr);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(Point), reinterpret_cast<void *>(offsetof(Point, intensity)));
            glDrawArrays(GL_POINTS, 0, counts[i]);
            buffers[i].release();
            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
        }
        shader.release(); vao.release();
        /** @brief 使用相同矩陣將原點與軸端投影到畫面；2D 疊圖不受點雲深度遮擋。
         *         Project the origin and axis endpoints with the same matrix; 2D overlays remain visible over cloud depth.
         */
        QPainter p(this); p.setRenderHint(QPainter::Antialiasing);
        p.setPen(QColor(90, 90, 90)); p.drawLine(width()/2, 0, width()/2, height());
        for (int i = 0; i < 2; ++i) {
            int left = i ? width()/2 : 0, paneWidth = i ? width()-left : width()/2;
            p.save(); p.setClipRect(left, 0, paneWidth, height());
            auto mvp = camera(paneWidth);
            /** @brief 將正交投影座標映射到該側 Qt 畫面座標，並反轉 Y 方向。
             *         Map orthographic projected coordinates into this pane's Qt coordinates, reversing Y.
             */
            auto project = [&](QVector3D v) {
                QVector3D q = mvp.map(v);
                return QPointF(left + (q.x()+1)*paneWidth/2, (1-q.y())*height()/2);
            };
            if(paneGround[i] && cloud && cloud->volume.enabled && cloud->volume.error.isEmpty() && cloud->volume.groundValid) {
                const auto &v=cloud->volume;
                auto groundPoint=[&](double x,double y) {return project({float(x),float(y),float(v.groundZ(x,y)+.03)});};
                // 地面參考網格每 10 m 一線，僅為疊圖；計算格距仍固定 0.1 m。
                // The 10 m ground overlay is visual only; the integration raster remains fixed at 0.1 m.
                p.setPen(QPen(QColor(100,200,255,130),1,Qt::DashLine));
                const auto &region=cloud->appliedFill.footprint;
                if(!region.empty()) {
                    double x0=1e30,x1=-1e30,y0=1e30,y1=-1e30;
                    for(const auto &q:region) {x0=std::min(x0,q[0]);x1=std::max(x1,q[0]);y0=std::min(y0,q[1]);y1=std::max(y1,q[1]);}
                    // 逐段檢查凹邊界，地面格線不穿出場域。 / Clip individual ground-grid segments to the concave site.
                    for(double x=std::ceil(x0/10)*10;x<=x1;x+=10) for(double y=y0;y<y1;y+=1)
                        if(insideFootprint(region,x,y+.5)) p.drawLine(groundPoint(x,y),groundPoint(x,std::min(y+1,y1)));
                    for(double y=std::ceil(y0/10)*10;y<=y1;y+=10) for(double x=x0;x<x1;x+=1)
                        if(insideFootprint(region,x+.5,y)) p.drawLine(groundPoint(x,y),groundPoint(std::min(x+1,x1),y));
                }
                for(const auto &edge:v.boundaries) {
                    p.setPen(QPen(QColor::fromHsv(int(edge.pile*67%360),220,255),1.5));
                    p.drawLine(groundPoint(edge.x1,edge.y1),groundPoint(edge.x2,edge.y2));
                }
                for(const auto &pile:v.piles) {
                    p.setPen(Qt::white);p.drawText(project({float(pile.centerX),float(pile.centerY),float(v.groundZ(pile.centerX,pile.centerY)+pile.meanHeight)}),
                        QString("#%1  %2 m³").arg(qulonglong(pile.id)).arg(pile.volume,0,'f',1));
                }
            }
            if(paneColumns[i] && cloud) {
                p.setPen(QColor(255,70,220));p.drawText(left+12,44,cloud->columns.finalized?"Final column model; colors=groups, magenta=measurement evidence":"Column evidence overlay (magenta); pending=yellow, confirmed=green, rejected=gray");
                for(const auto &c:cloud->columns.items) {
                    const auto color=c.state==ColumnState::Confirmed?(cloud->columns.finalized&&!c.group.isEmpty()?QColor::fromHsv(int(qHash(c.group,0)*67u%360),180,255):QColor(40,255,100)):(c.state==ColumnState::Rejected?QColor(150,150,150):QColor(255,230,30));
                    p.setPen(QPen(color,1.5));
                    for(size_t j=0;j<c.footprint.size();++j) {const auto &a=c.footprint[j],&b=c.footprint[(j+1)%c.footprint.size()];
                        p.drawLine(project({float(a[0]),float(a[1]),float(columnBottomZ(c,a[0],a[1]))}),project({float(b[0]),float(b[1]),float(columnBottomZ(c,b[0],b[1]))}));
                        p.drawLine(project({float(a[0]),float(a[1]),float(columnTopZ(c,a[0],a[1]))}),project({float(b[0]),float(b[1]),float(columnTopZ(c,b[0],b[1]))}));
                        if(c.shape!="circle"||j%8==0) p.drawLine(project({float(a[0]),float(a[1]),float(columnBottomZ(c,a[0],a[1]))}),project({float(a[0]),float(a[1]),float(columnTopZ(c,a[0],a[1]))}));
                    }
                    const auto b=columnBounds(c);p.drawText(project({float((b[0]+b[1])*.5),float((b[2]+b[3])*.5),float(columnTopZ(c,(b[0]+b[1])*.5,(b[2]+b[3])*.5))}),QString("%1C%2%3").arg(c.group.isEmpty()?QString():c.group+"/").arg(c.id).arg(c.provenance.startsWith("pattern_inferred")?"*":""));
                }
            }
            if(showAxes) {
                const QPointF origin=project({0,0,0});const QVector3D ends[]{{axisSize,0,0},{0,axisSize,0},{0,0,axisSize}};
                const QColor colors[]{QColor(255,70,70),QColor(70,255,70),QColor(80,140,255)};const char *names[]{"X","Y","Z"};
                for(int a=0;a<3;++a) {p.setPen(QPen(colors[a],2));const auto end=project(ends[a]);p.drawLine(origin,end);p.drawText(end+QPointF(4,-4),names[a]);}
                p.setPen(Qt::white);p.drawEllipse(origin,3,3);p.drawText(origin+QPointF(5,15),"O (0, 0, 0)");
            }
            const QStringList names{"原始點雲","處理後料面","清理後實測","補點","柱體量測證據","不顯示點雲"};
            int displayed=0;for(int layer=0;layer<4;++layer) if(visibleLayer(paneModes[i],layer)) displayed+=counts[layer];
            QString title=names[paneModes[i]]+QString(" | 抽樣 %1 點").arg(displayed);
            if(i==0&&hideLeftRoof&&cloud&&cloud->hasRoofPlanes) title+="（隱藏屋頂與構架）";
            p.fillRect(QRect(left+5,5,paneWidth-10,27), QColor(10,15,20,210));
            p.setPen(Qt::white);p.drawText(left+12,24,title); p.restore();
        }
        if (draggingROI) {
            p.setPen(QPen(Qt::yellow,2,Qt::DashLine)); p.setBrush(QColor(255,255,0,30));
            p.drawRect(QRect(roiStart,roiEnd).normalized());
        }
        /** @brief 與 GPU 相同的藍到紅線性色階，超出範圍時使用端點色。
         *         The same linear blue-to-red scale as the GPU; out-of-range values use endpoint colors.
         */
        if (cloud && colorMode != 0) {
            QRect bar(20,height()-46,240,12);
            QLinearGradient gradient(bar.topLeft(),bar.topRight());
            gradient.setColorAt(0,QColor::fromRgbF(.1,.25,1));
            gradient.setColorAt(1,QColor::fromRgbF(1,.15,.05));
            p.fillRect(bar,gradient); p.setPen(Qt::white);
            p.drawText(20,height()-55,colorMode==1 ? "Intensity (missing: gray)" : "Z height (m)");
            p.drawText(20,height()-16,QString::number(colorMin,'g',7));
            p.drawText(180,height()-16,QString::number(colorMax,'g',7));
        }
    }
    /** @brief 記錄拖曳起點。 / Record the drag start. @param event Qt 滑鼠事件。 / Qt mouse event. */
    void mousePressEvent(QMouseEvent *event) override {
        last = event->pos();
        if (selectingROI && event->button()==Qt::LeftButton) {
            roiStart=roiEnd=last; roiPane=last.x()>=width()/2 ? 1 : 0; draggingROI=true; update();
        }
    }
    /**
     * @brief 依按鍵更新共享旋轉或平移。 / Update shared rotation or pan according to held buttons.
     * @param event Qt 滑鼠事件，含位置、按鍵與修飾鍵。 / Qt mouse event with position, buttons, and modifiers.
     */
    void mouseMoveEvent(QMouseEvent *event) override {
        if (draggingROI) {
            const int left=roiPane ? width()/2 : 0, right=roiPane ? width()-1 : width()/2-1;
            roiEnd={std::clamp(event->pos().x(),left,right),std::clamp(event->pos().y(),0,height()-1)};
            update(); return;
        }
        QPoint delta = event->pos() - last; last = event->pos();
        if ((event->buttons() & Qt::RightButton) || (event->modifiers() & Qt::ShiftModifier)) {
            panX += delta.x() * 2 * radius * zoom / std::max(1,height());
            panY -= delta.y() * 2 * radius * zoom / std::max(1,height());
        } else if (event->buttons() & Qt::LeftButton) {
            yaw += delta.x()*.4f; pitch += delta.y()*.4f;
        }
        update();
    }
    /** @brief 俯視框選轉換成世界 XY 邊界；Z 由使用者另行限制。
     *         Convert a top-view rectangle into world XY bounds; users constrain Z separately.
     * @param event Qt 滑鼠釋放事件。 / Qt mouse-release event.
     */
    void mouseReleaseEvent(QMouseEvent *event) override {
        if (!draggingROI || event->button()!=Qt::LeftButton) return;
        draggingROI=false; selectingROI=false; unsetCursor();
        int left=roiPane ? width()/2 : 0, paneWidth=roiPane ? width()-left : width()/2;
        auto inverse=camera(paneWidth).inverted();
        auto world=[&](QPoint p) {return inverse.map(QVector3D(2.f*(p.x()-left)/paneWidth-1,1-2.f*p.y()/height(),0));};
        auto a=world(roiStart),b=world(roiEnd);
        if (roiCallback && (roiStart-roiEnd).manhattanLength()>5)
            roiCallback(std::min(a.x(),b.x()),std::max(a.x(),b.x()),std::min(a.y(),b.y()),std::max(a.y(),b.y()));
        update();
    }
    /** @brief 以滾輪等比縮放並限制範圍。 / Apply bounded multiplicative wheel zoom.
     * @param event 含滾輪角度增量的事件。 / Event containing the wheel angle delta.
     */
    void wheelEvent(QWheelEvent *event) override {
        zoom = std::clamp(zoom * float(std::pow(.85, event->angleDelta().y()/120.)), .0001f, 100.f); update();
    }
    void mouseDoubleClickEvent(QMouseEvent *) override {reset();}
public:
    /** @brief 進入 Z 朝上的俯視框選模式。 / Enter Z-up top-view ROI selection mode.
     * @param callback 接收所選 XY 邊界。 / Receives selected XY bounds.
     */
    void selectXY(std::function<void(double,double,double,double)> callback) {
        roiCallback=std::move(callback); yaw=pitch=0; selectingROI=true; setCursor(Qt::CrossCursor); update();
    }
    /** @brief 切換地平面與料堆邊界疊圖，不影響體積。 / Toggle ground/boundary overlays without changing volume.
     * @param enabled 是否顯示疊圖。 / Whether to show overlays.
     */
    /** @brief 顯示候選原量測證據及語意截面，不混入料堆資料。 / Show proposal evidence and sections without adding them to material data. */
    void setPaneDisplay(int pane,int mode,bool columns,bool ground) {if(pane<0||pane>1||mode<0||mode>5) return;paneModes[pane]=mode;paneColumns[pane]=columns;paneGround[pane]=ground;update();}
    void setLeftRoofDepth(float depth) {roofHideDepth=std::clamp(depth,.1f,15.f);update();}
    void setLeftRoofHidden(bool hidden) {hideLeftRoof=hidden;update();}
    void setAxesVisible(bool visible) {showAxes=visible;update();}
    /** @brief 放大選取柱體並在左右視窗比較原始結構與語意證據。 / Focus a selected column in synchronized raw/evidence panes. */
    void focusColumn(const ColumnSemantic &c,bool topView=false) {
        const auto b=columnBounds(c);const double bottom=columnBottomZ(c,(b[0]+b[1])*.5,(b[2]+b[3])*.5),top=columnTopZ(c,(b[0]+b[1])*.5,(b[2]+b[3])*.5);
        center=QVector3D(float((b[0]+b[1])*.5),float((b[2]+b[3])*.5),float((bottom+top)*.5));
        const double extent=topView?std::max({b[1]-b[0],b[3]-b[2],2.})*1.5:std::max({top-bottom,b[1]-b[0],b[3]-b[2],2.});
        zoom=std::clamp(float(extent/std::max(radius,.01f)),.0001f,100.f);panX=panY=0;yaw=topView?0:-35;pitch=topView?0:20;update();
    }
    /** @brief 切換補點橘色標示。 / Toggle orange synthetic-point highlighting. */
    void setFillHighlight(bool enabled) { highlightFill=enabled; update(); }
    /** @brief 只在畫面裁切 Z 範圍以便排除屋頂遮擋，不刪除量測點。
     *         Clip the displayed Z range to see beneath roofs without deleting measurements.
     * @param enabled 啟用預覽裁切。 / Enable preview clipping.
     * @param lower Z 下限（公尺）。 / Lower Z limit in meters.
     * @param upper Z 上限（公尺）。 / Upper Z limit in meters.
     */
    void setZSlice(bool enabled,float lower,float upper) { zSlice=enabled; sliceMin=lower; sliceMax=upper; update(); }
    /** @brief 建立檢視器並設定最小尺寸。 / Construct the viewer and set its minimum size.
     * @param parent Qt 父元件，可為空。 / Optional Qt parent widget.
     */
    explicit Viewer(QWidget *parent = nullptr) : QOpenGLWidget(parent) { setMinimumSize(640,400); }
    /** @brief 解除 context 回呼，避免物件解構後再次呼叫清理。 / Disconnect context callbacks before destruction and release resources. */
    ~Viewer() override {
        if (context()) disconnect(context(),nullptr,this,nullptr);
        cleanup();
    }
    /** @return GL 資源是否完成初始化。 / Whether GL resource initialization succeeded. */
    bool graphicsReady() const { return ready; }
    /**
     * @brief 只更新 GPU 色階參數，不重讀檔案或重新上傳點資料。 / Update GPU color uniforms without reloading or re-uploading points.
     * @param mode 0 單色、1 強度、2 高度。 / 0 uniform, 1 intensity, 2 height.
     * @param minimum 色階下限。 / Lower color limit.
     * @param maximum 色階上限；與下限相同時使用中間色。 / Upper color limit; equal limits use the midpoint color.
     */
    void setColorScale(int mode, float minimum, float maximum) {
        colorMode = mode; colorMin = minimum; colorMax = maximum; update();
    }
    /**
     * @brief 置換顯示資料並重設兩側相機。 / Replace display data and reset the shared camera.
     * @param value 已成功處理、不可為空的點雲結果。 / Non-null, successfully processed cloud result.
     * @param axis 座標軸長度，單位公尺。 / Axis length in meters.
     */
    void clearCloud() {cloud.reset();for(auto &count:counts) count=0;dirty=false;update();}
    void setCloud(std::shared_ptr<Cloud> value, float axis) {
        cloud = std::move(value); axisSize = axis; dirty = true; reset();
    }
    /** @brief 依點雲、原點與軸端重新取景並清除互動位移。 / Reframe the cloud, origin, and axis endpoints and clear interactive offsets. */
    void reset() {
        selectingROI=draggingROI=false; unsetCursor();
        if (cloud) {
            QVector3D lo(cloud->low.x,cloud->low.y,cloud->low.z), hi(cloud->high.x,cloud->high.y,cloud->high.z);
            hi.setX(std::max(hi.x(),axisSize)); hi.setY(std::max(hi.y(),axisSize)); hi.setZ(std::max(hi.z(),axisSize));
            center = (lo+hi)*.5f; radius = std::max(1.f,(hi-lo).length()*.65f);
        }
        zoom = 1; yaw = -35; pitch = 25; panX = panY = 0; update();
    }
};

/**
 * @brief 解析命令列、執行批次檢查或啟動 Qt 視覺化。 / Parse arguments and run batch checks or Qt visualization.
 * @param argc 命令列參數數量。 / Command-line argument count.
 * @param argv 命令列參數字串。 / Command-line argument strings.
 * @return 0 表示成功；非零表示參數、處理或視窗測試失敗。 / Zero on success; nonzero on argument, processing, or GUI-test failure.
 */
int main(int argc, char **argv) {
    /** @brief 先辨識 --check，讓無桌面模式使用 QCoreApplication，不建立 GUI/OpenGL context。
     *         Detect --check early so headless mode uses QCoreApplication without a GUI/OpenGL context.
     */
    bool check = false;
    for (int i=1;i<argc;++i) if (QString::fromLocal8Bit(argv[i]) == "--check") check = true;
    std::unique_ptr<QCoreApplication> app;
    QSurfaceFormat format; format.setVersion(3,3); format.setProfile(QSurfaceFormat::CoreProfile);
    format.setDepthBufferSize(24); format.setSamples(0); QSurfaceFormat::setDefaultFormat(format);
    if (check) app = std::make_unique<QCoreApplication>(argc,argv);
    else app = std::make_unique<QApplication>(argc,argv);
    QCoreApplication::setApplicationName("pointcloud_compare");
    /** @brief 註冊並驗證命令列選項；顯示上限獨立於完整量測計算。
     *         Register and validate command-line options; the display cap is independent of full-cloud processing.
     */
    QCommandLineParser parser; parser.setApplicationDescription("Qt5 PLY comparison; coordinates must be in meters");
    parser.addHelpOption(); parser.addPositionalArgument("directory", "PLY directory (default: current directory)", "[directory]");
    parser.addOption({"max-display-points", "Maximum displayed points per pane, 0 = all (default: 1000000)", "count", "1000000"});
    parser.addOption({"axis-size", "XYZ axis length in meters", "meters", "1"});
    parser.addOption({"output-dir", "Export each full processed cloud as PLY (requires --check)", "directory"});
    parser.addOption({"check", "Process every PLY and print counts without opening a window"});
    parser.addOption({"raw-coordinates", "Ignore calibration.json and inspect original scan coordinates"});
    parser.addOption({"no-columns", "Disable column proposals and reviewed semantic exclusions"});
    parser.addOption({"column-min-height", "Minimum persistent vertical span for column proposals", "meters", "4"});
    parser.addOption({"column-coverage", "Minimum occupied fraction of vertical bins", "fraction", "0.65"});
    parser.addOption({"column-max-width", "Maximum proposal cross-section width", "meters", "2.5"});
    parser.addOption({"no-surface-filter", "Disable final pre-volume local surface outlier filtering"});
    parser.addOption({"surface-outlier-tolerance", "Minimum local surface residual threshold in meters", "meters", "0.6"});
    parser.addOption({"no-mesh-volume", "Disable experimental mesh integration while retaining raster volume"});
    parser.addOption({"mesh-support-spacing", "Mesh support spacing, independent of volume raster", "meters", "0.25"});
    parser.addOption({"mesh-max-edge", "Maximum XY mesh edge/unsupported gap", "meters", "2.0"});
    parser.addOption({"mesh-max-slope", "Maximum mesh surface gradient", "ratio", "3.0"});
    parser.addOption({"no-volume", "Disable stockpile ground/volume analysis"});
    parser.addOption({"ground-plane", "Manual ground z=a*x+b*y+c (meters): a,b,c", "coefficients"});
    parser.addOption({"ground-tolerance", "Ground vertical inlier tolerance in meters", "meters", "0.15"});
    parser.addOption({"pile-threshold", "Minimum height above ground classified as material", "meters", "0.2"});
    parser.addOption({"pile-min-area", "Minimum connected pile area in square meters", "square-meters", "1.0"});
    parser.addOption({"no-floating-cleanup", "Disable detached/sparse cluster cleanup after structure removal"});
    parser.addOption({"floating-cell", "Post-structure 3D connectivity cell width in meters", "meters", "0.5"});
    parser.addOption({"floating-min-points", "Minimum measured points per post-structure cell", "count", "2"});
    parser.addOption({"floating-keep-all-components", "Keep detached dense components after structure removal"});
    parser.addOption({"no-structures", "Disable building structure removal"});
    parser.addOption({"wall-margin", "Interior wall setback in meters", "meters", "0.4"});
    parser.addOption({"roof-margin", "Roof exclusion distance in meters", "meters", "2.0"});
    parser.addOption({"surface-band", "Retained thickness above lowest XY-cell measurement", "meters", "0.8"});
    parser.addOption({"structure-cell", "Lower-surface XY cell width in meters", "meters", "0.5"});
    parser.addOption({"no-denoise", "Disable initial measurement denoising"});
    parser.addOption({"noise-cell-size", "Denoising cell size in meters, >= 0.1", "meters", "0.5"});
    parser.addOption({"noise-min-points", "Minimum raw measured points in a denoising cell", "count", "3"});
    parser.addOption({"keep-all-components", "Retain every sufficiently dense component"});
    parser.addOption({"fill-roi", "Enable stockpile filling inside xMin,xMax,yMin,yMax,zMin,zMax (meters)", "bounds"});
    parser.addOption({"no-fill", "Disable hole filling"});
    parser.addOption({"fill-spacing", "Synthetic XY point spacing in meters", "meters", "0.1"});
    parser.addOption({"fill-max-hole", "Maximum hole diagonal / TIN edge length in meters", "meters", "20.0"});
    parser.addOption({"fill-max-error", "Boundary fit/continuation tolerance in meters", "meters", "1.0"});
    parser.addOption({"fill-max-slope", "Maximum stockpile slope dz / horizontal distance", "slope", "3.0"});
    parser.addOption({"fill-method", "Filling method: tin (layered Delaunay) or plane (legacy)", "method", "tin"});
    parser.addOption({"fill-layer-height", "TIN overlapping height-band thickness in meters", "meters", "8"});
    parser.addOption({"fill-support-spacing", "TIN measured support spacing in meters", "meters", "0.5"});
    parser.addOption({"validate-fill", "Withhold up to 3 known circular patches and report prediction coverage/error"});
    parser.addOption({"smoke-test", "Close after first rendered dataset (GUI verification)"});
    parser.process(*app);
    bool valid; qulonglong limit = parser.value("max-display-points").toULongLong(&valid);
    if (!valid || limit > 100000000) { std::cerr << "Invalid display limit (0..100000000)\n"; return 1; }
    float axis = parser.value("axis-size").toFloat(&valid);
    if (!valid || !std::isfinite(axis) || axis <= 0) { std::cerr << "Invalid axis size\n"; return 1; }
    DenoiseOptions noise;
    noise.enabled = !parser.isSet("no-denoise");
    noise.cellSize = parser.value("noise-cell-size").toDouble(&valid);
    if (!valid || !std::isfinite(noise.cellSize) || noise.cellSize < .1 || noise.cellSize > 100) {
        std::cerr << "Invalid denoising cell size (0.1..100 meters)\n"; return 1;
    }
    noise.minPoints = parser.value("noise-min-points").toUInt(&valid);
    if (!valid || !noise.minPoints || noise.minPoints > 1000000) { std::cerr << "Invalid noise minimum points (1..1000000)\n"; return 1; }
    noise.largestComponent = !parser.isSet("keep-all-components");
    StructureOptions structures;
    structures.enabled=!parser.isSet("no-structures") && !parser.isSet("raw-coordinates");
    for(auto setting : {std::make_pair("wall-margin",&structures.wallMargin),std::make_pair("roof-margin",&structures.roofMargin),
                       std::make_pair("surface-band",&structures.surfaceBand),std::make_pair("structure-cell",&structures.cellSize)}) {
        bool valid=false;*setting.second=parser.value(setting.first).toDouble(&valid);
        if(!valid || !std::isfinite(*setting.second) || *setting.second<0 ||
           (QString(setting.first)=="surface-band" && *setting.second<=0) ||
           (QString(setting.first)=="structure-cell" && *setting.second<.1)) {std::cerr<<"Invalid structure parameter\n";return 1;}
    }
    structures.floating.enabled=!parser.isSet("no-floating-cleanup");
    structures.floating.largestComponent=!parser.isSet("floating-keep-all-components");
    bool floatingValid=false;
    structures.floating.cellSize=parser.value("floating-cell").toDouble(&floatingValid);
    if(!floatingValid || !std::isfinite(structures.floating.cellSize) || structures.floating.cellSize<.1 || structures.floating.cellSize>100) {
        std::cerr<<"Invalid floating cell size (0.1..100 m)\n";return 1;
    }
    structures.floating.minPoints=parser.value("floating-min-points").toULongLong(&floatingValid);
    if(!floatingValid || !structures.floating.minPoints || structures.floating.minPoints>1000000) {std::cerr<<"Invalid floating minimum points (1..1000000)\n";return 1;}
    structures.columns.enabled=!parser.isSet("no-columns")&&structures.enabled;
    for(auto item:{std::make_pair("column-min-height",&structures.columns.minHeight),std::make_pair("column-coverage",&structures.columns.minCoverage),std::make_pair("column-max-width",&structures.columns.maxWidth)}) {
        bool ok=false;*item.second=parser.value(item.first).toDouble(&ok);if(!ok||!std::isfinite(*item.second)||*item.second<=0) {std::cerr<<"Invalid column parameter\n";return 1;}
    }
    structures.surfaceFilter.enabled=!parser.isSet("no-surface-filter");
    bool surfaceValid=false;structures.surfaceFilter.tolerance=parser.value("surface-outlier-tolerance").toDouble(&surfaceValid);
    if(!surfaceValid||!std::isfinite(structures.surfaceFilter.tolerance)||structures.surfaceFilter.tolerance<=0) {std::cerr<<"Invalid surface filter tolerance\n";return 1;}
    VolumeOptions volume;volume.enabled=!parser.isSet("no-volume") && structures.enabled;
    volume.mesh.enabled=!parser.isSet("no-mesh-volume");
    for(auto item:{std::make_pair("mesh-support-spacing",&volume.mesh.supportSpacing),std::make_pair("mesh-max-edge",&volume.mesh.maxEdge),std::make_pair("mesh-max-slope",&volume.mesh.maxSlope)}) {
        bool ok=false;*item.second=parser.value(item.first).toDouble(&ok);if(!ok||!std::isfinite(*item.second)||*item.second<=0) {std::cerr<<"Invalid mesh parameter\n";return 1;}
    }
    if(volume.mesh.supportSpacing<.05||volume.mesh.maxEdge<volume.mesh.supportSpacing) {std::cerr<<"Invalid mesh support/edge spacing\n";return 1;}
    if(parser.isSet("ground-plane")) {
        const auto coefficients=parser.value("ground-plane").split(',');
        if(coefficients.size()!=3) {std::cerr<<"Ground plane requires a,b,c\n";return 1;}
        for(int i=0;i<3;++i) {bool ok=false;volume.ground[i]=coefficients[i].toDouble(&ok);if(!ok||!std::isfinite(volume.ground[i])) {std::cerr<<"Invalid ground plane\n";return 1;}}
        volume.automaticGround=false;
    }
    for(auto option:{std::make_pair("ground-tolerance",&volume.groundTolerance),std::make_pair("pile-threshold",&volume.pileThreshold),std::make_pair("pile-min-area",&volume.minimumArea)}) {
        bool ok=false;*option.second=parser.value(option.first).toDouble(&ok);
        if(!ok||!std::isfinite(*option.second)||*option.second<0||(QString(option.first)=="ground-tolerance"&&*option.second==0)) {std::cerr<<"Invalid volume parameter\n";return 1;}
    }
    FillOptions fill;
    if (parser.value("fill-method")!="tin" && parser.value("fill-method")!="plane") {std::cerr << "Invalid fill method\n"; return 1;}
    fill.method=parser.value("fill-method")=="tin" ? FillMethod::LayeredTin : FillMethod::LocalPlane;
    fill.validate=parser.isSet("validate-fill");
    fill.enabled=!parser.isSet("no-fill");
    fill.fullDomain=!parser.isSet("fill-roi");
    if (parser.isSet("fill-roi")) {
        const auto values=parser.value("fill-roi").split(',');
        if (values.size()!=6) {std::cerr << "Fill ROI requires xMin,xMax,yMin,yMax,zMin,zMax\n"; return 1;}
        double *bounds[]{&fill.xMin,&fill.xMax,&fill.yMin,&fill.yMax,&fill.zMin,&fill.zMax};
        for (int i=0;i<6;++i) {
            *bounds[i]=values[i].toDouble(&valid);
            if (!valid || !std::isfinite(*bounds[i])) {std::cerr << "Invalid fill ROI\n"; return 1;}
        }
    }
    for (auto option : {std::make_pair("fill-spacing",&fill.spacing),std::make_pair("fill-max-hole",&fill.maxHole),
                        std::make_pair("fill-max-error",&fill.maxError),std::make_pair("fill-max-slope",&fill.maxSlope),
                        std::make_pair("fill-layer-height",&fill.layerHeight),std::make_pair("fill-support-spacing",&fill.supportSpacing)}) {
        *option.second=parser.value(option.first).toDouble(&valid);
        if (!valid || !std::isfinite(*option.second) || *option.second<=0) {std::cerr << "Invalid filling parameter\n"; return 1;}
    }
    auto positional = parser.positionalArguments();
    if (positional.size()>1) parser.showHelp(1);
    /** @brief 以檔名排序列出目錄中的 PLY，不遞迴子目錄。
     *         List PLY files in filename order without recursing into subdirectories.
     */
    QDir dir(positional.isEmpty() ? QDir::currentPath() : positional[0]);
    QStringList files;
    for (const auto &name : dir.entryList(QDir::Files, QDir::Name)) if (name.endsWith(".ply",Qt::CaseInsensitive)) files << name;
    if (!dir.exists() || files.isEmpty()) { std::cerr << "No PLY files in directory\n"; return 1; }
    /** @brief 無視窗路徑：逐檔處理並輸出完整點數、樣本點數與錯誤。
     *         Headless path: process each file and report full counts, sample counts, and errors.
     */
    if(parser.isSet("output-dir") && !check) {std::cerr<<"Use --check with --output-dir, or the Export PLY button\n";return 1;}
    QString outputDirectory;
    if(parser.isSet("output-dir")) {
        outputDirectory=QFileInfo(parser.value("output-dir")).absoluteFilePath();
        if(QDir(outputDirectory).canonicalPath()==dir.canonicalPath()) {std::cerr<<"Choose an output directory different from the input directory\n";return 1;}
        if(!QDir().mkpath(outputDirectory)) {std::cerr<<"Cannot create output directory\n";return 1;}
    }
    if (check) {
        int status = 0;
        for (const auto &file : files) {
            const QString output=outputDirectory.isEmpty()?QString{}:QDir(outputDirectory).filePath(QFileInfo(file).completeBaseName()+"_processed_filled.ply");
            if(!output.isEmpty() && QFileInfo::exists(output)) {std::cerr<<"Output exists; skipped: "<<output.toStdString()<<"\n";status=1;continue;}
            auto c = processCloud(dir.filePath(file),limit,noise,fill,!parser.isSet("raw-coordinates"),structures,volume,output);
            if(!c.exportedPath.isEmpty()) std::cout<<"Exported full processed PLY: "<<c.exportedPath.toStdString()<<"\n";
            if(!c.exportError.isEmpty()) {std::cerr<<"PLY export failed: "<<c.exportError.toStdString()<<"\n";status=1;}
            std::cout<<"Coordinates: "<<(c.coordinateFrame.isEmpty()?"source":c.coordinateFrame.toStdString())<<"\n";
            std::cout << file.toStdString() << ": " << c.originalCount << " -> " << c.reducedCount
                      << " | displayed " << c.before.size() << "/" << c.after.size()+c.filled.size() << "\n";
            std::cout << "Pre-volume surface filter removed measured=" << c.surfaceFilter.removedMeasured << " estimated=" << c.surfaceFilter.removedEstimated << "\n";
            if(c.surfaceFilter.comparisonValid) std::cout << "Same-ground volume before="<<c.surfaceFilter.volumeBefore<<" after="<<c.surfaceFilter.volumeAfter<<" delta="<<c.surfaceFilter.volumeAfter-c.surfaceFilter.volumeBefore<<" m3\n";
            std::cout << "Input measurements: " << c.originalCount << " | noise removed: " << c.originalCount-c.reducedCount-c.structureRemoved-c.floatingRemoved-c.columns.removedMeasured-c.surfaceFilter.removedMeasured << " | measured structures removed: " << c.structureRemoved << " | floating/sparse removed: " << c.floatingRemoved << "\n";
            if(!c.footprintSource.isEmpty()) std::cout<<"Footprint: "<<c.footprintSource.toStdString()<<" | restored end measurements: "<<c.restoredEnds[0]<<" / "<<c.restoredEnds[1]<<"\n";
            std::cout << "Estimated fill: " << c.filledCount << " points in " << c.filledHoles
                << (fill.method==FillMethod::LayeredTin ? " contributing triangles\n" : " holes\n");
            if (fill.method==FillMethod::LayeredTin) {
                const auto &r=c.fillReport;
                std::cout << "TIN bands: " << r.bands << " | supports: " << r.supports << " | accepted triangles: " << r.triangles
                    << " | rejected distance/slope/consistency: " << r.rejectedDistance << "/" << r.rejectedSlope << "/" << r.rejectedConsistency
                    << " | rejected footprint: " << r.rejectedFootprint << " | output limit: " << r.outputLimited << "\n";
                if (fill.validate) {
                    std::cout << "Holdout patches: " << r.validationPatches << " | predicted/withheld: " << r.validationPredicted << "/" << r.validationTargets;
                    if (r.validationPredicted) std::cout << " | RMSE: " << std::sqrt(r.validationSquaredError/r.validationPredicted)
                        << " m | MAE: " << r.validationAbsoluteError/r.validationPredicted << " m | max error: " << r.validationMaxError << " m";
                    std::cout << "\n";
                    for (const auto &patch : r.validation) {
                        std::cout << "  Patch center " << patch.center.x << "," << patch.center.y << "," << patch.center.z
                            << " | radius " << patch.radius << " m | coverage " << patch.predicted << "/" << patch.targets;
                        if (patch.predicted) std::cout << " | RMSE " << patch.rmse << " m | max " << patch.maxError << " m\n";
                        else std::cout << " | Error: N/A (no prediction)\n";
                    }
                }
            }
            if(structures.columns.enabled) {
                std::cout<<"Column semantics: "<<c.columns.modelPath.toStdString()<<" | proposals/reviews "<<c.columns.items.size()<<" | removed measured/estimated "<<c.columns.removedMeasured<<"/"<<c.columns.removedEstimated<<"\n";
                for(const auto &column:c.columns.items) {const auto b=columnBounds(column);std::cout<<"  C"<<column.id<<" "<<(column.state==ColumnState::Confirmed?"confirmed":column.state==ColumnState::Rejected?"rejected":"pending")
                    <<" | XY "<<b[0]<<","<<b[1]<<","<<b[2]<<","<<b[3]<<" | observed Z "<<column.observedBottom<<".."<<column.observedTop<<" | coverage "<<column.verticalCoverage<<" | evidence "<<column.evidencePoints<<" | excluded "<<column.excludedVolume<<" m3\n";}
            }
            if(c.volume.enabled) {
                const auto &v=c.volume;
                if(!v.error.isEmpty()) {std::cerr<<"Volume unavailable: "<<v.error.toStdString()<<"\n";status=1;}
                else {
                    std::cout<<"Ground z = "<<v.ground[0]<<"*x + "<<v.ground[1]<<"*y + "<<v.ground[2]<<" | "<<(v.manualGround?"manual":"estimated")
                        <<" | measured inliers "<<v.groundInliers<<"/"<<v.groundSamples;
                    if(v.manualGround) std::cout<<" | RMSE N/A (manual)\n";else std::cout<<" | RMSE "<<v.groundRMSE<<" m\n";
                    std::cout<<"Volume raster: 0.10 x 0.10 m (minimum height) | covered "<<v.coveredCells<<"/"<<v.footprintCells<<" | unknown area "<<v.unknownCells*.01
                        <<" m2 | below-ground area "<<v.belowGroundCells*.01<<" m2 | piles "<<v.piles.size()<<" | supported volume "<<v.totalVolume<<" m3 | column intersection excluded "<<v.columnExcludedVolume<<" m3\n";
                    for(const auto &pile:v.piles) std::cout<<"  Pile "<<pile.id<<" | area "<<pile.area<<" m2 | mean height "<<pile.meanHeight<<" m | volume "<<pile.volume
                        <<" m3 | measured/fill volume "<<pile.measuredVolume<<"/"<<pile.estimatedVolume<<" m3 | enclosed unknown area "<<pile.enclosedUnknownCells*.01<<" m2\n";
                }
            }
            if(c.meshVolume.enabled) {
                const auto &m=c.meshVolume;
                if(!m.error.isEmpty()) {std::cerr<<"Mesh volume unavailable: "<<m.error.toStdString()<<"\n";status=1;}
                else std::cout<<"Mesh volume (experimental): "<<m.totalVolume<<" m3 | gross "<<m.grossVolume<<" | column excluded "<<m.columnExcludedVolume<<" | small components excluded "<<m.discardedSmallVolume<<" | surface area XY "<<m.coveredArea<<" / "<<m.footprintArea<<" | supports "<<m.vertices<<" at "<<m.supportSpacing<<" m | triangles "<<m.surface.size()<<" | mesh-raster "<<m.totalVolume-c.volume.totalVolume<<" m3\n";
            }
            std::cout << "Fill grid: " << c.appliedFill.spacing << " m | XY: " << c.appliedFill.xMin << "," << c.appliedFill.xMax
                << "," << c.appliedFill.yMin << "," << c.appliedFill.yMax << "\n";
            if (c.error.isEmpty()) std::cout << "Z range: " << c.zMin << " .. " << c.zMax
                << " | Intensity: " << (c.hasIntensity ? "available" : "missing")
                << " | range: " << c.intensityMin << " .. " << c.intensityMax << "\n";
            if (c.error.isEmpty()) std::cout << "Auto P5-P95 | Z: " << c.zP05 << " .. " << c.zP95
                << " | Intensity: " << c.intensityP05 << " .. " << c.intensityP95 << "\n";
            if (!c.error.isEmpty()) { std::cerr << c.error.toStdString() << "\n"; status = 1; }
        }
        return status;
    }
    /** @brief GUI 路徑：建立檔案選單、共享相機檢視器、狀態列與操作提示。
     *         GUI path: create the file selector, shared-camera viewer, status text, and controls hint.
     */
    /** @brief 精簡主畫面只管理顯示、匯出與體積報告；處理設定由命令列固定。
     *         The compact window controls display/export and reports volume; processing options come from the CLI.
     */
    QWidget window;window.setWindowTitle("料場點雲與體積");
    auto *layout=new QVBoxLayout(&window);
    auto *fileRow=new QHBoxLayout;auto *select=new QComboBox;select->addItems(files);
    auto *exportPly=new QPushButton("匯出處理後 PLY");exportPly->setEnabled(false);
    exportPly->setToolTip("匯出完整處理後實測與補點；與左右顯示選擇無關");
    auto *exportMesh=new QPushButton("匯出料面 Mesh");exportMesh->setEnabled(false);exportMesh->setToolTip("匯出實驗性 TIN 的開放上表面供檢查；柱體以解析交集另扣，非封閉實體");
    fileRow->addWidget(exportMesh);fileRow->addWidget(new QLabel("點雲檔案"));fileRow->addWidget(select,1);fileRow->addWidget(exportPly);layout->addLayout(fileRow);
    auto *displayRow=new QHBoxLayout;
    QComboBox *paneMode[2];QCheckBox *paneColumns[2],*paneGround[2];
    for(int pane=0;pane<2;++pane) {
        displayRow->addWidget(new QLabel(pane?"右側顯示":"左側顯示"));paneMode[pane]=new QComboBox;
        paneMode[pane]->addItems({"原始點雲","處理後料面（實測＋補點）","清理後實測","補點","柱體量測證據","不顯示點雲"});paneMode[pane]->setCurrentIndex(pane?1:0);
        paneColumns[pane]=new QCheckBox("柱體框");paneGround[pane]=new QCheckBox("地面與料堆邊界");paneGround[pane]->setChecked(true);
        displayRow->addWidget(paneMode[pane]);displayRow->addWidget(paneColumns[pane]);displayRow->addWidget(paneGround[pane]);if(!pane) displayRow->addSpacing(20);
    }
    auto *hideRoof=new QCheckBox("左側隱藏屋頂與構架");hideRoof->setEnabled(false);
    hideRoof->setToolTip("隱藏屋頂與屋頂下方指定深度內的點，只影響左側；可調整深度以露出料堆");
    auto *roofDepth=new QDoubleSpinBox;roofDepth->setRange(.1,15);roofDepth->setDecimals(1);roofDepth->setSingleStep(.5);roofDepth->setValue(3);roofDepth->setSuffix(" m");roofDepth->setKeyboardTracking(false);roofDepth->setEnabled(false);
    roofDepth->setToolTip("由斜屋頂垂直向下的顯示裁切深度；增加可隱藏較低構架，也可能遮住較高料面");
    displayRow->insertWidget(4,hideRoof);displayRow->insertWidget(5,roofDepth);displayRow->addStretch();layout->addLayout(displayRow);
    auto *colorRow=new QHBoxLayout;auto *colorSelect=new QComboBox;colorSelect->addItems({"單色","強度","Z 高度"});colorSelect->setCurrentIndex(2);
    auto *minimum=new QDoubleSpinBox,*maximum=new QDoubleSpinBox;
    for(auto *spin:{minimum,maximum}) {spin->setRange(-double(std::numeric_limits<float>::max()),double(std::numeric_limits<float>::max()));spin->setDecimals(3);spin->setKeyboardTracking(false);}
    auto *autoRange=new QPushButton("自動色階 P5–P95");auto *highlight=new QCheckBox("補點標橘色");highlight->setChecked(true);
    auto *axes=new QCheckBox("座標原點與軸");axes->setChecked(true);
    colorRow->addWidget(new QLabel("共用著色"));colorRow->addWidget(colorSelect);colorRow->addWidget(new QLabel("下限"));colorRow->addWidget(minimum);colorRow->addWidget(new QLabel("上限"));colorRow->addWidget(maximum);colorRow->addWidget(autoRange);colorRow->addWidget(highlight);colorRow->addWidget(axes);colorRow->addStretch();layout->addLayout(colorRow);
    auto *viewer=new Viewer;viewer->setToolTip("左鍵旋轉；右鍵或 Shift 拖曳平移；滾輪縮放；雙擊重新取景。左右視角同步。");layout->addWidget(viewer,1);
    auto *status=new QLabel;status->setWordWrap(true);layout->addWidget(status);
    auto *volumeSummary=new QLabel;volumeSummary->setWordWrap(true);volumeSummary->setTextInteractionFlags(Qt::TextSelectableByMouse);layout->addWidget(volumeSummary);
    auto *volumeTable=new QTableWidget(0,6);volumeTable->setHorizontalHeaderLabels({"柵格料堆","淨體積 m³","投影面積 m²","實測支撐 m³","補點支撐 m³","內部缺測 m²"});
    volumeTable->setEditTriggers(QAbstractItemView::NoEditTriggers);volumeTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);volumeTable->setMaximumHeight(145);layout->addWidget(volumeTable);
    std::shared_ptr<Cloud> activeCloud;QFutureWatcher<std::shared_ptr<Cloud>> watcher;
    QObject::connect(hideRoof,&QCheckBox::toggled,viewer,[viewer,hideRoof,roofDepth](bool on){viewer->setLeftRoofHidden(on);roofDepth->setEnabled(on&&hideRoof->isEnabled());});
    QObject::connect(roofDepth,QOverload<double>::of(&QDoubleSpinBox::valueChanged),viewer,[viewer](double depth){viewer->setLeftRoofDepth(float(depth));});
    auto applyDisplay=[&] {for(int pane=0;pane<2;++pane) viewer->setPaneDisplay(pane,paneMode[pane]->currentIndex(),paneColumns[pane]->isChecked(),paneGround[pane]->isChecked());};
    for(int pane=0;pane<2;++pane) {
        QObject::connect(paneMode[pane],QOverload<int>::of(&QComboBox::currentIndexChanged),&window,[&](int){applyDisplay();});
        QObject::connect(paneColumns[pane],&QCheckBox::toggled,&window,[&](bool){applyDisplay();});QObject::connect(paneGround[pane],&QCheckBox::toggled,&window,[&](bool){applyDisplay();});
    }
    auto applyColor=[&] {const bool scalar=colorSelect->currentIndex()!=0;minimum->setEnabled(scalar);maximum->setEnabled(scalar);autoRange->setEnabled(scalar);
        if(minimum->value()<=maximum->value()) viewer->setColorScale(colorSelect->currentIndex(),float(minimum->value()),float(maximum->value()));};
    auto resetRange=[&] {if(activeCloud) {const bool intensity=colorSelect->currentIndex()==1;const QSignalBlocker a(minimum),b(maximum);minimum->setValue(intensity?activeCloud->intensityP05:activeCloud->zP05);maximum->setValue(intensity?activeCloud->intensityP95:activeCloud->zP95);}applyColor();};
    QObject::connect(colorSelect,QOverload<int>::of(&QComboBox::currentIndexChanged),&window,[&](int){resetRange();});QObject::connect(autoRange,&QPushButton::clicked,&window,resetRange);
    for(auto *spin:{minimum,maximum}) QObject::connect(spin,QOverload<double>::of(&QDoubleSpinBox::valueChanged),&window,[&](double){applyColor();});
    QObject::connect(highlight,&QCheckBox::toggled,viewer,[viewer](bool on){viewer->setFillHighlight(on);});QObject::connect(axes,&QCheckBox::toggled,viewer,[viewer](bool on){viewer->setAxesVisible(on);});
    /** @brief 固定處理選項；匯出沿用同一程序且不受顯示選擇影響。
     *         Hold processing settings fixed; exports use the same pipeline independently of pane selections.
     */
    auto load=[&](const QString &destination=QString()) {
        if(watcher.isRunning()) return;
        select->setEnabled(false);exportMesh->setEnabled(false);exportPly->setEnabled(false);volumeSummary->clear();volumeTable->setRowCount(0);status->setText(destination.isEmpty()?"正在處理點雲與計算體積…":"正在匯出完整處理後點雲…");
        const QString path=dir.filePath(select->currentText());const bool calibrated=!parser.isSet("raw-coordinates");
        watcher.setFuture(QtConcurrent::run([path,limit,noise,fill,calibrated,structures,volume,destination]{return std::make_shared<Cloud>(processCloud(path,limit,noise,fill,calibrated,structures,volume,destination));}));
    };
    QObject::connect(&watcher,&QFutureWatcher<std::shared_ptr<Cloud>>::finished,&window,[&] {
        auto c=watcher.result();select->setEnabled(true);
        if(!c->error.isEmpty()) {activeCloud.reset();exportMesh->setEnabled(false);hideRoof->setEnabled(false);roofDepth->setEnabled(false);viewer->clearCloud();exportPly->setEnabled(false);volumeTable->setRowCount(0);volumeTable->hide();volumeSummary->clear();status->setText("處理失敗："+c->error);if(parser.isSet("smoke-test")) app->exit(2);return;}
        activeCloud=c;exportMesh->setEnabled(c->meshVolume.enabled&&c->meshVolume.error.isEmpty()&&!c->meshVolume.surface.empty());hideRoof->setEnabled(c->hasRoofPlanes);roofDepth->setEnabled(c->hasRoofPlanes&&hideRoof->isChecked());exportPly->setEnabled(true);viewer->setCloud(c,axis);resetRange();applyDisplay();
        if(!c->exportError.isEmpty()) status->setText("匯出失敗："+c->exportError);else status->setText(c->exportedPath.isEmpty()?QString():"已匯出："+c->exportedPath);
        const auto &v=c->volume;volumeTable->setRowCount(0);volumeTable->setVisible(v.enabled&&v.error.isEmpty());
        if(!v.enabled) volumeSummary->setText("體積計算未啟用。");
        else if(!v.error.isEmpty()) volumeSummary->setText("無法計算體積："+v.error);
        else {
            QString text=QString("柵格淨體積 %1 m³　｜　%2 個料堆　｜　柱體扣除 %3 m³\n10 × 10 cm 柵格・最低離地高度　｜　資料覆蓋 %4%　｜　場內缺測 %5 m²（未計入體積）")
                .arg(v.totalVolume,0,'f',1).arg(qulonglong(v.piles.size())).arg(v.columnExcludedVolume,0,'f',1).arg(v.footprintCells?100.*v.coveredCells/v.footprintCells:0,0,'f',1).arg(v.unknownCells*.01,0,'f',1);
            text+=QString("\n地面 z = %1x + %2y + %3　｜　%4").arg(v.ground[0],0,'f',6).arg(v.ground[1],0,'f',6).arg(v.ground[2],0,'f',3).arg(v.manualGround?"手動地面":QString("自動地面 RMSE %1 m").arg(v.groundRMSE,0,'f',3));
            if(c->surfaceFilter.comparisonValid) text+=QString("　｜　末段過濾體積差 %1 m³（同地面）").arg(c->surfaceFilter.volumeAfter-c->surfaceFilter.volumeBefore,0,'f',1);
            if(c->meshVolume.enabled) {
                const auto &m=c->meshVolume;
                if(!m.error.isEmpty()) text+="\nMesh 試算無法完成："+m.error;
                else text+=QString("\nMesh 試算 %1 m³　｜　與柵格差 %2 m³　｜　柱體扣除 %3 m³　｜　TIN 覆蓋 %4%　｜　%5 個三角形（支點 %6 m）")
                    .arg(m.totalVolume,0,'f',1).arg(m.totalVolume-v.totalVolume,0,'f',1).arg(m.columnExcludedVolume,0,'f',1).arg(m.footprintArea?100*m.coveredArea/m.footprintArea:0,0,'f',1).arg(qulonglong(m.surface.size())).arg(m.supportSpacing,0,'f',2);
            }
            volumeSummary->setText(text);
            for(const auto &p:v.piles) {const int row=volumeTable->rowCount();volumeTable->insertRow(row);const double values[]{double(p.id),p.volume,p.area,p.measuredVolume,p.estimatedVolume,p.enclosedUnknownCells*.01};for(int col=0;col<6;++col) volumeTable->setItem(row,col,new QTableWidgetItem(QString::number(values[col],'f',col?2:0)));}
        }
        if(parser.isSet("smoke-test")) QTimer::singleShot(1500,app.get(),[&] {
            bool ok=viewer->graphicsReady();for(int pane=0;pane<2&&ok;++pane) for(int mode=0;mode<6&&ok;++mode) {paneMode[pane]->setCurrentIndex(mode);ok=!viewer->grabFramebuffer().isNull();}
            std::cout<<"Compact pane selections: "<<(ok?"PASS":"FAIL")<<std::endl;app->exit(ok?0:2);
        });
    });
    QObject::connect(select,QOverload<int>::of(&QComboBox::currentIndexChanged),&window,[&](int){load();});
    QObject::connect(exportMesh,&QPushButton::clicked,&window,[&] {
        if(watcher.isRunning()||!activeCloud) return;
        const auto path=QFileDialog::getSaveFileName(&window,"匯出計算用開放料面 Mesh",dir.filePath(QFileInfo(select->currentText()).completeBaseName()+"_surface_mesh.ply"),"PLY mesh (*.ply)");if(path.isEmpty()) return;
        try {exportMeshSurface(path,activeCloud->meshVolume);status->setText("已匯出料面 Mesh："+path);}catch(const std::exception &e){status->setText("Mesh 匯出失敗："+QString::fromUtf8(e.what()));}
    });
    QObject::connect(exportPly,&QPushButton::clicked,&window,[&] {
        if(watcher.isRunning()||!activeCloud) return;
        const QString initial=dir.filePath(QFileInfo(select->currentText()).completeBaseName()+"_processed_filled.ply");
        const auto path=QFileDialog::getSaveFileName(&window,"匯出完整處理後點雲",initial,"PLY (*.ply)");if(path.isEmpty()) return;
        if(!QDir().mkpath(QFileInfo(path).absolutePath())) {status->setText("無法建立輸出目錄。");return;}load(path);
    });
    applyDisplay();applyColor();window.resize(1400,900);window.show();load();
    const int code=app->exec();watcher.waitForFinished();return code;
}
