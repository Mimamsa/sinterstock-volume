# SinterStock Volume

專為**三四階燒結礦室內料場**開發的點雲料堆體積計算程式。採用 C++／Qt5，整合座標校準、雜點與建築結構移除、料面補洞，以及柵格與 Mesh 體積計算，提供處理前後點雲視覺化與 PLY 匯出功能。

GitHub 儲存庫名稱：`sinterstock-volume`。

完整說明入口：[docs/README.md](docs/README.md)。本專案只保留現行版本文件，生效模型在`config/` JSON，不使用文件內的舊快照。

## 建置與執行

需求：CMake ≥ 3.16、C++17、Qt5 ≥ 5.12、Eigen3、Qhull、Boost，以及桌面 OpenGL 3.3。

Ubuntu／Debian 套件：`build-essential cmake qtbase5-dev libgl1-mesa-dev libqhull-dev libeigen3-dev libboost-dev`。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
./build/pointcloud_compare /path/to/ply_directory
```

讀取指定目錄的 PLY，不搜尋子目錄；每份獨立處理，不融合可能變動的料堆。座標與長度為公尺。

```bash
# 水平手動地面：共同建築座標 Z = -24.31 m
./build/pointcloud_compare /path/to/scans --ground-plane 0,0,-24.31
# 保留自動擬合的 a、b，指定原點到地面的法線距離 d（不可與 --ground-plane 併用）
./build/pointcloud_compare /path/to/scans --ground-distance -24.31
# 無桌面：逐檔處理並輸出報告
./build/pointcloud_compare /path/to/scans --check
# 無桌面查看完整參數
./build/pointcloud_compare --check --help
# 完整處理後點雲匯出，不受顯示點數上限影響
./build/pointcloud_compare /path/to/scans --check --output-dir /path/to/new_exports
```

原始 PLY 及既有輸出不覆寫。原始量測不做前置三維體素質心化；顯示抽樣不影響計算與完整匯出。

## 操作與模型

左右視圖共用相機，左鍵旋轉、右鍵或 Shift 拖曳平移、滾輪縮放、雙擊重新取景。每側可選原始／處理後／實測／補點／結構證據，並開關結構框與地面；兩側預設都顯示共同地平面網格與料堆邊界。左側可隱藏屋頂及構架；著色支援強度或 Z、自動 P5–P95。可切換「拍照模式」隱藏角落文字與色階，並按「擷取畫面」將左右視窗存成單張 PNG。見[介面操作](docs/compact_interface.md)。

| 生效檔 | 用途 |
|---|---|
| `config/site_model.json` | 永久共同建築座標與固定結構參考 |
| 來源目錄 `calibration.json` | 每份掃描的註冊矩陣及 SHA-256 |
| `config/processing_region.json` | 倒凸字型凹多邊形場域 |
| `config/column_semantics.json` | 已確認柱、梁、箱體截面及上下界 |
| `config/column_catalog.json` | 群組與成員清單 |
| `config/column_groups.json`、`config/column_rows.json`、`config/structural_beams.json` | 四柱組、柱列及橫樑尺寸關聯 |

目前結構模型共 **203 個實體、14 群**，定案後不自動恢復候選。S1–S8 柱高 14 m；S5 最低層為高 4.5 m 的 S5_BOX，無最低梁，第二層梁底 9.55 m；S9 柱高 11.5 m、三層梁。完整尺寸與 ID 見[結構模型](docs/columns_final.md)及[梁柱尺寸](docs/structural_beams.md)。

## 體積原則

清理、結構排除、補洞及末段離群清理後，僅由實測擬合 `Z=aX+bY+c`。可透過 `--ground-distance d` 或 GUI「調整原點到地面距離 d」保留擬合 a、b，指定原點到地面的最短距離 d（帶正負號，公尺），換算 `c = d√(1+a²+b²)`，d 為負時地面位於原點下方，為正時位於上方；GUI 需按「套用地面並重新計算」才生效，取消勾選後套用可恢復原擬合。完整手動平面仍可用 `--ground-plane a,b,c`，兩種指定方式不可併用。全流程見[處理流程](docs/processing_pipeline.md#8-估計地平面或採用手動地面)。

- 柵格固定 10 cm，每格使用最低離地高度，實測優先、無實測才用補點，未知格不計。
- 只在倒凸有效場域內積分，保留兩端延伸料地與預設 0.4 m 牆邊內縮。
- 只扣[已確認結構聯集與料體的交集](docs/structural_union_volume.md)，露出部分不扣、梁下空間保留、重疊只扣一次。
- Mesh 使用獨立 TIN 積分，共用資料／地面／結構。匯出是開放料面，並非可直接求相同淨體積的封閉實體。見[Mesh 方法與選項](docs/mesh_volume_trial.md)。

盤點參考 47,000–50,000 m³ 與幾何估算的差異尚未解決；測試通過不等於現場庫存精度已驗證。

## 常用調整

```bash
./build/pointcloud_compare /path/to/scans --surface-outlier-tolerance 0.8
./build/pointcloud_compare /path/to/scans --no-mesh-volume
./build/pointcloud_compare /path/to/scans --mesh-support-spacing 0.5 --mesh-max-edge 2
```

[離群清理](docs/pre_volume_filter.md)說明門檻與限制。[固定座標校準](docs/calibration_extension_walls.md)說明未來掃描、延伸端牆及 `Automatic calibration rejected scan` 排查；地面參數不修正配準失敗。
