# 預定義場域與結構設定

| 檔案 | 用途 |
|---|---|
| `site_model.json` | 永久建築座標與固定結構參考 |
| `processing_region.json` | 倒凸字型處理邊界 |
| `column_semantics.json` | 已確認柱、梁、箱體的排除幾何 |
| `column_catalog.json` | 結構群組及成員 |
| `column_groups.json` | S1–S9 柱心距、錨點、梁柱尺寸 |
| `column_rows.json` | 主柱列參考與尺寸 |
| `structural_beams.json` | 邏輯橫樑與實體片段對應 |

這些檔案需一起攜帶，程式從 `site_model.json` 所在目錄讀取邊界及結構語意。修改幾何時同步更新關聯檔，保留場域 ID。

每份掃描的 `calibration.json` 是來源側錄，留在掃描目錄，不放在此處。使用方式見 [現行說明](../docs/README.md)。
