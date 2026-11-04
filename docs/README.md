# 現行版本說明

本目錄只保留目前程式與生效模型的說明。場域及結構設定以專案`config/` JSON 為準，文件不另存模型快照。

| 文件 | 內容 |
|---|---|
| [處理流程](processing_pipeline.md) | 單份 PLY 至淨體積、地面與倒凸邊界 |
| [介面操作](compact_interface.md) | 左右視圖、著色、隱藏屋頂、地面常數調整與匯出 |
| [結構模型](columns_final.md) | 生效檔、目前分群、ID 與檢查方式 |
| [梁柱尺寸](structural_beams.md) | S1–S9、S5_BOX 的現行尺寸與檢查圖 |
| [結構聯集扣除](structural_union_volume.md) | 只扣埋入料堆部分，重疊只扣一次 |
| [Mesh 體積](mesh_volume_trial.md) | TIN 積分、參數及開放上表面匯出限制 |
| [體積前離群清理](pre_volume_filter.md) | 局部穩健表面濾波與門檻 |
| [固定座標校準](calibration_extension_walls.md) | 未來掃描註冊、延伸端牆與失敗排查 |

[舊審查 ID 對照](column_id_mapping.csv) 僅用於追溯舊柱體編號，不是現行實體清單。現行清單請查 `../config/column_semantics.json` 與 `../config/column_catalog.json`。

歷史備份、舊版圖表、診斷 PLY、試算 JSON 及已被取代的方法紀錄已移除。幾何計算測試通過不代表已完成現場庫存精度驗證。
