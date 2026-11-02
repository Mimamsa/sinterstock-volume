# 固定建築座標校準與失敗排查

`config/site_model.json` 是永久參考。未來每份掃描獨立找固定屋頂／牆面並註冊，不用變動料堆對齊，也不重新建立原點。屋頂已知斜率 0.2，Z 朝上、X 朝北、Y 朝西。

## 使用方式

新掃描開啟時，`pointcloud_compare` 可在背景呼叫同一建置目錄的 `calibrate_clouds`，通過後將矩陣寫入來源目錄的 `calibration.json`；原始 PLY 不改。側錄綁定場域 ID 與來源 SHA-256，檔案改變後不得重用舊矩陣。

```bash
# 只檢查，不寫入
./build/calibrate_clouds /path/to/scan.ply --site "$PWD/config/site_model.json" --analyze-only
# 只保存矩陣，不輸出整份校準 PLY
./build/calibrate_clouds /path/to/scan.ply --site "$PWD/config/site_model.json" --matrix-only
```

預定義模型集中在專案 `config/`；查找順序為有效側錄路徑、來源目錄的 `config/site_model.json`／舊版 `site_model.json`，再找工作目錄及執行檔上一層的 `config/`。舊側錄路徑失效時會繼續查找，不需搬動掃描檔。

移動執行檔時，要一起攜帶 `calibrate_clouds` 與 `config/`，或以命令列明確指定模型。已校準匯出的 PLY 有識別標記，避免再次套用轉換。

## 延伸料區端牆

若最大牆面密度峰把延伸區端牆當成主建築端牆，量到的跨度可能不符永久模型。程式在主峰跨度不符時會排除屋面樣本、找多個外側牆峰，挑選符合固定跨度且兩側都有支持的平面對，再執行剛體精修與保留樣本驗證。

既有尺寸允差不放寬：寬 0.7 m、長 1 m；結構匹配 RMSE ≤ 0.18 m、P95 ≤ 0.28 m、匹配率 ≥ 25%。180° 方向不明時仍拒絕，可在有方向證據時使用 `--orientation forward/reverse`。

主牆未掃到、固定結構不足、方向不明或品質不足，應保留拒絕結果，不能只依給定尺寸補造證據。`Automatic calibration rejected scan` 後面的詳細訊息可用上述唯讀命令重現。`--ground-plane` 在校準後才用於地面與體積，不會修正配準失敗。

測試涵蓋不同姿態、180°、延伸端牆、變動料堆、來源 SHA 失效、固定模型不變與結構不足拒絕。平面殘差是關聯範圍內的配準指標，不是測量級絕對精度保證。
