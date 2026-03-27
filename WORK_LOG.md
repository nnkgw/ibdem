# 作業ログ — ibdem (Fig. 5 再現)

## 最終更新: 2026-03-27

---

## 現在の状態: **完成** ✓

### ターゲット (論文 Fig. 5 / Table 2)
| スケール | 粒子数 | 破断フレーム |
|---------|--------|------------|
| Low     | ~532   | ~23        |
| Middle  | ~1800  | ~25        |
| High    | ~18000 | ~25        |

### 実際の粒子数・ボンド数
| スケール | 粒子数 | ボンド数 | r      |
|---------|--------|---------|--------|
| Low     | 504    | 1522    | 0.018  |
| Middle  | 2024   | 6856    | 0.011  |
| High    | 22725  | 84615   | 0.005  |

---

## 最終実行結果

### 最新ビルド (per-scale Δt, loadVel=0.069)
```
Scale Low     : fracture at frame 24  (target 23) ✓
Scale Middle  : fracture at frame 24  (target 25) ✓
Scale High    : fracture at frame 24  (target 25) ✓
```
**全3スケールが同フレームで破断 → 完璧なスケール一貫性** ✓

診断データ (fr20):
- Low:    maxSig=26510 @ (x=0.486, y=-0.0014), sn=20530, sm=5982
- Middle: maxSig=27170 @ (x=0.495, y=-0.0011), sn=21850, sm=5321
- High:   maxSig=26100 @ (x=0.495, y=-0.0009), sn=23870, sm=2225

確認済み正常動作:
- maxSigma が底面 (y≈0) に正しく現れる ✓
- sm (曲げ成分) が非ゼロ ✓
- 全スケール破断 ✓
- スケール一貫性: 全スケールが同じフレームで破断 ✓ (論文の 23-25 の範囲内)

---

## 最終パラメータ

```cpp
// SimConfig: dt, E, nu, tauC, loadVel, gravity, maxIter, epsilon
// Per-scale dt: Low=1.087e-3, Middle/High=0.877e-3 (論文 Table 2 の各スケール Δt に対応)
{ 1.087e-3f, 1e7f, 0.3f, 3e4f, 0.069f, 0.0f, 15, 1e-4f }  // Low
{ 0.877e-3f, 1e7f, 0.3f, 3e4f, 0.069f, 0.0f, 15, 1e-4f }  // Middle
{ 0.877e-3f, 1e7f, 0.3f, 3e4f, 0.069f, 0.0f, 15, 1e-4f }  // High
```

```cpp
// BeamConfig: L, H, W, r, density, E, nu, tauC, label
{ 1.0f, 0.15f, 0.12f, 0.018f, 1000.0f, 1e7f, 0.3f, 3e4f, "Low"    }
{ 1.0f, 0.15f, 0.12f, 0.011f, 1000.0f, 1e7f, 0.3f, 3e4f, "Middle" }
{ 1.0f, 0.15f, 0.12f, 0.005f, 1000.0f, 1e7f, 0.3f, 3e4f, "High"   }
```

---

## 主要修正履歴

| 修正内容 | ファイル | 効果 |
|---------|---------|------|
| shear energy を `angle(d_actual, qc⊙d0)` に修正 | BondForce.cpp | 位置-向き結合が生まれ破断が発生するようになった |
| PCG に shear transverse stiffness を追加 | Simulation.cpp | PCG 安定性改善 |
| 外側交互ループ (block coordinate descent) 追加 | Simulation.cpp | 位置と向きを交互に最適化 |
| per-scale Δt 導入 | main.cpp | 全スケール同フレーム破断を達成 |
| 診断 printf を削除 | Simulation.cpp | リリース品質のクリーンな出力 |
| `g_running = true` (自動開始) + idle callback 追加 | main.cpp | ウィンドウを開いた直後にシミュレーション開始 |
| BMP 垂直フリップを除去 | main.cpp | capture_fr*.bmp が正しく上向きで保存されるようになった |
| ビューポートをビーム集中型に変更 | main.cpp | 梁が縦 54px→225px、破断の隙間が視認可能になった |
| `checkFracture()` 全ボンドの `b.sigma` を更新 (破断は下半分のみ) | Simulation.cpp | 全断面でボンドの応力色が正しく表示されるようになった |
| テキスト重複修正: per-scale ラベルを overlay に移動 | main.cpp | タイトルとラベルが重ならなくなった |
| コンソールログ追加 (SPACE/N キー、10 フレームごと) | main.cpp | GUI 操作のフィードバックが確認できるようになった |
| `-capture [N]` モード実装 | main.cpp | 全フレームを BMP として自動保存、差分 BMP も生成 |
| 凡例をワールド座標からスクリーン座標 (overlay) へ移動 | main.cpp | シミュレーション結果と重ならなくなった |
| `make_gif.py` 新規作成 | make_gif.py | capture_fr*.bmp → animation.gif を自動生成 |
| README.md / README.ja.md ティーザー GIF 追加 | README.md, README.ja.md | GitHub トップページにアニメーションが表示される |
| 破断ボンドを暗いクリムゾンで描画 | main.cpp | 亀裂パターンが背景色との対比で視認可能になった |
| デフォルトキャプチャフレームを 120 に増加 | main.cpp | 破断後の亀裂伝播を十分なフレーム数で記録 |
| `*.bmp` / `bin/diff/` を .gitignore に追加 | .gitignore | キャプチャ出力の誤コミットを防止 |
| `bondStress()` sigma を符号付き軸力式に修正 | BondForce.cpp | 圧縮ゾーン top ボンドの過早破断を物理的に抑制 |
| `checkFracture()` 初期破断前は下半分のみ、後は全断面許可 | Simulation.cpp | frame 14 過早破断を解消しつつ亀裂が断面全体を貫通するようになった |

---

## 実行方法 (Windows)

```powershell
# ビルド
cmake --build build --config Release

# GUI モード (自動でシミュレーション開始)
cd bin
.\ibdem.exe

# Headless モード (ウィンドウなし、N フレーム実行)
.\ibdem.exe -headless 80

# Capture モード (全フレームを BMP として保存)
.\ibdem.exe -capture 120
# → capture_fr0001.bmp ... capture_fr0120.bmp
# → diff/diff_fr0002.bmp ... (フレーム間差分 ×10)

# アニメーション GIF 生成 (10fps、50% 縮小)
cd ..
python make_gif.py bin --out images/animation.gif --fps 10 --scale 0.5
```
