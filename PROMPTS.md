# Session Log

## 2026-03-18

---

### Prompt 1

3711852.pdf の実装を行いたいです。C++ とOpenGL(freeglut/glut), glm を用いてFig.5 の再現を行いたいです。CMake を用い、Windows, macos, ubuntu で同一ソースからクロスプラットフォームな実装にしたいです。その際にWindows, ubuntu は、freeglut を用い、macos はApple 標準のGLUT を用いてください。まずは、ローカルにある論文を精読し実装計画の .md を作ってください。また私のプロンプトはすべて .md ファイルに記録してください。

### Response 1 (summary)

- PyMuPDF (fitz) を pip install し、PDF 全文を抽出して精読
- 論文: "Implicit Bonded Discrete Element Method with Manifold Optimization"
  ACM Trans. Graph. 44(1), Article 9 (January 2025), DOI: 10.1145/3711852
- Fig. 5 = 三点曲げ試験を 3 スケール (Low/Middle/High) で行い、BDEM のスケール整合性を示す図
- `IMPLEMENTATION_PLAN.md` を作成
  - 数学的基盤 (質量行列・ポテンシャル・破断判定・多様体最適化) を全方程式込みで記載
  - 6 フェーズの実装計画
  - CMake クロスプラットフォーム設計 (macOS=Apple GLUT Framework, Win/Ubuntu=freeglut)
  - 依存ライブラリのインストール手順
- `PROMPTS.md` を作成 (本ファイル)

---

### Prompt 2

Phase1 を実行してください。c++ のソースコードは/src 以下にまとめて CMake ファイルはルートにおいてください。またソースコード内のコメントに関してはすべて英語で記載して、インデントはスペース２で実装してください。

### Response 2 (summary)

以下ファイルを新規作成しビルド・動作確認まで完了。

#### 作成ファイル

| ファイル | 内容 |
|---------|------|
| `CMakeLists.txt` | クロスプラットフォームビルド設定。macOS は `find_library` で Apple GLUT/OpenGL Framework、Windows/Ubuntu は `find_package(GLUT)` で freeglut を検出。GLM は `find_package(glm CONFIG)` → ヘッダ検索フォールバック |
| `CMakePresets.json` | `windows-vcpkg` / `unix-default` プリセット定義 |
| `src/Particle.h` | 粒子データ構造 (`pos`, `rot`, `vel`, `angVel`, `mass`, `radius`, `fixed`, `isSupport`, `isLoad`) |
| `src/Bond.h` | ボンドデータ構造 (`i`, `j`, `l0`, `d0`, `q0`, `broken`) |
| `src/HexPacking.h` | `BeamConfig` 構造体、`generateBeam` / `generateBonds` / `tagBoundaryParticles` 宣言 |
| `src/HexPacking.cpp` | Y-Z 断面に 2D 三角格子 (HCP近似) を生成し X 方向に繰り返す粒子生成。空間ハッシュグリッド O(n) でのボンド生成。三点曲げ用境界条件タグ付け |
| `src/main.cpp` | GLUT ウィンドウ (1200×450px)、3 分割ビューポート、GL_POINTS + GL_LINES による静的描画。灰=通常粒子、青=サポート、赤=荷重点 |

#### ビルド手順 (Windows)

```
cd ibdem
cmake --preset windows-vcpkg
cmake --build build --config Release
copy C:/vcpkg/installed/x64-windows/bin/freeglut.dll build/Release/
./build/Release/ibdem.exe
```

#### 実行結果

```
Generating scale 0 (Low)    r=0.0180 ...  particles: 504    bonds: 1,522
Generating scale 1 (Middle) r=0.0110 ...  particles: 2,024  bonds: 6,856
Generating scale 2 (High)   r=0.0050 ...  particles: 22,725 bonds: 84,615
Initialisation complete.
Press ESC or Q to quit.
```

ウィンドウが開き 3 スケールの梁が並んで表示されることを確認。

---

### Prompt 3

これまでの出力と私のプロンプトを .md に記録して終了してください。

### Response 3 (summary)

本ファイル (`PROMPTS.md`) を全セッションのログとして更新して終了。

---

## 2026-03-19 (Session 2 — Phase 2: Physics Solver)

---

### Prompt 4

残りの作業をすべて行ってください。私はこれから寝るので問い合わせることなく実装を行い実行して画像を確認し検証も行ってください。今のところ表示文字列が重なって表示されているのでそのようなところも一通り修正してください。残りの作業をタイムアウトするまで勧めてください。

### Response 4 (summary) — Phase 2 物理ソルバー実装

#### 新規・変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `src/BondForce.h` | `BondGradient`, `bondEnergy`, `bondGradient`, `bondStress` の宣言 |
| `src/BondForce.cpp` | 伸び・せん断・曲げ-ねじりエネルギー、勾配、応力の実装 |
| `src/Simulation.h` | `SimConfig`, `Simulation` クラス宣言 |
| `src/Simulation.cpp` | 陰的時間積分 (勾配降下 + Armijo ライン探索 + S³多様体更新) |
| `src/main.cpp` | シミュレーションループ、荷重制御、GUI描画、`-headless N` テストモード追加 |

#### Phase 2 で実装した主要アルゴリズム

1. **ボンドエネルギー** (BondForce.cpp)
   - 伸びエネルギー: `V_stretch = 0.5*kn*(dist-l0)^2`
   - せん断エネルギー: `V_shear = 0.5*kt*theta^2`, `theta = angle(d_cur, R(qc)*d0)`
     - `d_cur = normalize(pj-pi)` で位置-姿勢結合を実現
   - 曲げ-ねじりエネルギー: `V_bt = 0.5*t^T*K*t`, `t = R0*Gl(qj)*qi`
   - 各エネルギーの位置勾配・四元数勾配を解析的に計算

2. **対角前処理付き勾配ステップ** (Simulation.cpp)
   - `h_p[i] = m/dt^2 + Σ(kn+kt)`
   - `h_q[i] = Iq/dt^2 + Σ(EI/l0)`
   - Armijo バックトラッキング・ライン探索
   - S³上の測地線更新: `q_new = q*cos(|dq|) + (dq/|dq|)*sin(|dq|)`

3. **ヘッドレステスト** (`-headless N`): ウィンドウなしで N フレーム実行し破断フレームを出力

#### パラメータ調整の試行錯誤

| 試行 | dt | loadVel | 結果 | 問題 |
|------|-----|---------|------|------|
| 1 | 1e-2 | 0.15 | frame 0 破断 | loadVel 過大 (1フレームで3.3mm超) |
| 2 | 1e-2 | 0.015 | frame 0-1 破断 | dt 大 → Jacobi 発散, 姿勢結合なし → スケール非整合 |
| 3 | 論文値(2e-5等) | 7.5/15.0/2.5 | frame 0 破断 | τ が数千倍以上 → 別の問題あり |

#### 現在判明している問題 (未解決)

- 全スケールで frame 0 に破断
- 診断出力結果:
  ```
  Low    maxTau=6.5e+06 (tauC=3e4) — 約217倍
  Middle maxTau=2.4e+07            — 約787倍
  High   maxTau=6.3e+11            — 天文学的値 (数値爆発の疑い)
  ```
- High スケールは bond 84,553/84,615 本が破断 → ほぼ全滅
- せん断エネルギーの位置-姿勢結合 (`d_cur = normalize(pj-pi)`) は実装済み
- 四元数が単位球から外れていないはずだが `t = R0*Gl(qj)*qi` が異常に大きい

---

### Prompt 5

ここまでの作業を .md ファイルに記録してください

### Response 5 (summary)

本ファイル (`PROMPTS.md`) を現セッション分まで更新。

---

## 2026-03-19 (Session 3 — PCG ソルバーのデバッグ・スケール整合性調査)

---

### Prompt 6

引き続きタイムアウトになるまで作業を進めてみてください

### Response 6 (summary) — PCG 修正・ブロック座標降下の実装

前セッションで発見した「High スケールの応力が極端に小さい (frame 0: 58.88 Pa, 期待値 ~1200 Pa)」問題の調査と修正。

#### 発見した主要バグ

**Bug 1: PCG 収束判定の絶対値/相対値誤り**

```cpp
// 旧 (誤): 閾値が n_free に比例 → High スケールで即座に収束と誤判定
if (n_free > 0 && gnorm < epsilon * n_free) break;

// 新 (正): 相対収束 → 初期残差比で判定
if (rz0 > 0.0f && rz < epsilon * epsilon * rz0) break;
```

High スケールは粒子数 22725 のため閾値が 22.7 倍緩くなり、PCG が 0 回反復で終了していた。
修正後: PCG が 20〜26 回反復で収束 (rz_ratio ≈ 7e-5) していることを確認。

**Bug 2: 回転 Jacobi ループが位置を汚染**

```cpp
// 旧 (誤): gradientStep が位置・回転を両方更新
for (int iter = 0; ...) {
    gradientStep(1.0f);   // 位置も動かしてしまう
    ...
}

// 新 (正): PCG 最適位置を保存・復元し、回転のみ更新
std::vector<glm::vec3> pos_save(N);
for (int i = 0; i < N; i++) pos_save[i] = particles[i].pos;
for (int iter = 0; ...) {
    gradientStep(1.0f);
    for (int i = 0; i < N; i++) particles[i].pos = pos_save[i]; // 位置復元
    ...
}
```

**Bug 3: PCG+Jacobi のブロック分離が不十分**

PCG (位置最適化) と Jacobi (回転最適化) を交互に繰り返すブロック座標降下が必要。
1 パスだけでは、PCG が固定した回転を前提に位置を解いても、回転が更新されると位置の最適性が崩れる。

→ **5 回の外部ループ**で PCG+Jacobi を交互に実行するよう変更。

```cpp
const int outerLoops = 5;
const int innerPCG   = maxIter / outerLoops;  // 200/5 = 40 回 PCG
const int innerJac   = 20;                    // 回転 Jacobi は収束が速い

for (int outer = 0; outer < outerLoops; ++outer) {
    solvePCG(is_load, load_target, innerPCG, epsilon);
    // 境界条件クランプ
    for (int iter = 0; iter < innerJac; ++iter) {
        // 位置保存 → gradientStep → 位置復元 (回転のみ更新)
    }
}
```

#### 各修正後のヘッドレステスト結果

| 修正内容 | Low | Middle | High | 目標 |
|---------|-----|--------|------|------|
| PCG 相対収束 + maxIter=200 | 20 | 9 | 18 | 23, 25, 25 |
| + 位置復元バグ修正 | 23 ✓ | 7 | 17 | 23, 25, 25 |
| + 5 外部ループ | 23 ✓ | 7 | 17 | 23, 25, 25 |

Low スケールが目標 23 フレームに到達。Middle と High は依然として不整合。

#### Frame 0 応力の比較

| スケール | frame 0 maxSigma | frame 0 maxTau | 期待値 | 状況 |
|---------|-----------------|----------------|--------|------|
| Low (r=0.018) | 2844 Pa | 2390 Pa | ~1875 Pa | 1.5× 過大 |
| Middle (r=0.011) | 1845 Pa | 346 Pa | ~1875 Pa | ≈ 正常 |
| High (r=0.005) | 141 Pa | 55 Pa | ~1875 Pa | 13× 過小 |

#### 診断: PCG は実際に収束している

PCG の反復回数と収束率の診断出力を追加:
```
[sc0 frame0 PCG] iters=22 rz_ratio=6.756e-05   ← Low
[sc1 frame0 PCG] iters=23 rz_ratio=7.159e-05   ← Middle
[sc2 frame0 PCG] iters=24 rz_ratio=7.096e-05   ← High
```
PCG は全スケールで 20〜26 回反復で収束 (rz_ratio ≈ 7e-5) している。位置は収束しているはずだが、応力は 13× の差異がある。

#### 現在判明している問題 (未解決)

**問題 1: isLoad 粒子数のスケール不整合**

`tagBoundaryParticles` で `searchRadius = 3*r` を使用。半径が r に比例するため、実際に選ばれる isLoad 粒子数が整合しない:

| スケール | isLoad 粒子数 (推定) | N × kn 比 |
|---------|-------------------|-----------|
| Low | ~9 | 0.162 |
| Middle | ~6 | 0.066 |
| High | ~6 | 0.030 |

N × kn (有効荷重) が 5.4× 差異。これが応力差の主因の可能性。

**問題 2: 応力が非線形に増大 (Middle)**

Middle スケールの frame ごとの応力増分:
```
frame 0: 1845 Pa
frame 1: 3795 Pa (+1950/frame)
frame 2: 5773 Pa (+1978/frame)
frame 3: 7965 Pa (+2192/frame)
frame 4: 13670 Pa (+5705/frame ← 突然3倍に増大)
frame 5: 19870 Pa (+6200/frame)
frame 6: 26150 Pa (+6280/frame)
frame 7: 31360 Pa ← 破断
```
frame 3-4 で stress の伸び率が突然 3× に増大するのは振動 (動的増幅) の可能性。

**問題 3: High スケールの応力が構造的に過小**

PCG が収束しているにもかかわらず応力が 13× 小さい原因:
- shear energy の勾配が回転 (q_hat ≈ identity) 前提で計算されているため、位置 PCG が shear-dominated (せん断支配) 解に収束している可能性
- 真の beam bending 解には Euler-Bernoulli 曲げ剛性が必要だが、回転が identity の間は shear stiffness のみが有効
- 5 外部ループ + 20 Jacobi では回転の収束が不十分な可能性 (大域モード: L/r ≈ 200 で条件数が悪化)

#### 現在のパラメータ (main.cpp)

```cpp
static const SimConfig SIM_CFG[NUM_SCALES] = {
  //  dt          E       nu    tauC   loadVel  grav  maxIter  eps
  { 5.35e-4f, 1e7f, 0.3f, 3e4f, 0.2491f, 0.0f, 200, 1e-4f }, // Low
  { 3.27e-4f, 1e7f, 0.3f, 3e4f, 0.4075f, 0.0f, 200, 1e-4f }, // Middle
  { 1.49e-4f, 1e7f, 0.3f, 3e4f, 0.8946f, 0.0f, 200, 1e-4f }, // High
};
```

#### 次に試すべきこと

1. **isLoad 粒子数の均一化**: `searchRadius` を絶対値 or `sqrt(r)` スケールに変更
2. **より多い外部ループ**: outerLoops を 5 → 20 に増やして収束確認
3. **ダンピング追加**: 速度に人工減衰係数 α を掛けて振動を抑制
4. **荷重制御切り替え**: 変位制御 → 力制御 (prescribed force) に変更
5. **bend-twist 項を PCG Hessian に追加**: 現在は stretch + shear(pos) のみ

---

### Prompt 7

ここまでの作業を .md ファイルに書き出してください

### Response 7 (summary)

本ファイル (`PROMPTS.md`) を Session 3 の分まで更新。

