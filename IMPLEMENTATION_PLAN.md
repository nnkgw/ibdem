# 実装計画: Implicit BDEM with Manifold Optimization — Fig. 5 再現

## 論文情報
- **タイトル**: Implicit Bonded Discrete Element Method with Manifold Optimization
- **著者**: Jia-Ming Lu, Geng-Chen Cao, Chenfeng Li, Shi-Min Hu
- **掲載**: ACM Trans. Graph. 44, 1, Article 9 (January 2025)
- **DOI**: https://doi.org/10.1145/3711852

---

## 再現対象: Fig. 5

**内容**: 三点曲げ試験（Three-Point Bending）における3手法（FEM, MPM, BDEM）のスケール整合性比較。
本実装では **BDEMパート（右列）** のみを再現する。

- **3スケール**（Low/Middle/High）で同一荷重をかけ、ほぼ同一フレームで破断することを示す
- 各スケールで粒子を hexagonal close packing で充填したビームを生成
- 下部2点でサポート、上部中央1点で荷重（三点曲げ）
- 破断フレームがスケール間で一致することを視覚的に確認

### Table 2 より（BDEMパート）

| Scale  | # Particles | E       | ν   | τ_c   | Fracture Frame |
|--------|-------------|---------|-----|-------|----------------|
| Low    | 532         | 1E7 Pa  | 0.3 | 3E4 Pa| 23             |
| Middle | 1,800       | 1E7 Pa  | 0.3 | 3E4 Pa| 25             |
| High   | 18,000      | 1E7 Pa  | 0.3 | 3E4 Pa| 25             |

---

## 数学的基盤

### 粒子状態
各粒子 i の状態: `x_i = {p_i, q_i}` （位置 R³ + 回転クォータニオン S³）

### 拡張質量行列（Eq. 1）
```
Ms = diag(mI₃,  (8/5)mr²I₄)
```
- `m`: 粒子質量
- `r`: 球半径
- I₃, I₄: 3×3, 4×4 単位行列

### 最適化目的関数（Eq. 6、一次精度版）
```
Φ(x_{n+1}) = (1/2Δt²)(x_{n+1} - x̂)ᵀ Mₛ (x_{n+1} - x̂) + V(x_{n+1})
x̂ = xₙ + vₙΔt
```

### ポテンシャルエネルギー

#### 伸縮エネルギー（Eq. 8）
```
V_stretch = (1/2) kₙ (‖pᵢ - pⱼ‖ - l₀)²
kₙ = ES / l₀,  S = π((rᵢ+rⱼ)/2)²
```

#### せん断エネルギー（Eq. 9）
```
V_shear = (1/2) kₜ θ²
q_c = (qᵢ + qⱼ) / ‖qᵢ + qⱼ‖
d = q_c ⊙ d₀
θ = angle(d₀, d)
kₜ = GS / l₀
```

#### 曲げ・ねじりエネルギー（Eq. 10, 11）
```
V_bendtwist = (1/2) tᵀ K t
t = R₀ Gₗ(qⱼ) qᵢ
R₀ = Gₗ(q₀) Gᵣ(q₀)
K = diag(GJ/l₀, EI/l₀, EI/l₀)   ← ねじり, 曲げx, 曲げy
I = πr₀⁴/4,  J = πr₀⁴/2
```

#### ずり剛性・ヤング率の関係
```
G = E / (2(1+ν))
```

### 破断判定（Eq. 12, 13）
```
σ = ‖Fₙ‖/S + ‖M_b‖r₀/I   (引張応力)
τ = ‖Fₛ‖/S + ‖M_t‖r₀/J   (せん断応力)
破断条件: σ > σ_c  または  τ > τ_c
```
τ_c = 3E4 (Table 2), σ_c も同値とする（論文に明記なし）

内力の計算（Eq. 14–17）：
```
Fₙ = ∇ₚ V_stretch
Fₛ = ∇ₚ V_shear
‖M_b‖ = ‖⟨Kt, e₁⟩e₁ + ⟨Kt, e₂⟩e₂‖
‖M_t‖ = ‖⟨Kt, e₀⟩e₀‖
```

### クォータニオン操作（Appendix A）

#### 左・右乗算行列（Eq. 41, 42）
```
Qₗ(q) = [ sI₃ + t×   t ]    (4×4)
          [ -tᵀ        s ]

Qᵣ(q) = [ sI₃ - t×   t ]    (4×4)
          [ -tᵀ        s ]
```

#### ヌル空間演算子（Eq. 47）
```
Gₗ(q) = [ sI₃ - t× | -t ]   (3×4)
Gᵣ(q) = [ sI₃ + t× | -t ]   (3×4)
```

### 多様体上の勾配（Eq. 18, 19）
```
grad_p f = grad_p f̄
grad_q f = (I - qqᵀ) grad_q f̄
```

### 多様体上のヘッセ行列（Eq. 23–27）
```
Hess_{pᵢ,pⱼ} f = Hess_{pᵢ,pⱼ} f̄
Hess_{pᵢ,qⱼ} f = Hess_{pᵢ,qⱼ} f̄ (I - qⱼqⱼᵀ)
Hess_{qᵢ,pⱼ} f = (I - qᵢqᵢᵀ) Hess_{qᵢ,pⱼ} f̄
Hess_{qᵢ,qⱼ} f = (I - qᵢqᵢᵀ) Hess_{qᵢ,qⱼ} f̄ (I - qⱼqⱼᵀ)   i≠j
Hess_{qᵢ,qᵢ} f += −(qᵢᵀ grad_{qᵢ} f̄)(I - qᵢqᵢᵀ)
```

### ヌル空間縮約（Eq. 29）
```
G(q) H G(q)ᵀ y = z
y = G(q)Δx,   z = G(q) b
Δx = G(q)ᵀ y
```

---

## アルゴリズム（Algorithm 3: Stable Integrator）

```
入力: Δt, 収束閾値 ε
x₀ ← xₙ + vₙΔt
k ← 0
while true:
    b = grad f(xₖ) on manifold M
    if ‖b‖ < ε: break
    H = Hess f(xₖ) on manifold M
    H' = SPD projection of H
    εᵣ = min(0.5, √|grad f(xₖ)|)
    Solve H' Δx = -b in tangent space (ヌル空間縮約 + PCG)
    α = line search on manifold
    xₖ₊₁ ← xₖ + α Δx (geodesic update on S³)
    xₖ₊₁ の q を正規化
    k ← k + 1
xₙ₊₁ = xₖ₊₁
vₙ₊₁ = (xₙ₊₁ - xₙ) / Δt
衝突後速度補正（摩擦コーン投影）
ボンド破断判定
```

---

## タイムステップ設定（Section 5.1）
```
T = 2 / √λ_max
λ_max = (k/m) * N,   N ≈ 16 (hexagonal packing)
k ≈ (π/2) E r
Δt_explicit = 0.2 T
Δt_implicit = 20〜100 × T  (目安; Table 2 参照)
```

Table 2 実際値:
- Low: Δt = 2.0E-5 s
- Middle: Δt = 9.1E-6 s
- High: Δt = 5.8E-5 s

---

## ジオメトリ設定（三点曲げ）

### Hexagonal Close Packing
- 粒子半径 r: スケールに応じて設定
- ビーム形状: 矩形断面の梁（幅 W, 高さ H, 長さ L）
- スケール比: Low:Middle:High ≈ 1:2:6 程度（粒子数比から推定）

### 境界条件
- 下部2点（x=0, x=L付近）でサポート粒子（固定または変位拘束）
- 中央上部で一定速度/変位の強制荷重（下向き）

### ボンド生成
- 隣接粒子（距離 ≤ 2.1r 程度）を自動ボンド
- 初期方向 d₀ = (pⱼ - pᵢ) / ‖pⱼ - pᵢ‖

---

## プロジェクト構成

```
ibdem/
├── CMakeLists.txt
├── IMPLEMENTATION_PLAN.md
├── PROMPTS.md
├── src/
│   ├── main.cpp                  # エントリーポイント、GLUTループ
│   ├── Simulation.h/.cpp         # シミュレーション全体管理
│   ├── Particle.h                # 粒子データ構造
│   ├── Bond.h/.cpp               # ボンドデータ + エネルギー/勾配/ヘッセ
│   ├── HexPacking.h/.cpp         # Hexagonal Close Packing 生成
│   ├── math/
│   │   ├── Quat.h                # クォータニオン演算 (glmラップ + 独自演算子)
│   │   └── SparseMatrix.h        # 疎行列 (Eigen ライク、自前実装 or Eigen)
│   ├── solver/
│   │   ├── ManifoldOptimizer.h/.cpp  # Algorithm 3
│   │   └── PCG.h/.cpp               # 前処理共役勾配法
│   └── renderer/
│       ├── Renderer.h/.cpp       # OpenGL 描画
│       └── Camera.h/.cpp         # カメラ操作
└── third_party/                  # (オプション) Eigen ヘッダのみ
```

---

## CMake クロスプラットフォーム設計

```cmake
cmake_minimum_required(VERSION 3.16)
project(ibdem CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(glm REQUIRED)

if(APPLE)
    # macOS: Framework GLUT (Apple 標準)
    find_library(OPENGL_LIB OpenGL)
    find_library(GLUT_LIB GLUT)
else()
    # Windows / Ubuntu: freeglut
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)   # freeglut が CMake GLUT として見える
endif()
```

### 依存ライブラリ

| OS       | OpenGL       | GLUT            | GLM                    |
|----------|-------------|-----------------|------------------------|
| Windows  | 標準         | freeglut (vcpkg or manual) | vcpkg or header-only |
| macOS    | Framework    | Apple GLUT (Framework) | Homebrew or header-only |
| Ubuntu   | libGL        | freeglut-dev    | libglm-dev             |

### インストール手順
- **Windows**: `vcpkg install freeglut glm`
- **macOS**: `brew install glm` (GLUT は macOS 標準)
- **Ubuntu**: `sudo apt install freeglut3-dev libglm-dev`

---

## 可視化方針（OpenGL + freeglut/GLUT）

### 表示内容
1. **3ビューポート** (Low / Middle / High) を横並び表示
2. 各ビューポートに:
   - 粒子球（アクティブ粒子: 灰色 / サポート粒子: 青 / 荷重粒子: 赤）
   - ボンド（アクティブ: 緑 / 破断: 非表示 or 赤）
   - 地面ライン
3. フレーム番号・破断フレームを HUD 表示

### カメラ
- 固定斜め視点（3点曲げを側面から俯瞰）
- マウスドラッグで回転（オプション）

### アニメーション
- `glutTimerFunc` で物理ステップ → 再描画
- スペースキーで一時停止/再開
- `r` キーでリセット

---

## 実装フェーズ

### Phase 1: 基盤・粒子生成（優先度: 高）
- [ ] CMakeLists.txt クロスプラットフォーム設定
- [ ] `Particle.h` / `Bond.h` データ構造
- [ ] `HexPacking.cpp` — 矩形梁のHCP充填
- [ ] OpenGL ウィンドウ + 球描画（静的表示確認）

### Phase 2: エネルギー・勾配・ヘッセ（優先度: 高）
- [ ] `Bond.cpp` — V_stretch, V_shear, V_bendtwist の実装
- [ ] クォータニオン演算ユーティリティ (`Quat.h`)
- [ ] 勾配 ∇V（解析微分）
- [ ] ヘッセ行列 ∇²V（解析微分）
- [ ] 疎行列アセンブリ

### Phase 3: ソルバー（優先度: 高）
- [ ] PCG ソルバー（多様体ヌル空間縮約）
- [ ] SPD 射影（対角シフト or 負固有値クリッピング）
- [ ] 多様体上の直線探索
- [ ] Algorithm 3 メインループ
- [ ] 速度更新・摩擦コーン投影

### Phase 4: 破断・衝突（優先度: 中）
- [ ] 破断判定（σ, τ 計算）
- [ ] ボンド削除
- [ ] 簡易反発ポテンシャル（自己衝突・床）

### Phase 5: 可視化・UI（優先度: 中）
- [ ] 3分割ビューポート
- [ ] アニメーションループ（glutTimerFunc）
- [ ] キーボード操作（スペース/r/ESC）
- [ ] フレームカウンタ HUD

### Phase 6: 検証（優先度: 高）
- [ ] Low スケール破断フレーム ≈ 23
- [ ] Middle スケール破断フレーム ≈ 25
- [ ] High スケール破断フレーム ≈ 25
- [ ] スケール間で破断タイミングが一致することを確認

---

## 重要な実装上の注意点

1. **クォータニオン規約**: 論文は `q = {u,v,w,s}` (xyz虚部, w実部)。glm は `{x,y,z,w}` と同じ。
2. **質量行列の回転成分**: `(8/5)mr²` は球の慣性モーメント `I = (2/5)mr²` から `4I = (8/5)mr²`。
3. **ヘッセのSPD射影**: 各ボンドのブロックごとに固有値分解し負値を0にクリップ、またはグローバルに対角シフト。
4. **ヌル空間縮約**: `G(q)` は 3×4 行列。縮約後の系は各粒子あたり 3+3=6 自由度（位置3+回転3）。
5. **HCP充填**: 粒子中心間隔 = 2r。六方最密充填の格子ベクトルを用いる。
6. **境界条件**: サポート粒子の位置を固定（勾配・速度をゼロにクランプ）。

---

## 参考: 論文中の記号対応

| 記号     | 意味                          |
|---------|-------------------------------|
| p       | 粒子位置 (R³)                  |
| q       | 粒子回転クォータニオン (S³)     |
| l₀      | ボンド初期長                  |
| d₀      | ボンド初期方向                |
| q₀      | d₀ を {1,0,0} に移す回転      |
| Gₗ, Gᵣ  | 左・右ヌル空間演算子 (3×4)     |
| Qₗ, Qᵣ  | 左・右クォータニオン乗算行列 (4×4) |
| E       | ヤング率                      |
| G       | せん断弾性係数 = E/(2(1+ν))   |
| ν       | ポアソン比                    |
| σ_c, τ_c | 引張・せん断強度             |
| Δt      | タイムステップ                |
| ε       | 収束判定閾値                  |
| εᵣ      | PCG 相対誤差 = min(0.5,√‖grad f‖) |
