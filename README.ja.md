[English](README.md)

# ibdem — 陰的 BDEM 三点曲げ試験（Fig. 5 再現）

以下の論文における三点曲げスケール一貫性実験の C++ 実装です：

> **"Implicit Bonded Discrete Element Method with Manifold Optimization"**
> ACM Transactions on Graphics (2023), DOI: 10.1145/3592427
> （論文 ID: 3711852）

3 種類のスケール（Low / Middle / High）を同時にシミュレーションし、すべてのスケールが同じ荷重レベル（フレーム ~23–25）で破断することを確認します。これは論文 Table 2 および Fig. 5 の結果に対応します。

ボンドは応力レベルで色付けされます：**青**（低応力）→ **緑** → **赤**（破断直前、σ ≈ τ_c）。

---

## ビルド

### Windows

[vcpkg](https://github.com/microsoft/vcpkg) と Visual Studio 2022 をインストール後：

```powershell
vcpkg install freeglut:x64-windows glm:x64-windows
vcpkg integrate install

mkdir build
cd build
cmake ../ -G "Visual Studio 17" -A x64 -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_INSTALLATION_ROOT/scripts/buildsystems/vcpkg.cmake"
cmake --build . --config Release
```

### macOS

```bash
brew install glm

mkdir build && cd build
cmake ../ -G Xcode
cmake --build . --config Release
```

### Ubuntu / Debian

```bash
sudo apt-get install \
  libxi-dev libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev \
  libxrandr-dev libxxf86vm-dev freeglut3-dev libglm-dev

mkdir build && cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

---

## 実行

### GUI モード（通常）

```
./bin/ibdem
```

3 つのビューポート（Low / Middle / High）が横並びの 1200×450 ウィンドウが開きます。
**シミュレーションは自動的に開始**されます。応力が高まるにつれてボンドの色が青から赤に変化し、フレーム 24 前後で破断したボンドが消えます。

**操作方法：**

| キー | 操作 |
|------|------|
| `SPACE` | 一時停止 / 再開 |
| `N` | 1 フレームだけ進める（一時停止中） |
| `R` | 全スケールを初期状態にリセット |
| `S` | BMP スクリーンショット保存（`capture_frNNNN.bmp`） |
| `Q` / `ESC` | 終了 |

### ヘッドレスモード

ウィンドウなしで実行し、各スケールの破断フレームをコンソールに出力して終了します：

```
./bin/ibdem -headless [N]
```

デフォルト N = 60。期待される出力：

```
Scale Low     : fracture at frame 24  (target ~23)
Scale Middle  : fracture at frame 24  (target ~25)
Scale High    : fracture at frame 24  (target ~25)
```

### キャプチャーモード

全フレームを描画して BMP ファイルとして保存し、完了後に自動終了します：

```
./bin/ibdem -capture [N]
```

デフォルト N = 40。カレントディレクトリに `capture_fr0001.bmp` 〜 `capture_frNNNN.bmp` が生成されます。
ウィンドウを操作せずにアニメーションをオフラインで確認するのに便利です。
※ High スケール（~22K 粒子）は 1 フレームあたり数秒かかるため、40 フレームの保存には数分を要します。

---

## パラメータ

すべてのパラメータは論文 Table 2 に対応しています。

| パラメータ | 値 |
|-----------|---|
| 梁の寸法 | L = 1.0 m, H = 0.15 m, W = 0.12 m |
| ヤング率 E | 1 × 10⁷ Pa |
| ポアソン比 ν | 0.3 |
| 破断しきい値 τ_c | 3 × 10⁴ Pa |
| 荷重速度 | 0.069 m/フレーム |
| 粒子半径 r — Low | 0.018 m（約 504 粒子） |
| 粒子半径 r — Middle | 0.011 m（約 2024 粒子） |
| 粒子半径 r — High | 0.005 m（約 22725 粒子） |

---

## プロジェクト構成

```
ibdem/
├── src/
│   ├── main.cpp          # GLUT ウィンドウ・描画・キー入力・実行モード
│   ├── Simulation.cpp/h  # 陰的ソルバー（PCG + 多様体降下法）
│   ├── BondForce.cpp/h   # ボンド応力・破断判定
│   ├── HexPacking.cpp/h  # 六方最密充填による粒子生成
│   ├── Particle.h        # 粒子状態（位置・向き・フラグ）
│   └── Bond.h            # ボンド幾何・キャッシュ済み応力
├── .github/workflows/    # CI: windows.yml, macos.yml, ubuntu.yml
├── CMakeLists.txt
├── CMakePresets.json
└── README.md
```
