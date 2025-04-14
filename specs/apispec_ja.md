# OpenMM Python API 仕様書

このドキュメントはOpenMMのPython APIの概要、主要クラス、使用方法、およびユースケースについて説明しています。各関数・メソッドの引数の詳細も含まれています。

## 目次

1. [イントロダクション](#イントロダクション)
2. [モジュール構造](#モジュール構造)
3. [主要コンポーネント](#主要コンポーネント)
   - [System](#system)
   - [Force](#force)
   - [Integrator](#integrator)
   - [Context](#context)
   - [Platform](#platform)
   - [Modeller](#modeller-openmmappmodeller)
4. [シミュレーションの設定](#シミュレーションの設定)
5. [ファイル入出力](#ファイル入出力)
6. [解析と可視化](#解析と可視化)
7. [高度な機能](#高度な機能)
8. [コード例](#コード例)
9. [トラブルシューティング](#トラブルシューティング)
10. [関数リファレンス](#関数リファレンス)

## イントロダクション

OpenMMは分子シミュレーションのためのツールキットです。スタンドアロンアプリケーションとして使用することも、自分のコードから呼び出すライブラリとしても使用できます。カスタムフォースやインテグレータを通じた極端な柔軟性、オープン性、高いパフォーマンス（特に最近のGPUにおいて）を組み合わせて提供し、シミュレーションコードの中でも真にユニークなものとなっています。

## モジュール構造

OpenMMのPython APIは、以下の主要なモジュールで構成されています：

- `openmm` - コアライブラリへのインターフェース
- `openmm.app` - 分子システムの構築、ファイル入出力、シミュレーション実行のためのツール
- `openmm.unit` - 単位付き物理量の計算

## 主要コンポーネント

### System

`System`クラスは分子系の物理的記述を表します。粒子（原子）、力場（forces）、および制約で構成されます。

#### 主要メソッド

**addParticle**
```python
System.addParticle(mass)
```
- `mass` (量子化された質量) - 粒子の質量（単位付き）

**addConstraint**
```python
System.addConstraint(particle1, particle2, distance)
```
- `particle1` (int) - 最初の粒子のインデックス
- `particle2` (int) - 2番目の粒子のインデックス
- `distance` (量子化された長さ) - 制約される距離（単位付き）

**addForce**
```python
System.addForce(force)
```
- `force` (Force) - システムに追加する力場オブジェクト

**getNumParticles**
```python
System.getNumParticles()
```
- 戻り値: システム内の粒子数

**getNumConstraints**
```python
System.getNumConstraints()
```
- 戻り値: システム内の制約数

**getNumForces**
```python
System.getNumForces()
```
- 戻り値: システム内の力場の数

**getForce**
```python
System.getForce(index)
```
- `index` (int) - 力場のインデックス
- 戻り値: 指定されたインデックスの力場オブジェクト

**setDefaultPeriodicBoxVectors**
```python
System.setDefaultPeriodicBoxVectors(a, b, c)
```
- `a`, `b`, `c` (Vec3) - 周期的境界条件の箱ベクトル（単位付き）

### Force

`Force`は分子系内の相互作用を表す抽象基底クラスです。OpenMMには多数の組み込み力場クラスがあります：

#### HarmonicBondForce

```python
force = HarmonicBondForce()
```

**addBond**
```python
HarmonicBondForce.addBond(particle1, particle2, length, k)
```
- `particle1` (int) - 最初の粒子のインデックス
- `particle2` (int) - 2番目の粒子のインデックス
- `length` (量子化された長さ) - 平衡結合長（単位付き）
- `k` (量子化されたエネルギー/長さ^2) - 力の定数（単位付き）

#### HarmonicAngleForce

```python
force = HarmonicAngleForce()
```

**addAngle**
```python
HarmonicAngleForce.addAngle(particle1, particle2, particle3, angle, k)
```
- `particle1` (int) - 最初の粒子のインデックス
- `particle2` (int) - 中央の粒子のインデックス
- `particle3` (int) - 3番目の粒子のインデックス
- `angle` (量子化された角度) - 平衡角度（単位付き、ラジアン）
- `k` (量子化されたエネルギー/角度^2) - 力の定数（単位付き）

#### PeriodicTorsionForce

```python
force = PeriodicTorsionForce()
```

**addTorsion**
```python
PeriodicTorsionForce.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
```
- `particle1`, `particle2`, `particle3`, `particle4` (int) - 4つの粒子のインデックス
- `periodicity` (int) - 周期性
- `phase` (量子化された角度) - 位相角度（単位付き、ラジアン）
- `k` (量子化されたエネルギー) - 力の定数（単位付き）

#### NonbondedForce

```python
force = NonbondedForce()
```

**addParticle**
```python
NonbondedForce.addParticle(charge, sigma, epsilon)
```
- `charge` (量子化された電荷) - 粒子の電荷（単位付き）
- `sigma` (量子化された長さ) - Lennard-Jonesシグマパラメータ（単位付き）
- `epsilon` (量子化されたエネルギー) - Lennard-Jonesイプシロンパラメータ（単位付き）

**addException**
```python
NonbondedForce.addException(particle1, particle2, chargeProd, sigma, epsilon, replace=False)
```
- `particle1`, `particle2` (int) - 2つの粒子のインデックス
- `chargeProd` (量子化された電荷^2) - 電荷の積（単位付き）
- `sigma` (量子化された長さ) - Lennard-Jonesシグマパラメータ（単位付き）
- `epsilon` (量子化されたエネルギー) - Lennard-Jonesイプシロンパラメータ（単位付き）
- `replace` (bool) - 既存の例外を置き換えるかどうか

**setNonbondedMethod**
```python
NonbondedForce.setNonbondedMethod(method)
```
- `method` (int) - 非結合相互作用の計算方法（`NoCutoff`, `CutoffNonPeriodic`, `CutoffPeriodic`, `Ewald`, `PME`, `LJPME`）

**setCutoffDistance**
```python
NonbondedForce.setCutoffDistance(distance)
```
- `distance` (量子化された長さ) - カットオフ距離（単位付き）

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

**addPerBondParameter**
```python
CustomBondForce.addPerBondParameter(name)
```
- `name` (str) - パラメータ名

**addBond**
```python
CustomBondForce.addBond(particle1, particle2, parameters)
```
- `particle1`, `particle2` (int) - 2つの粒子のインデックス
- `parameters` (list) - パラメータ値のリスト

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

**addPerParticleParameter**
```python
CustomNonbondedForce.addPerParticleParameter(name)
```
- `name` (str) - パラメータ名

**addParticle**
```python
CustomNonbondedForce.addParticle(parameters)
```
- `parameters` (list) - パラメータ値のリスト

**addExclusion**
```python
CustomNonbondedForce.addExclusion(particle1, particle2)
```
- `particle1`, `particle2` (int) - 排除する2つの粒子のインデックス

### Integrator

`Integrator`クラスは運動方程式の数値積分を実装します。主なインテグレータには以下のものがあります：

#### VerletIntegrator

```python
integrator = VerletIntegrator(stepSize)
```
- `stepSize` (量子化された時間) - 積分のタイムステップ（単位付き）

#### LangevinIntegrator

```python
integrator = LangevinIntegrator(temperature, frictionCoeff, stepSize)
```
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frictionCoeff` (量子化された逆時間) - 摩擦係数（単位付き）
- `stepSize` (量子化された時間) - 積分のタイムステップ（単位付き）

#### BrownianIntegrator

```python
integrator = BrownianIntegrator(temperature, frictionCoeff, stepSize)
```
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frictionCoeff` (量子化された逆時間) - 摩擦係数（単位付き）
- `stepSize` (量子化された時間) - 積分のタイムステップ（単位付き）

#### VariableVerletIntegrator

```python
integrator = VariableVerletIntegrator(errorTol)
```
- `errorTol` (float) - 許容誤差（無次元）

#### VariableLangevinIntegrator

```python
integrator = VariableLangevinIntegrator(temperature, frictionCoeff, errorTol)
```
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frictionCoeff` (量子化された逆時間) - 摩擦係数（単位付き）
- `errorTol` (float) - 許容誤差（無次元）

#### MTSIntegrator

```python
integrator = MTSIntegrator(timestep, loops)
```
- `timestep` (list) - 階層的なタイムステップのリスト（単位付き）
- `loops` (list) - 各レベルで行うループの数のリスト

#### AMDIntegrator

```python
integrator = AMDIntegrator(temperature, frictionCoeff, stepSize, alpha, E)
```
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frictionCoeff` (量子化された逆時間) - 摩擦係数（単位付き）
- `stepSize` (量子化された時間) - 積分のタイムステップ（単位付き）
- `alpha` (量子化されたエネルギー) - ブースト・パラメータ（単位付き）
- `E` (量子化されたエネルギー) - エネルギーしきい値（単位付き）

### Context

`Context`はシミュレーションの実行状態を保持します。特定のSystem、Integrator、およびPlatformを組み合わせて作成されます。

```python
context = Context(system, integrator, platform=None, properties=None)
```
- `system` (System) - シミュレーションするシステム
- `integrator` (Integrator) - 使用するインテグレータ
- `platform` (Platform, optional) - 使用するプラットフォーム
- `properties` (dict, optional) - プラットフォーム固有のプロパティ

#### 主要メソッド

**getState**
```python
Context.getState(getPositions=False, getVelocities=False, getForces=False, 
                getEnergy=False, getParameters=False, getParameterDerivatives=False, 
                getIntegratorParameters=False, enforcePeriodicBox=False, groups=-1)
```
- `getPositions` (bool) - 位置を含めるかどうか
- `getVelocities` (bool) - 速度を含めるかどうか
- `getForces` (bool) - 力を含めるかどうか
- `getEnergy` (bool) - エネルギーを含めるかどうか
- `getParameters` (bool) - コンテキストパラメータを含めるかどうか
- `getParameterDerivatives` (bool) - パラメータ導関数を含めるかどうか
- `getIntegratorParameters` (bool) - インテグレータパラメータを含めるかどうか
- `enforcePeriodicBox` (bool) - 位置を周期的境界条件内に強制するかどうか
- `groups` (int) - 含めるフォースグループ（ビットマスク）

**setState**
```python
Context.setState(state)
```
- `state` (State) - 設定する状態

**reinitialize**
```python
Context.reinitialize(preserveState=False)
```
- `preserveState` (bool) - 現在の状態を保持するかどうか

**setPositions**
```python
Context.setPositions(positions)
```
- `positions` (list) - 粒子の位置のリスト（単位付き）

**setVelocities**
```python
Context.setVelocities(velocities)
```
- `velocities` (list) - 粒子の速度のリスト（単位付き）

**setParameter**
```python
Context.setParameter(name, value)
```
- `name` (str) - パラメータ名
- `value` (float) - パラメータ値

**getParameter**
```python
Context.getParameter(name)
```
- `name` (str) - パラメータ名
- 戻り値: パラメータ値

**setPeriodicBoxVectors**
```python
Context.setPeriodicBoxVectors(a, b, c)
```
- `a`, `b`, `c` (Vec3) - 周期的境界条件の箱ベクトル（単位付き）

### Platform

`Platform`クラスはシミュレーションを実行するハードウェアプラットフォーム（CPU、CUDA、OpenCL、参照実装）を表します。

#### 主要メソッド

**getName**
```python
Platform.getName()
```
- 戻り値: プラットフォーム名

**getPropertyNames**
```python
Platform.getPropertyNames()
```
- 戻り値: このプラットフォームで定義されているプロパティ名のリスト

**getPropertyValue**
```python
Platform.getPropertyValue(context, property)
```
- `context` (Context) - コンテキスト
- `property` (str) - プロパティ名
- 戻り値: プロパティの値

**setPropertyValue**
```python
Platform.setPropertyValue(context, property, value)
```
- `context` (Context) - コンテキスト
- `property` (str) - プロパティ名
- `value` (str) - 設定する値

#### 静的メソッド

**getPlatformByName**
```python
Platform.getPlatformByName(name)
```
- `name` (str) - プラットフォーム名
- 戻り値: 名前付きプラットフォーム

**getNumPlatforms**
```python
Platform.getNumPlatforms()
```
- 戻り値: 利用可能なプラットフォームの数

**findPlatform**
```python
Platform.findPlatform(capabilities)
```
- `capabilities` (list) - 必要な機能のリスト
- 戻り値: 条件を満たすプラットフォーム

### Modeller (openmm.app.Modeller)

`Modeller`クラスは分子モデルの編集ツールを提供し、水の追加や水素原子の追加などの操作を可能にします。

```python
modeller = Modeller(topology, positions)
```
- `topology` (Topology) - 初期分子系のトポロジー
- `positions` (list) - 初期原子位置（単位付き）

#### 基本的なプロパティとメソッド

**topology**
```python
modeller.topology
```
- システムの構造を記述するTopologyオブジェクト

**positions**
```python
modeller.positions
```
- 原子位置のリスト（単位付き）

**getTopology**
```python
modeller.getTopology()
```
- 戻り値: モデルのTopology

**getPositions**
```python
modeller.getPositions()
```
- 戻り値: 原子位置のリスト（単位付き）

#### 分子構造の編集

**add**
```python
modeller.add(addTopology, addPositions)
```
- `addTopology` (Topology) - 追加する分子系のトポロジー
- `addPositions` (list) - 追加する分子系の原子位置（単位付き）
- 説明: 既存のモデルに新しい分子構造を追加します

**delete**
```python
modeller.delete(toDelete)
```
- `toDelete` (list) - 削除する原子、残基、またはチェーンのリスト
- 説明: 指定した原子、残基、チェーンをモデルから削除します

**deleteWater**
```python
modeller.deleteWater()
```
- 説明: モデルからすべての水分子を削除します

**convertWater**
```python
modeller.convertWater(model='tip3p')
```
- `model` (str) - 変換先の水モデル（'tip3p', 'tip4pew', 'tip5p', 'spce', 'swm4ndp'など）
- 説明: 水分子を指定したモデルに変換します。異なる仮想サイト構造を持つモデル間の変換をサポートします

#### 溶媒和とイオン添加

**addSolvent**
```python
modeller.addSolvent(forcefield, model='tip3p', boxSize=None, boxVectors=None, 
                   padding=None, numAdded=None, boxShape='cube', 
                   positiveIon='Na+', negativeIon='Cl-', 
                   ionicStrength=0*molar, neutralize=True, residueTemplates=dict())
```
- `forcefield` (ForceField) - 使用する力場
- `model` (str) - 使用する水モデル（'tip3p', 'spce', 'tip4pew', 'tip5p', 'swm4ndp'）
- `boxSize` (Vec3, optional) - 直方体ボックスのサイズ（単位付き）
- `boxVectors` (tuple, optional) - 周期的境界条件の箱ベクトル
- `padding` (quantity, optional) - 溶質の周りの水のパディング距離（単位付き）
- `numAdded` (int, optional) - 追加する溶媒分子（水+イオン）の総数
- `boxShape` (str) - ボックスの形状（'cube', 'dodecahedron', 'octahedron'）
- `positiveIon` (str) - 正イオンのタイプ（'Cs+', 'K+', 'Li+', 'Na+', 'Rb+'）
- `negativeIon` (str) - 負イオンのタイプ（'Cl-', 'Br-', 'F-', 'I-'）
- `ionicStrength` (quantity) - イオン強度（単位付き、molar）
- `neutralize` (bool) - システムを中和するかどうか
- `residueTemplates` (dict) - 特定の残基に使用するテンプレートの指定
- 説明: システムに水分子とイオンを追加して溶媒和します。様々なボックス形状と水モデルをサポートします

#### 水素原子の追加

**addHydrogens**
```python
modeller.addHydrogens(forcefield=None, pH=7.0, variants=None, platform=None, residueTemplates=dict())
```
- `forcefield` (ForceField, optional) - 使用する力場
- `pH` (float) - pHの値（プロトン化状態の決定に使用）
- `variants` (list, optional) - 各残基のバリアント指定（HIE, HID, HIPなど）
- `platform` (Platform, optional) - 水素原子位置計算に使用するプラットフォーム
- `residueTemplates` (dict) - 特定の残基に使用するテンプレートの指定
- 説明: モデルに水素原子を追加します。pHに基づいてプロトン化状態を自動的に選択できます
- 戻り値: 各残基に選択されたバリアントのリスト

#### 特殊粒子の操作

**addExtraParticles**
```python
modeller.addExtraParticles(forcefield, ignoreExternalBonds=False, residueTemplates=dict())
```
- `forcefield` (ForceField) - 使用する力場
- `ignoreExternalBonds` (bool) - 外部結合を無視するかどうか
- `residueTemplates` (dict) - 特定の残基に使用するテンプレートの指定
- 説明: 力場で定義された仮想サイトなどの追加粒子をモデルに追加します

#### 膜システムの構築

**addMembrane**
```python
modeller.addMembrane(forcefield, lipidType='POPC', membraneCenterZ=0*nanometer, 
                    minimumPadding=1*nanometer, positiveIon='Na+', negativeIon='Cl-',
                    ionicStrength=0*molar, neutralize=True, residueTemplates=dict(), platform=None)
```
- `forcefield` (ForceField) - 使用する力場
- `lipidType` (str) - 使用する脂質タイプ（例：'POPC'）
- `membraneCenterZ` (quantity) - 膜のZ座標中心位置（単位付き）
- `minimumPadding` (quantity) - 膜端と溶質分子間の最小距離（単位付き）
- `positiveIon`, `negativeIon` (str) - 使用するイオンのタイプ
- `ionicStrength` (quantity) - イオン強度（単位付き、molar）
- `neutralize` (bool) - システムを中和するかどうか
- `residueTemplates` (dict) - 特定の残基に使用するテンプレートの指定
- `platform` (Platform, optional) - エネルギー最小化に使用するプラットフォーム
- 説明: 分子系に脂質二重膜を追加します。タンパク質膜複合体の構築に有用です

#### ヘルパーメソッド

**loadHydrogenDefinitions**
```python
Modeller.loadHydrogenDefinitions(file)
```
- `file` (str) - 水素定義XMLファイルのパス
- 説明: カスタム水素定義を読み込みます（クラスメソッド）

### Force

`Force`は分子系内の相互作用を表す抽象基底クラスです。OpenMMには多数の組み込み力場クラスがあります：

#### HarmonicBondForce

```python
force = HarmonicBondForce()
```

**addBond**
```python
HarmonicBondForce.addBond(particle1, particle2, length, k)
```
- `particle1` (int) - 最初の粒子のインデックス
- `particle2` (int) - 2番目の粒子のインデックス
- `length` (量子化された長さ) - 平衡結合長（単位付き）
- `k` (量子化されたエネルギー/長さ^2) - 力の定数（単位付き）

#### HarmonicAngleForce

```python
force = HarmonicAngleForce()
```

**addAngle**
```python
HarmonicAngleForce.addAngle(particle1, particle2, particle3, angle, k)
```
- `particle1` (int) - 最初の粒子のインデックス
- `particle2` (int) - 中央の粒子のインデックス
- `particle3` (int) - 3番目の粒子のインデックス
- `angle` (量子化された角度) - 平衡角度（単位付き、ラジアン）
- `k` (量子化されたエネルギー/角度^2) - 力の定数（単位付き）

#### PeriodicTorsionForce

```python
force = PeriodicTorsionForce()
```

**addTorsion**
```python
PeriodicTorsionForce.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
```
- `particle1`, `particle2`, `particle3`, `particle4` (int) - 4つの粒子のインデックス
- `periodicity` (int) - 周期性
- `phase` (量子化された角度) - 位相角度（単位付き、ラジアン）
- `k` (量子化されたエネルギー) - 力の定数（単位付き）

#### NonbondedForce

```python
force = NonbondedForce()
```

**addParticle**
```python
NonbondedForce.addParticle(charge, sigma, epsilon)
```
- `charge` (量子化された電荷) - 粒子の電荷（単位付き）
- `sigma` (量子化された長さ) - Lennard-Jonesシグマパラメータ（単位付き）
- `epsilon` (量子化されたエネルギー) - Lennard-Jonesイプシロンパラメータ（単位付き）

**addException**
```python
NonbondedForce.addException(particle1, particle2, chargeProd, sigma, epsilon, replace=False)
```
- `particle1`, `particle2` (int) - 2つの粒子のインデックス
- `chargeProd` (量子化された電荷^2) - 電荷の積（単位付き）
- `sigma` (量子化された長さ) - Lennard-Jonesシグマパラメータ（単位付き）
- `epsilon` (量子化されたエネルギー) - Lennard-Jonesイプシロンパラメータ（単位付き）
- `replace` (bool) - 既存の例外を置き換えるかどうか

**setNonbondedMethod**
```python
NonbondedForce.setNonbondedMethod(method)
```
- `method` (int) - 非結合相互作用の計算方法（`NoCutoff`, `CutoffNonPeriodic`, `CutoffPeriodic`, `Ewald`, `PME`, `LJPME`）

**setCutoffDistance**
```python
NonbondedForce.setCutoffDistance(distance)
```
- `distance` (量子化された長さ) - カットオフ距離（単位付き）

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

**addPerBondParameter**
```python
CustomBondForce.addPerBondParameter(name)
```
- `name` (str) - パラメータ名

**addBond**
```python
CustomBondForce.addBond(particle1, particle2, parameters)
```
- `particle1`, `particle2` (int) - 2つの粒子のインデックス
- `parameters` (list) - パラメータ値のリスト

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

**addPerParticleParameter**
```python
CustomNonbondedForce.addPerParticleParameter(name)
```
- `name` (str) - パラメータ名

**addParticle**
```python
CustomNonbondedForce.addParticle(parameters)
```
- `parameters` (list) - パラメータ値のリスト

**addExclusion**
```python
CustomNonbondedForce.addExclusion(particle1, particle2)
```
- `particle1`, `particle2` (int) - 排除する2つの粒子のインデックス

## シミュレーションの設定

OpenMMのシミュレーションを設定する最も一般的な方法は、`openmm.app.Simulation`クラスを使用することです。

```python
simulation = Simulation(topology, system, integrator, platform=None, platformProperties=None, state=None)
```
- `topology` (Topology) - 分子系のトポロジー
- `system` (System) - シミュレーションするシステム
- `integrator` (Integrator) - 使用するインテグレータ
- `platform` (Platform, optional) - 使用するプラットフォーム
- `platformProperties` (dict, optional) - プラットフォーム固有のプロパティ
- `state` (State, optional) - 初期状態

### 主要メソッド

**context**
```python
Simulation.context
```
- シミュレーションのContextオブジェクト

**minimizeEnergy**
```python
Simulation.minimizeEnergy(tolerance=10*kilojoules/mole, maxIterations=0)
```
- `tolerance` (量子化されたエネルギー) - 収束許容誤差（単位付き）
- `maxIterations` (int) - 最大反復回数（0は無制限）

**step**
```python
Simulation.step(steps)
```
- `steps` (int) - 実行するステップ数

**saveState**
```python
Simulation.saveState(file)
```
- `file` (str または file-like object) - 状態を保存するファイル

**loadState**
```python
Simulation.loadState(file)
```
- `file` (str または file-like object) - 状態を読み込むファイル

**reporters**
```python
Simulation.reporters
```
- レポーターのリスト

シミュレーションの基本的なセットアップは以下のステップで行います：

1. トポロジー（`Topology`）と位置情報を定義する
2. 力場（`ForceField`）を選択する
3. システム（`System`）を作成する
4. インテグレータ（`Integrator`）を選択する
5. シミュレーション（`Simulation`）オブジェクトを作成する
6. 位置と速度を設定する
7. レポーターを追加する
8. シミュレーションを実行する

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# 1. PDBファイルからトポロジーと位置情報を読み込む
pdb = PDBFile('input.pdb')

# 2. 力場ファイルを指定
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# 3. システムを作成
system = forcefield.createSystem(
    pdb.topology, 
    nonbondedMethod=PME,  # 周期的境界条件を使用
    nonbondedCutoff=1*nanometer,  # 非結合相互作用のカットオフ
    constraints=HBonds    # 水素結合の長さを固定
)

# 4. インテグレータを選択（ここではランジュバン動力学）
integrator = LangevinIntegrator(
    300*kelvin,       # 温度
    1/picosecond,     # 摩擦係数
    0.002*picoseconds # タイムステップ
)

# 5. シミュレーションを作成
simulation = Simulation(pdb.topology, system, integrator)

# 6. 初期位置を設定
simulation.context.setPositions(pdb.positions)

# 7. レポーターを追加
simulation.reporters.append(PDBReporter('output.pdb', 1000))  # 1000ステップごとに構造を出力
simulation.reporters.append(StateDataReporter(
    'output.csv', 1000, step=True, temperature=True, potentialEnergy=True))

# 8. エネルギー最小化を実行
simulation.minimizeEnergy()

# 9. シミュレーションを実行（例：100psの分子動力学）
simulation.step(50000)  # 0.002 ps × 50,000 steps = 100 ps
```

## ファイル入出力

OpenMMは様々な分子ファイル形式の読み込みと書き込みをサポートしています：

### 入力ファイル

#### PDBFile

```python
pdb = PDBFile(file)
```
- `file` (str または file-like object) - PDBファイルのパスまたはオブジェクト

#### AmberPrmtopFile / AmberInpcrdFile

```python
prmtop = AmberPrmtopFile(file)
inpcrd = AmberInpcrdFile(file)
```
- `file` (str または file-like object) - AMBERファイルのパスまたはオブジェクト

#### GromacsTopFile / GromacsGroFile

```python
top = GromacsTopFile(file, periodicBoxVectors=None, includeDir=None)
gro = GromacsGroFile(file)
```
- `file` (str または file-like object) - GROMACSファイルのパスまたはオブジェクト
- `periodicBoxVectors` (tuple, optional) - 周期的境界条件の箱ベクトル
- `includeDir` (str, optional) - インクルードファイルのディレクトリ

#### CharmmPsfFile / CharmmCrdFile

```python
psf = CharmmPsfFile(file)
crd = CharmmCrdFile(file)
```
- `file` (str または file-like object) - CHARMMファイルのパスまたはオブジェクト

### 出力/レポーター

#### PDBReporter

```python
reporter = PDBReporter(file, reportInterval, enforcePeriodicBox=True)
```
- `file` (str または file-like object) - 出力ファイルのパスまたはオブジェクト
- `reportInterval` (int) - 報告間隔（ステップ数）
- `enforcePeriodicBox` (bool) - 粒子を周期的箱内に収めるかどうか

#### DCDReporter

```python
reporter = DCDReporter(file, reportInterval, append=False, enforcePeriodicBox=True)
```
- `file` (str または file-like object) - 出力ファイルのパスまたはオブジェクト
- `reportInterval` (int) - 報告間隔（ステップ数）
- `append` (bool) - 既存のファイルに追加するかどうか
- `enforcePeriodicBox` (bool) - 粒子を周期的箱内に収めるかどうか

#### XTCReporter

```python
reporter = XTCReporter(file, reportInterval, append=False, enforcePeriodicBox=True)
```
- `file` (str または file-like object) - 出力ファイルのパスまたはオブジェクト
- `reportInterval` (int) - 報告間隔（ステップ数）
- `append` (bool) - 既存のファイルに追加するかどうか
- `enforcePeriodicBox` (bool) - 粒子を周期的箱内に収めるかどうか

#### StateDataReporter

```python
reporter = StateDataReporter(file, reportInterval, step=False, time=False, 
                            potentialEnergy=False, kineticEnergy=False, 
                            totalEnergy=False, temperature=False, volume=False, 
                            density=False, progress=False, remainingTime=False, 
                            speed=False, elapsedTime=False, separator=',', 
                            systemMass=None, totalSteps=None)
```
- `file` (str, file-like object またはNone) - 出力ファイルまたはコンソール出力（None）
- `reportInterval` (int) - 報告間隔（ステップ数）
- 他の引数 - 報告する物理量や出力形式を指定

#### CheckpointReporter

```python
reporter = CheckpointReporter(file, reportInterval)
```
- `file` (str) - 出力ファイルのパス
- `reportInterval` (int) - 報告間隔（ステップ数）

## 解析と可視化

OpenMMでは、シミュレーション結果の解析に使用できるツールがいくつか提供されています：

- `StateDataReporter`を使用して、エネルギー、温度などの物理量の時系列を記録
- トラジェクトリをPDB、DCD、XTC形式で出力し、VMD、PyMOL、MDAnalysisなどの外部ツールで解析
- シミュレーションの`State`オブジェクトから位置、速度、力などの情報を直接取得して分析

### State オブジェクト

`State`オブジェクトはシミュレーションの瞬間的な状態を表します。

```python
state = context.getState(getPositions=False, getVelocities=False, getForces=False, 
                        getEnergy=False, getParameters=False, getParameterDerivatives=False, 
                        enforcePeriodicBox=False, groups=-1)
```

#### 主要メソッド

**getPositions**
```python
State.getPositions(asNumpy=False)
```
- `asNumpy` (bool) - NumPy配列として返すかどうか
- 戻り値: 粒子の位置（単位付き）

**getVelocities**
```python
State.getVelocities(asNumpy=False)
```
- `asNumpy` (bool) - NumPy配列として返すかどうか
- 戻り値: 粒子の速度（単位付き）

**getForces**
```python
State.getForces(asNumpy=False)
```
- `asNumpy` (bool) - NumPy配列として返すかどうか
- 戻り値: 粒子に作用する力（単位付き）

**getPotentialEnergy**
```python
State.getPotentialEnergy()
```
- 戻り値: ポテンシャルエネルギー（単位付き）

**getKineticEnergy**
```python
State.getKineticEnergy()
```
- 戻り値: 運動エネルギー（単位付き）

**getPeriodicBoxVectors**
```python
State.getPeriodicBoxVectors()
```
- 戻り値: 周期的境界条件の箱ベクトル（単位付き）

## 高度な機能

### カスタム力場

OpenMMの最も強力な機能の一つは、カスタム力場を定義する能力です。これらを使用すると、あらゆる物理的相互作用をほぼ表現できます：

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム結合力は、2つの原子間の任意の関数形式のポテンシャルを定義します。数式は変数`r`（結合の長さ）とユーザー定義のパラメータを使用できます。

#### CustomAngleForce

```python
force = CustomAngleForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム角度力は、3つの原子間の任意の関数形式のポテンシャルを定義します。数式は変数`theta`（角度）とユーザー定義のパラメータを使用できます。

#### CustomTorsionForce

```python
force = CustomTorsionForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム二面角力は、4つの原子間の任意の関数形式のポテンシャルを定義します。数式は変数`theta`（二面角）とユーザー定義のパラメータを使用できます。

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム非結合力は、粒子ペア間の任意の関数形式のポテンシャルを定義します。数式は変数`r`（粒子間の距離）とユーザー定義のパラメータを使用できます。

#### CustomExternalForce

```python
force = CustomExternalForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム外部力は、個々の粒子に作用する任意の関数形式の外部ポテンシャルを定義します。数式は変数`x`、`y`、`z`（粒子の座標）とユーザー定義のパラメータを使用できます。

### 拘束シミュレーション

OpenMMには、さまざまな種類の拘束を適用するためのツールが含まれています：

#### 結合長の拘束

結合長拘束は`System.createSystem()`メソッドの`constraints`パラメータで指定します：

```python
system = forcefield.createSystem(topology, constraints=HBonds)
```

- `constraints` - 拘束のレベル（`None`, `HBonds`, `AllBonds`, `HAngles`）

#### CustomCVForce

```python
force = CustomCVForce(energy)
```
- `energy` (str) - ポテンシャルエネルギーを定義する数式

カスタム集合変数力は、複雑な集合変数（CV）に基づくポテンシャルを定義します。

#### RMSDForce

```python
force = RMSDForce(referencePositions, particles=None)
```
- `referencePositions` (list) - 参照位置のリスト（単位付き）
- `particles` (list, optional) - 含める粒子のインデックスのリスト

RMSDに基づく拘束力は、現在の構造と参照構造の間のRMSD（二乗平均平方根偏差）に基づくポテンシャルを定義します。

### 高度なサンプリング方法

OpenMMには、通常の分子動力学シミュレーションでは効率的に探索できない系のサンプリングを改善するための様々な手法が実装されています。以下に主要な高度なサンプリング手法とその使用例を示します。

#### MTSIntegrator（マルチタイムステップ積分）

```python
integrator = MTSIntegrator(timestep, loops)
```
- `timestep` (list) - 階層的なタイムステップのリスト（単位付き）
- `loops` (list) - 各レベルで行うループの数のリスト

マルチタイムステップ積分器は、異なる力をさまざまな頻度で評価することで計算効率を向上させます。例えば、速く変化する力（結合力など）と遅く変化する力（遠距離相互作用など）を区別して扱います。

**使用例**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# PDBファイルから構造を読み込む
pdb = PDBFile('input.pdb')

# 力場を定義
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# システムを作成し、力を2つのグループに分ける
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# NonbondedForceを見つけて、グループ1に設定
for i in range(system.getNumForces()):
    force = system.getForce(i)
    if isinstance(force, NonbondedForce):
        force.setForceGroup(1)
    else:
        # 他の力はグループ0に設定
        force.setForceGroup(0)

# MTSIntegratorを設定
# グループ0（結合力など）は0.5 fsで評価
# グループ1（非結合力など）は2.0 fsで評価
inner_ts = 0.5*femtoseconds
outer_ts = 2.0*femtoseconds

integrator = MTSIntegrator([inner_ts, outer_ts], [4, 1])
# グループ0はinner_tsで4回評価され、グループ1はouter_tsで1回評価される

# シミュレーションを準備
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# エネルギー最小化
simulation.minimizeEnergy()

# シミュレーション実行
simulation.reporters.append(StateDataReporter('mts_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))
simulation.step(50000)  # 50,000ステップのシミュレーション実行
```

#### AMDIntegrator（加速分子動力学）

```python
integrator = AMDIntegrator(temperature, frictionCoeff, stepSize, alpha, E)
```
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frictionCoeff` (量子化された逆時間) - 摩擦係数（単位付き）
- `stepSize` (量子化された時間) - 積分のタイムステップ（単位付き）
- `alpha` (量子化されたエネルギー) - ブースト・パラメータ（単位付き）
- `E` (量子化されたエネルギー) - エネルギーしきい値（単位付き）

加速分子動力学（AMD）積分器は、ポテンシャルエネルギー地形を修正することで、分子系の構造空間の探索を加速します。特に、エネルギー障壁が高く、通常のシミュレーションでは遷移が起こりにくい系に有効です。

**使用例**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np

# PDBファイルから構造を読み込む
pdb = PDBFile('protein.pdb')

# 力場を定義
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# システムを作成
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# まず標準的なシミュレーションを短時間実行して、エネルギーを見積もる
temp_integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
temp_simulation = Simulation(pdb.topology, system, temp_integrator)
temp_simulation.context.setPositions(pdb.positions)
temp_simulation.minimizeEnergy()
temp_simulation.context.setVelocitiesToTemperature(300*kelvin)
temp_simulation.step(5000)  # 短いシミュレーションを実行

# 平均ポテンシャルエネルギーを計算
state = temp_simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

# AMDパラメータを設定
# E_threshはシステムの平均ポテンシャルエネルギーよりも少し高く設定
E_thresh = potential_energy * 1.05 * kilojoule_per_mole
# alphaはE_threshの約1/5程度に設定
alpha = 0.2 * E_thresh

# AMDIntegratorを使用
amd_integrator = AMDIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, alpha, E_thresh)

# 新しいシミュレーションを準備
simulation = Simulation(pdb.topology, system, amd_integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

# レポーターを追加
simulation.reporters.append(DCDReporter('amd_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('amd_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# シミュレーション実行
simulation.step(500000)  # 500,000ステップのシミュレーション実行（1 ns）
```

#### MetaDynamics（メタダイナミクス）

```python
meta = MetaDynamics(system, variables, temperature, biasFactor, height, frequency)
```
- `system` (System) - シミュレーションするシステム
- `variables` (list) - 集合変数のリスト
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `biasFactor` (float) - バイアス因子
- `height` (量子化されたエネルギー) - Gaussianヒルのエネルギー高さ（単位付き）
- `frequency` (int) - ヒルを追加する頻度（ステップ数）

メタダイナミクスは、すでに訪問した構造空間の領域にエネルギーペナルティを適用することで、分子系のサンプリングを改善します。特に、タンパク質の折りたたみや配座遷移のような複雑な構造変化の研究に有用です。

**使用例**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# PDBファイルから構造を読み込む
pdb = PDBFile('protein.pdb')

# 力場を定義
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# システムを作成
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# 基本インテグレータを作成
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# 特定の原子間距離（集合変数）を定義
# 例: 残基10と残基20の間のCA原子の距離
atoms = list(pdb.topology.atoms())
ca_atoms = [atom.index for atom in atoms if atom.name == 'CA']

# 残基10と残基20のCA原子を見つける
res10_ca = None
res20_ca = None
for atom in atoms:
    if atom.name == 'CA':
        if atom.residue.id == 10:
            res10_ca = atom.index
        elif atom.residue.id == 20:
            res20_ca = atom.index

# 集合変数（CV）を定義
cv_force = CustomBondForce('r')
cv_force.addBond(res10_ca, res20_ca, [])
cv_force.setUsesPeriodicBoundaryConditions(True)

# メタダイナミクス用の変数を設定
cvs = [cv_force]
temperature = 300*kelvin
bias_factor = 10.0  # よく使われる値（5-20）
hill_height = 1.0*kilojoule_per_mole  # ヒルの高さ
hill_frequency = 100  # 100ステップごとにヒルを追加

# メタダイナミクスを初期化
meta = MetaDynamics(system, cvs, temperature, bias_factor, hill_height, hill_frequency)

# シミュレーションを準備
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# エネルギー最小化
simulation.minimizeEnergy()

# レポーターを追加
simulation.reporters.append(DCDReporter('metadynamics_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('metadynamics_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# カスタムレポーターで集合変数の値を記録
class CVReporter(object):
    def __init__(self, file, reportInterval):
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        self._file.write('Step,CV_Value,Bias_Energy\n')
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, True)
    
    def report(self, simulation, state):
        cv_value = meta.getCollectiveVariableValues(simulation.context)[0]
        bias_energy = meta.getBiasEnergy(simulation.context)
        self._file.write(f'{simulation.currentStep},{cv_value},{bias_energy}\n')

simulation.reporters.append(CVReporter('cv_values.csv', 100))

# シミュレーション実行
simulation.step(500000)  # 500,000ステップのシミュレーション実行（1 ns）
```

#### SimulatedTempering（シミュレーテッドテンパリング）

```python
st = SimulatedTempering(system, temperatures, scaleForces=True)
```
- `system` (System) - シミュレーションするシステム
- `temperatures` (list) - 温度のリスト（単位付き）
- `scaleForces` (bool) - 力をスケーリングするかどうか

シミュレーテッドテンパリングは、シミュレーション中に温度を周期的に変更することで、エネルギー障壁を超えられるようにします。これにより、通常の一定温度のシミュレーションでは到達できない構造空間を探索することができます。

**使用例**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np

# PDBファイルから構造を読み込む
pdb = PDBFile('protein.pdb')

# 力場を定義
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# システムを作成
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# 基本インテグレータを作成
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# 温度範囲を設定（例: 300K から 400Kの間で8レベル）
min_temp = 300.0
max_temp = 400.0
num_temps = 8
temperatures = [min_temp + (max_temp - min_temp) * i / (num_temps-1) for i in range(num_temps)]
temperatures = [t*kelvin for t in temperatures]

# シミュレーテッドテンパリングを初期化
st = SimulatedTempering(system, temperatures)

# 自由エネルギー重みを初期化（最初は等しくするか、事前に計算した値を使用）
weights = np.zeros(num_temps)
st.setTemperatureWeights(weights)

# シミュレーションを準備
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# エネルギー最小化
simulation.minimizeEnergy()

# レポーターを追加
simulation.reporters.append(DCDReporter('simtemp_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('simtemp_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# カスタムレポーターで温度の状態を記録
class TemperatureStateReporter(object):
    def __init__(self, file, reportInterval):
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        self._file.write('Step,TemperatureIndex,Temperature,AcceptanceRate\n')
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, True)
    
    def report(self, simulation, state):
        temp_index = st.getCurrentTemperature(simulation.context)
        temp = temperatures[temp_index].value_in_unit(kelvin)
        acceptance = st.getAcceptanceRate(temp_index)
        self._file.write(f'{simulation.currentStep},{temp_index},{temp},{acceptance}\n')

simulation.reporters.append(TemperatureStateReporter('temperature_states.csv', 100))

# 重みの適応的更新を行う関数
def update_weights(simulation, st, interval=50000):
    if simulation.currentStep % interval == 0 and simulation.currentStep > 0:
        counts = st.getStateVisitCounts()
        if min(counts) > 0:  # すべての状態が訪問されたことを確認
            # 対数スケールで状態の訪問頻度に基づいて重みを更新
            weights = [-np.log(count/sum(counts)) for count in counts]
            # 最小値を0に設定
            weights = [w - min(weights) for w in weights]
            st.setTemperatureWeights(weights)
            print(f"Step {simulation.currentStep}: Updated weights: {weights}")

# シミュレーション実行（適応的重み更新を含む）
for i in range(20):  # 合計1,000,000ステップ（2 ns）を20回に分ける
    simulation.step(50000)  # 50,000ステップずつ実行
    update_weights(simulation, st)  # 重みを更新
```

## コード例

### 明示的溶媒系のセットアップ例

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# PDBファイルから構造を読み込む
pdb = PDBFile('protein.pdb')

# 力場を定義
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# タンパク質を水で溶媒和
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=1.0*nanometer, model='tip3p')

# システムを構築
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# 温度と圧力のコントロールを追加（NPTアンサンブル）
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# インテグレータを設定
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# シミュレーションを準備
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# エネルギー最小化
simulation.minimizeEnergy()

# 平衡化
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(10000)  # 20 ps の平衡化

# 本番シミュレーション
simulation.reporters.append(DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter(
    'data.csv', 1000, step=True, time=True,
    potentialEnergy=True, temperature=True, density=True
))
simulation.step(500000)  # 1 ns の本番シミュレーション
```

### カスタム力場の例

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# カスタム二面角ポテンシャルの例
custom_torsion = CustomTorsionForce('k*(1+cos(n*theta-gamma))')
custom_torsion.addPerTorsionParameter('k')
custom_torsion.addPerTorsionParameter('n')
custom_torsion.addPerTorsionParameter('gamma')

# カスタム非結合力の例（modified Lennard-Jones）
custom_nb = CustomNonbondedForce(
    '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)'
)
custom_nb.addPerParticleParameter('sigma')
custom_nb.addPerParticleParameter('epsilon')
```

## トラブルシューティング

### 一般的な問題と解決策

1. **シミュレーションが不安定で発散する**
   - エネルギー最小化を適切に行う
   - タイムステップを小さくする（0.5-1 fs）
   - 拘束条件を見直す（特に水素原子を含む結合）
   - 初期構造の立体障害を確認する

2. **パフォーマンスの問題**
   - 適切なプラットフォーム（CUDA/OpenCL）を使用しているか確認
   - 非結合相互作用のカットオフを適切に設定
   - プラットフォーム固有の最適化オプションを設定

3. **プラットフォーム選択の問題**
   - 特定のプラットフォームを明示的に選択：
     ```python
     platform = Platform.getPlatformByName('CUDA')
     properties = {'CudaPrecision': 'mixed'}
     simulation = Simulation(topology, system, integrator, platform, properties)
     ```

4. **単位関連のエラー**
   - OpenMMでは物理量には常に単位が必要
   - 単位の変換と計算には`openmm.unit`モジュールを使用

## 関数リファレンス

### openmm.unit モジュール

単位付き物理量を扱うためのモジュールです。

```python
from openmm.unit import *

# 長さの単位
nanometer      # ナノメートル（1e-9 m）
angstrom       # オングストローム（1e-10 m）
picometer      # ピコメートル（1e-12 m）
bohr           # ボーア半径（量子力学的な長さの単位）

# 時間の単位
picosecond     # ピコ秒（1e-12 s）
femtosecond    # フェムト秒（1e-15 s）
nanosecond     # ナノ秒（1e-9 s）

# エネルギーの単位
kilojoule_per_mole  # キロジュール/モル（生化学的標準）
kilocalorie_per_mole  # キロカロリー/モル
hartree        # ハートリー（量子化学的エネルギー単位）
electronvolt   # 電子ボルト

# 温度の単位
kelvin         # ケルビン

# 質量の単位
dalton         # ダルトン（原子質量単位）
atomic_mass_unit  # 原子質量単位

# 電荷の単位
elementary_charge  # 素電荷

# 圧力の単位
bar            # バール
atmosphere     # 大気圧

# 単位の変換
length = 5.0 * nanometer  # 5 nmの長さ
length_in_angstrom = length.value_in_unit(angstrom)  # オングストロームに変換
```

### Modeller クラス

`Modeller`クラスは分子系を修正するためのツールを提供します。

```python
modeller = Modeller(topology, positions)
```
- `topology` (Topology) - 分子系のトポロジー
- `positions` (list) - 粒子の位置のリスト（単位付き）

#### 主要メソッド

**addSolvent**
```python
Modeller.addSolvent(forcefield, model='tip3p', boxSize=None, boxVectors=None,
                   padding=None, neutralize=True, positiveIon='Na+',
                   negativeIon='Cl-', ionicStrength=0*molar)
```
- `forcefield` (ForceField) - 使用する力場
- `model` (str) - 水モデル（'tip3p', 'tip4pew', 'tip5p'など）
- `boxSize` (Vec3, optional) - ボックスサイズ（単位付き）
- `boxVectors` (tuple, optional) - 周期的境界条件の箱ベクトル
- `padding` (量子化された長さ, optional) - 溶質の周りの水のパディング（単位付き）
- `neutralize` (bool) - システムを中和するかどうか
- `positiveIon`, `negativeIon` (str) - 使用するイオンの種類
- `ionicStrength` (量子化された濃度) - イオン強度（単位付き）

**addHydrogens**
```python
Modeller.addHydrogens(forcefield=None, pH=7.0, variants=None)
```
- `forcefield` (ForceField, optional) - 使用する力場
- `pH` (float) - pHの値
- `variants` (list, optional) - 残基バリアントのリスト

**delete**
```python
Modeller.delete(atoms)
```
- `atoms` (list) - 削除する原子のリスト

### ForceField クラス

`ForceField`クラスは分子力場の定義と適用を処理します。

```python
forcefield = ForceField(*files)
```
- `files` (str) - 力場定義ファイルのリスト

#### 主要メソッド

**createSystem**
```python
ForceField.createSystem(topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer,
                       constraints=None, rigidWater=True, removeCMMotion=True,
                       hydrogenMass=None, residueTemplates=dict(), verbose=False,
                       ewaldErrorTolerance=0.0005)
```
- `topology` (Topology) - 分子系のトポロジー
- `nonbondedMethod` (int) - 非結合相互作用の計算方法
- `nonbondedCutoff` (量子化された長さ) - 非結合相互作用のカットオフ距離（単位付き）
- `constraints` (int) - 拘束のレベル（`None`, `HBonds`, `AllBonds`, `HAngles`）
- `rigidWater` (bool) - 水分子を剛体として扱うかどうか
- `removeCMMotion` (bool) - 重心の運動を除去するかどうか
- `hydrogenMass` (量子化された質量, optional) - 水素原子の質量（単位付き）
- `residueTemplates` (dict) - 残基テンプレートのマッピング
- `verbose` (bool) - 詳細な出力を表示するかどうか
- `ewaldErrorTolerance` (float) - Ewald法の誤差許容値

**getUnmatchedResidues**
```python
ForceField.getUnmatchedResidues(topology)
```
- `topology` (Topology) - 分子系のトポロジー
- 戻り値: マッチしない残基のリスト

**getMatchingTemplates**
```python
ForceField.getMatchingTemplates(topology, ignoreExternalBonds=False)
```
- `topology` (Topology) - 分子系のトポロジー
- `ignoreExternalBonds` (bool) - 外部結合を無視するかどうか
- 戻り値: マッチするテンプレートのリスト

### MonteCarloBarostat

```python
barostat = MonteCarloBarostat(pressure, temperature, frequency=25)
```
- `pressure` (量子化された圧力) - 目標圧力（単位付き）
- `temperature` (量子化された温度) - シミュレーションの温度（単位付き）
- `frequency` (int) - モンテカルロの試行頻度（ステップ数）

### RMSDForce

```python
force = RMSDForce(referencePositions, particles=None)
```
- `referencePositions` (list) - 参照位置のリスト（単位付き）
- `particles` (list, optional) - 含める粒子のインデックスのリスト

OpenMMの詳細なドキュメントは[公式サイト](http://docs.openmm.org/)で参照できます。