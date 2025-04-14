# OpenMM Documentation MCP Server

OpenMMの分子動力学シミュレーションドキュメントのための検索サーバーです。Model Context Protocol (MCP)に準拠しており、LLMとの連携に最適化されています。

## 概要

このサーバーはOpenMMのドキュメントをベクトル化し、検索可能な形で提供します。ユーザーはクエリを送信することで、関連するドキュメントセクションを検索できます。FAISSベクトルデータベースを使用して高速な検索を実現しています。

OpenMMは生物学的分子系のための高性能シミュレーションライブラリであり、このサーバーを通じてドキュメントの検索が容易になります。特に以下のようなケースに最適です：

- 特定の分子動力学手法に関するドキュメントの検索
- OpenMMの関数やクラスの使用方法の検索
- シミュレーションのパラメータ設定に関する情報の取得

## 機能

- ドキュメント検索: 自然言語クエリに基づく関連ドキュメントの検索
- インデックス情報: ベクトルデータベースの情報を取得
- MCP準拠: AI/LLMとの統合に最適
- 高速検索: FAISSを利用した効率的なベクトル検索
- 柔軟なクエリ: 日本語および英語のクエリに対応

## 技術スタック

このプロジェクトは以下の主要技術を使用しています：

- **FastAPI**: REST APIの実装
- **FAISS**: 効率的なベクトル検索インデックス
- **LangChain**: 文書処理とLLM連携
- **HuggingFace Models**: 文書埋め込み生成
- **Pytest**: テスト自動化

## セットアップ

### 前提条件

- Python 3.9以上
- `uv` パッケージマネージャー (推奨) または `pip`
- 最低8GB RAM（インデックス作成時に16GB以上推奨）

### 必要なパッケージ

主な依存パッケージ：

```
fastapi==0.110.0
uvicorn==0.28.0
faiss-cpu==1.8.0
langchain==0.1.12
langchain-community==0.0.29
langchain-text-splitters==0.0.1
sentence-transformers==2.5.1
pydantic==2.6.4
pytest==8.3.5
```

すべての依存関係は`requirements.txt`ファイルに記載されています。

### インストール

```bash
# リポジトリをクローン
git clone https://github.com/yourusername/openMM-Doc-MCP.git
cd openMM-Doc-MCP

# 仮想環境の作成と有効化（uv使用）
uv venv

# uvでパッケージをインストール (推奨)
uv pip install -r requirements.txt

# または、pipを使用する場合
# python -m venv .venv
# source .venv/bin/activate  # Linuxの場合
# .venv\Scripts\activate     # Windowsの場合
# pip install -r requirements.txt
```

### 環境設定

必要に応じて環境変数を設定できます:

```bash
# サーバーポートの設定（デフォルトは8080）
export MCP_SERVER_PORT=8888

# インデックスディレクトリの設定（オプション）
export MCP_INDEX_DIR="/path/to/custom/index"
```

Windowsでの環境変数設定:

```
set MCP_SERVER_PORT=8888
set MCP_INDEX_DIR=C:\path\to\custom\index
```

### インデックスの作成

ドキュメントのインデックスを作成するには、以下のコマンドを実行します：

```bash
uv run python create_faiss_index.py
```

オプションパラメータ:

```bash
# カスタムドキュメントディレクトリを指定
uv run python create_faiss_index.py --docs_dir "/path/to/docs"

# カスタム出力ディレクトリを指定
uv run python create_faiss_index.py --output_dir "/path/to/output"

# カスタム埋め込みモデルを指定
uv run python create_faiss_index.py --embedding_model "intfloat/multilingual-e5-large"
```

インデックス作成には時間がかかる場合があります。処理中はシステムリソース（特にメモリ）を多く使用します。デフォルトでは、`data/indices/docs/`ディレクトリにインデックスファイルが作成されます。

## 使用方法

### サーバーの起動

```bash
uv run python server.py
```

デフォルトでは、サーバーは http://localhost:8080 でリッスンします。

### クエリの実行

コマンドラインから直接検索するには：

```bash
uv run python search_molecular_simulation.py "分子動力学の基本原理"
```

または、HTTPリクエストを送信することもできます：

```
POST http://localhost:8080/query
Content-Type: application/json

{
  "query": "分子動力学の基本原理",
  "top_k": 5
}
```

## テスト

このプロジェクトには、ベクトルデータベースとMCPサーバーのテストが含まれています。

### テストの実行

すべてのテストを実行するには：

```bash
uv run -m pytest
```

特定のテストを実行するには：

```bash
# サーバーのテストのみを実行
uv run -m pytest tests/test_server.py

# ベクトルデータベースのテストを実行
uv run -m pytest src/vector_db/tests/
```

詳細な出力を確認するには `-v` オプションを追加します：

```bash
uv run -m pytest tests/test_server.py -v
```

### テストの構造

- `tests/test_server.py`: MCPサーバーの機能テスト
- `src/vector_db/tests/`: ベクトルデータベース関連のテスト
  - `test_indexer.py`: インデクサーのテスト
  - `test_retriever.py`: レトリーバーのテスト

## API仕様

### MCP連携

このサーバーはModel Context Protocol (MCP)を実装し、以下のツールを提供しています：

1. **search_documents**
   - クエリ文字列に基づいて類似ドキュメントを検索します
   - パラメータ:
     - `query`: 検索クエリテキスト（必須）
     - `top_k`: 返される結果の数（デフォルト5）
     - `index_path`: FAISSインデックスへのカスタムパス（オプション）
   - 戻り値: 関連するドキュメントセクションを含む検索結果の辞書

2. **get_index_info**
   - 現在ロードされているベクトルデータベースインデックスの情報を取得します
   - パラメータ:
     - `index_path`: FAISSインデックスへのカスタムパス（オプション）
   - 戻り値: インデックス情報を含む辞書

詳細なAPI仕様は、以下のファイルで提供されています：
- 日本語: [specs/apispec_ja.md](specs/apispec_ja.md)
- 英語: [specs/apispec_en.md](specs/apispec_en.md)

## Claude Desktop連携

### Claude Desktopでの設定方法

このMCPサーバーをClaude Desktopと連携して、ドキュメント検索機能を強化できます：

1. **Claude Desktopのインストール**: [Anthropicのウェブサイト](https://www.anthropic.com/claude)からダウンロードしてインストールします

2. **MCPサーバーの起動**:
   ```bash
   uv run python server.py
   ```

3. **Claude Desktopの設定**:
   - Claude Desktopを開きます
   - 設定（歯車アイコン）> 詳細設定 > MCPに移動します
   - MCP連携を有効にします
   - 以下の詳細でサーバーを追加します：
     - 名前: OpenMMドキュメンテーション
     - URL: http://localhost:8080
   - 「追加」と「保存」をクリックします

4. **Claudeでの使用方法**:
   - Claude Desktopで新しい会話を開始します
   - ClaudeがOpenMMドキュメントにアクセスできるようになります
   - 「OpenMMのフォースフィールドパラメータの使い方を説明して」や「OpenMMで分子動力学シミュレーションをどのようにセットアップするか」といった質問ができます

### Claude Desktop接続のトラブルシューティング

Claude DesktopがMCPサーバーに接続できない場合：

1. サーバーが実行中であることを確認します（`uv run python server.py`）
2. URLがClaude Desktop設定に正しく入力されていることを確認します
3. ファイアウォールが接続をブロックしていないことを確認します
4. Claude Desktopの再起動を試みてください

## トラブルシューティング

### 一般的な問題

1. **サーバーが起動しない**
   - 依存パッケージが正しくインストールされているか確認してください
   - ポートが既に使用されていないか確認してください
   - ログを確認して具体的なエラーメッセージを確認してください

2. **検索結果が返ってこない**
   - インデックスファイルが正しく作成されているか確認してください
   - クエリが空でないことを確認してください
   - インデックスパスが正しいか確認してください

3. **メモリエラーが発生する**
   - インデックス作成時には十分なメモリ（16GB以上推奨）が必要です
   - 大きなモデルを使用している場合は、より小さいモデルに切り替えてみてください

### デバッグ方法

詳細なログを有効にするには、環境変数を設定します：

```bash
export DEBUG=true
uv run python server.py
```

### テスト実行時の問題

1. **NumPyエラー**
   - `libgfortran.5.dylib`のエラーが発生する場合、`uv run`でテストを実行することで解決できます
   - 環境によっては、`conda install -c conda-forge libgfortran`でライブラリをインストールする必要があります

2. **テストがタイムアウトする**
   - テストタイムアウトを増やす: `uv run -m pytest --timeout=30`

## パフォーマンスチューニング

### レイテンシ最適化

- より小さいインデックスを使用する
- `faiss-gpu`を使用してGPU高速化を有効にする（対応GPUを搭載した環境の場合）
- インデクシングパラメータを調整する

### メモリ使用量の最適化

大量のドキュメントを処理する場合、メモリ使用量を最適化するために以下の設定を変更できます：

```python
# create_faiss_index.pyでのチャンクサイズの調整
--chunk_size 256  # デフォルトは512
--chunk_overlap 20  # デフォルトは50
```

## 開発者向け情報

### コーディング規約

このプロジェクトでは以下のコーディング規約に従ってください：

- PEP 8スタイルガイド
- 関数やクラスには適切なdocstringを記述する
- テストカバレッジを維持する

### 貢献方法

1. リポジトリをフォークする
2. 機能ブランチを作成する (`git checkout -b feature/amazing-feature`)
3. 変更をコミットする (`git commit -m 'Add some amazing feature'`)
4. ブランチにプッシュする (`git push origin feature/amazing-feature`)
5. プルリクエストを作成する

### 拡張ガイド

#### 新しいベクトルストアの追加

現在はFAISSを使用していますが、他のベクトルストアを追加することも可能です：

1. `src/vector_db/indexer.py`に新しいインデクサークラスを追加
2. `src/vector_db/retriever.py`に対応するレトリーバークラスを追加
3. 新しいインデクサー/レトリーバーのテストを作成
4. 環境変数またはコマンドライン引数で選択できるようにする

## ディレクトリ構造

```
openMM-Doc-MCP/
├── create_faiss_index.py   # インデックス作成スクリプト
├── search_molecular_simulation.py # CLIからの検索
├── server.py               # MCPサーバー実装
├── data/
│   └── indices/
│       └── docs/
│           ├── index.faiss # FAISSインデックスファイル
│           └── index.pkl   # Pickleインデックスファイル
├── specs/
│   ├── apispec_en.md       # API仕様書（英語）
│   └── apispec_ja.md       # API仕様書（日本語）
├── src/
│   └── vector_db/          # ベクトルデータベース関連モジュール
│       ├── indexer.py      # インデクサー実装
│       ├── retriever.py    # レトリーバー実装
│       └── tests/          # ベクトルDB関連のテスト
│           ├── conftest.py
│           ├── test_indexer.py
│           └── test_retriever.py
└── tests/
    └── test_server.py      # サーバーのテスト
```

## ライセンス

このプロジェクトは[LICENSE](LICENSE)の下で提供されています。