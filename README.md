# openMM-Doc-MCP

OpenMMドキュメントをベクターDBで検索できるようにするModel Context Protocol (MCP) サーバーです。ベクターデータベースを使用して分子シミュレーションドキュメントの類似検索を行うことができます。

## セットアップ

### 依存パッケージのインストール

このプロジェクトは依存関係の管理に`uv`を使用しています。以下のコマンドで必要なパッケージをインストールできます。

```bash
uv pip install -r requirements.txt
```

MCPを使用するには、以下のコマンドで追加パッケージをインストールします。

```bash
uv add "mcp>=1.2.0"
```

## ベクターDBの作成

マークダウンドキュメントからベクターDBを作成するには、以下のコマンドを実行します。

```bash
uv run python create_faiss_index.py
```

このコマンドにより、`specs`ディレクトリ内のマークダウンファイルが処理され、`data/indices/docs`ディレクトリにFAISSインデックスが生成されます。

### 主なオプション

`create_faiss_index.py`スクリプトには以下のオプションがあります:

- `--docs-path`: マークダウンファイルが格納されているディレクトリパス（デフォルト: `specs`）
- `--output-dir`: 出力ディレクトリ（デフォルト: `data/indices`）
- `--index-name`: インデックスの名前（デフォルト: `docs`）
- `--chunk-size`: チャンクサイズ（文字数）（デフォルト: `100`）
- `--chunk-overlap`: チャンクオーバーラップ（文字数）（デフォルト: `10`）
- `--embedding-model`: 埋め込みモデル（デフォルト: `intfloat/multilingual-e5-large`）

例えば、チャンクサイズを変更する場合は以下のように実行します。

```bash
uv run python create_faiss_index.py --chunk-size 200 --chunk-overlap 20
```

## MCPサーバーの実行

MCPサーバーを起動するには、以下のコマンドを実行します。

```bash
uv run --with mcp[cli] mcp run server.py
```

これにより、分子シミュレーションドキュメントの検索機能を提供するMCPサーバーが起動します。

## MCPサーバーの利用方法

MCPサーバーは以下の機能を提供します：

### search_documents

クエリテキストに基づいて類似ドキュメントを検索します。

**パラメータ**:
- `query`: 検索クエリテキスト（必須）
- `top_k`: 返される結果の数（デフォルト: 5）
- `index_path`: カスタムインデックスパス（オプション）

**使用例**:
```
分子シミュレーションのパラメータ設定について教えてください
```

### get_index_info

ベクターデータベースインデックスの情報を取得します。

**パラメータ**:
- `index_path`: カスタムインデックスパス（オプション）

**使用例**:
```
インデックスの情報を表示してください
```

## ドキュメント検索スクリプト

コマンドラインでドキュメント検索を行うには、以下のスクリプトを使用できます。

```bash
uv run python search_molecular_simulation.py --query "分子シミュレーションのパラメータ"
```

## テストの実行

プロジェクト内のテストは`pytest`を使用して実行できます。以下のコマンドでベクターDBの作成と読み込みに関するテストを実行できます。

### 全テストの実行

```bash
uv run python -m pytest src/vector_db/tests -v
```

### インデクサーのテスト実行

ベクターDB作成処理のテストのみを実行する場合：

```bash
uv run python -m pytest src/vector_db/tests/test_indexer.py -v
```

### リトリーバーのテスト実行

ベクターDB検索処理のテストのみを実行する場合：

```bash
uv run python -m pytest src/vector_db/tests/test_retriever.py -v
```

### 個別のテストケース実行

特定のテストケースのみを実行する場合：

```bash
uv run python -m pytest src/vector_db/tests/test_indexer.py::test_markdown_chunker -v
```

## Claude Desktopとの統合

Claude Desktopで本MCPサーバーを使用するには、Claude Desktopの設定ファイルを編集します。

設定ファイルの場所:
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`

以下の設定を追加してください（絶対パスに置き換えてください）:

```json
{
  "mcpServers": {
    "OpenMM Document Search Service": {
      "command": "uv",
      "args": [
        "run",
        "--with",
        "mcp[cli]",
        "--with",
        "langchain",
        "--with",
        "langchain-community",
        "--with",
        "sentence-transformers",
        "--with",
        "faiss-cpu",
        "mcp",
        "run",
        "/Users/yuzotakagi/dev/drug/openMM-Doc-MCP/server.py"
      ]
    }
  }
}
```

**重要**: 
- 上記の例では絶対パス (`/Users/yuzotakagi/dev/drug/openMM-Doc-MCP/server.py`) を使用していますが、これはあなたの環境に合わせて変更してください。
- `sentence-transformers`と`faiss-cpu`が含まれていることを確認してください。これらはベクトル埋め込みと検索に必要なパッケージです。
- 初めて実行する際はモデルのダウンロードが行われるため、少し時間がかかる場合があります。

## プロジェクト構成

```
.
├── create_faiss_index.py   # ベクターDBインデックス作成スクリプト
├── search_molecular_simulation.py # コマンドライン検索スクリプト
├── server.py               # MCPサーバー実装
├── data
│   └── indices             # ベクターDBの保存先
│       └── docs            # 生成されたFAISSインデックス
│           ├── index.faiss # FAISSインデックスファイル
│           └── index.pkl   # メタデータファイル
├── specs                   # マークダウンドキュメント
│   ├── apispec_en.md       # 英語APIドキュメント
│   └── apispec_ja.md       # 日本語APIドキュメント
└── src
    └── vector_db           # ベクターDB関連のコード
        ├── __init__.py     # パッケージ初期化ファイル
        ├── indexer.py      # インデックス作成モジュール
        ├── retriever.py    # 検索モジュール
        └── tests           # テストコード
            ├── conftest.py # テスト用フィクスチャ
            ├── test_indexer.py  # インデクサーのテスト
            └── test_retriever.py # リトリーバーのテスト