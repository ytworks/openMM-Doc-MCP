# openMM-Doc-MCP

OpenMMドキュメントをベクターDBで検索できるようにするModel Context Protocol (MCP) サーバーです。

## セットアップ

### 依存パッケージのインストール

このプロジェクトは依存関係の管理に`uv`を使用しています。以下のコマンドで必要なパッケージをインストールできます。

```bash
uv pip install -r requirements.txt
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

## プロジェクト構成

```
.
├── create_faiss_index.py   # ベクターDBインデックス作成スクリプト
├── data
│   └── indices             # ベクターDBの保存先
│       └── docs            # 生成されたFAISSインデックス
├── specs                   # マークダウンドキュメント
│   ├── apispec_en.md       # 英語APIドキュメント
│   └── apispec_ja.md       # 日本語APIドキュメント
└── src
    └── vector_db           # ベクターDB関連のコード
        ├── indexer.py      # インデックス作成モジュール
        ├── retriever.py    # 検索モジュール
        └── tests           # テストコード
            ├── conftest.py # テスト用フィクスチャ
            ├── test_indexer.py  # インデクサーのテスト
            └── test_retriever.py # リトリーバーのテスト