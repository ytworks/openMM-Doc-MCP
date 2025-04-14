# OpenMM Documentation MCP Server

OpenMMの分子動力学シミュレーションドキュメントのための検索サーバーです。Model Context Protocol (MCP)に準拠しており、LLM（大規模言語モデル）との連携に最適化されています。

> 🌐 **Language/言語**: [English](README.md) | [日本語](README_ja.md)

## 概要

このMCPサーバーは、OpenMMのドキュメントを自然言語で検索する機能を提供します。最新の言語モデルを使用してドキュメントの内容をベクトル埋め込みに変換し、FAISSベクトルデータベースに保存して効率的な検索を実現しています。クエリを受け取ると、意味的に最も関連性の高いドキュメントセクションを検索して返します。特に以下のような用途に最適です：

- 分子動力学手法に関連するドキュメントの検索
- OpenMMの関数やクラスの使用方法の検索
- シミュレーションパラメータや設定に関する情報の取得
- コード例や実装の詳細の参照

## 機能

- **意味検索**: キーワードだけでなく、意味に基づいてドキュメントを検索
- **MCP連携**: Claude Desktopやその他のMCP対応アプリケーションと完全互換
- **多言語サポート**: 日本語と英語の両方のクエリに対応
- **効率的な検索**: FAISSを使用した高性能なベクトル類似性検索
- **カスタマイズ可能**: 埋め込みモデルや検索パラメータの設定が可能

## セットアップ

### 前提条件

- Python 3.9以上
- `uv` パッケージマネージャー（推奨）または `pip`
- 最低8GB RAM（インデックス作成時には16GB以上を推奨）

### インストール

```bash
# リポジトリをクローン
git clone https://github.com/yourusername/openMM-Doc-MCP.git
cd openMM-Doc-MCP

# uvを使用して仮想環境を作成して有効化
uv venv

# uvでパッケージをインストール（推奨）
uv pip install -r requirements.txt

# pipを使用する場合
# python -m venv .venv
# source .venv/bin/activate  # Linux/macOSの場合
# .venv\Scripts\activate     # Windowsの場合
# pip install -r requirements.txt
```

### インデックスの作成

サーバーを使用する前に、OpenMMのドキュメントのベクトルインデックスを作成する必要があります：

```bash
uv run python create_faiss_index.py
```

オプションパラメータ：

```bash
# カスタムドキュメントディレクトリを指定
uv run python create_faiss_index.py --docs_dir "/path/to/docs"

# カスタム出力ディレクトリを指定
uv run python create_faiss_index.py --output_dir "/path/to/output"

# カスタム埋め込みモデルを指定
uv run python create_faiss_index.py --embedding_model "intfloat/multilingual-e5-large"
```

インデックス作成にはある程度の時間とメモリが必要です。デフォルトでは、インデックスファイルは `data/indices/docs/` ディレクトリに作成されます。

### 設定

環境変数を使用してサーバーを設定できます：

```bash
# サーバーポートを設定（デフォルトは8080）
export MCP_SERVER_PORT=8888

# インデックスディレクトリを設定（オプション）
export MCP_INDEX_DIR="/path/to/custom/index"
```

Windowsの場合：

```
set MCP_SERVER_PORT=8888
set MCP_INDEX_DIR=C:\path\to\custom\index
```

## 使用方法

### サーバーの起動

```bash
uv run python server.py
```

デフォルトでは、サーバーは http://localhost:8080 でリッスンします。

### コマンドライン検索

コマンドラインから直接検索することができます：

```bash
uv run python search_molecular_simulation.py "水箱シミュレーションの設定方法"
```

### HTTPリクエスト

サーバーにHTTPリクエストを送信することもできます：

```
POST http://localhost:8080/query
Content-Type: application/json

{
  "query": "水箱シミュレーションの設定方法",
  "top_k": 5
}
```

## Claude Desktop連携

### Claude Desktopでの設定方法

Claude Desktopの設定ファイルを編集して、このMCPサーバーを追加します。設定ファイルのパスは：

- **macOS**:
  `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**:
  `%APPDATA%\Claude\claude_desktop_config.json`

以下のJSON設定を（既存の`mcpServers`オブジェクト内に）追加します：

```json
{
  "mcpServers": {
    "OpenMMドキュメンテーション": {
      "command": "uv",
      "args": [
        "run",
        "--with",
        "mcp[cli]",
        "--with",
        "faiss-cpu",
        "--with",
        "langchain",
        "--with",
        "sentence-transformers",
        "mcp",
        "run",
        "/path/to/openMM-Doc-MCP/server.py"
      ]
    }
  }
}
```

注意点：
- `uv`コマンドが環境パスに含まれていない場合は、絶対パスを使用してください（例：`/path/to/uv`）。
- `/path/to/openMM-Doc-MCP/server.py`をこのスクリプトの絶対パスに置き換えてください。
- 常に相対パスではなく絶対パスを使用してください。

### Claude Desktop接続のトラブルシューティング

Claude DesktopがMCPサーバーに接続できない場合：

1. 設定ファイル内の`server.py`へのパスが正しいことを確認してください（絶対パス）
2. `uv`が正しくインストールされてアクセス可能であることを確認してください
3. システムログでエラーを確認してください
4. 設定変更後にClaude Desktopを再起動してみてください

## API リファレンス

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
- [specs/apispec_en.md](specs/apispec_en.md)
- [specs/apispec_ja.md](specs/apispec_ja.md)

## テスト

このプロジェクトには、ベクトルデータベースとMCPサーバーのテストが含まれています。

### テストの実行

```bash
# すべてのテストを実行
uv run -m pytest

# 特定のテストを実行
uv run -m pytest tests/test_server.py
uv run -m pytest src/vector_db/tests/

# 詳細な出力
uv run -m pytest tests/test_server.py -v
```

## トラブルシューティング

### 一般的な問題

1. **サーバーが起動しない**
   - 依存パッケージが正しくインストールされているか確認してください
   - ポートが既に使用されていないか確認してください
   - ログで具体的なエラーメッセージを確認してください

2. **検索結果が返ってこない**
   - インデックスファイルが正しく作成されているか確認してください
   - クエリが空でないことを確認してください
   - インデックスパスが正しいか確認してください

3. **メモリエラーが発生する**
   - インデックス作成時には十分なメモリ（16GB以上推奨）が必要です
   - 大きなモデルを使用している場合は、より小さいモデルに切り替えてみてください

### デバッグ

詳細なログを有効にするには：

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

- より小さい埋め込みモデルを使用して推論を高速化
- `faiss-gpu`を使用してGPU高速化を有効にする（対応GPUを搭載した環境の場合）
- インデックス作成時のチャンクサイズとオーバーラップパラメータを調整

### メモリ使用量の最適化

大量のドキュメントを処理する場合：

```bash
# create_faiss_index.pyでチャンクサイズを調整
uv run python create_faiss_index.py --chunk_size 256 --chunk_overlap 20
```

## ディレクトリ構造

```
openMM-Doc-MCP/
├── create_faiss_index.py   # インデックス作成スクリプト
├── search_molecular_simulation.py # コマンドライン検索ユーティリティ
├── server.py               # MCPサーバー実装
├── data/
│   └── indices/
│       └── docs/           # インデックスファイルのデフォルト場所
│           ├── index.faiss # FAISSインデックスファイル
│           └── index.pkl   # メタデータPickleファイル
├── specs/
│   ├── apispec_en.md       # API仕様書（英語）
│   └── apispec_ja.md       # API仕様書（日本語）
├── src/
│   └── vector_db/          # ベクトルデータベース関連モジュール
│       ├── indexer.py      # インデクサー実装
│       ├── retriever.py    # レトリーバー実装
│       └── tests/          # ベクトルDBテスト
└── tests/
    └── test_server.py      # サーバーテスト
```

## ライセンス

このプロジェクトは[LICENSE](LICENSE)の下で提供されています。