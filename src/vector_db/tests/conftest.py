"""
pytestのテスト用の共通設定やフィクスチャを提供するモジュール
"""
import os
import tempfile
import pytest
import numpy as np
from unittest.mock import MagicMock
from langchain_community.embeddings import HuggingFaceEmbeddings
from src.vector_db.indexer import Embedder

class MockEmbeddings:
    """テスト用のモック埋め込みモデル"""
    def __init__(self):
        pass
    
    def embed_documents(self, texts):
        """文書を埋め込みベクトルに変換する（モック版）"""
        # 文書ごとに固定サイズの一意なベクトルを返す
        return [np.random.rand(384).tolist() for _ in texts]
    
    def embed_query(self, text):
        """クエリを埋め込みベクトルに変換する（モック版）"""
        # 固定サイズのランダムベクトルを返す
        return np.random.rand(384).tolist()

@pytest.fixture(scope="session")
def temp_dir():
    """テスト用の一時ディレクトリを提供するフィクスチャ"""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir

@pytest.fixture(scope="session")
def test_markdown_content():
    """テスト用のマークダウンコンテンツを提供するフィクスチャ"""
    return """# OpenMM ドキュメント

## 概要
OpenMMは分子シミュレーションのための高性能ツールキットです。

## 使い方
まず、システムをセットアップします。
次に、シミュレーションを実行します。

## 高度な機能
OpenMMはカスタム力場もサポートしています。
"""

@pytest.fixture(scope="session")
def test_markdown_file(temp_dir, test_markdown_content):
    """テスト用のマークダウンファイルを作成して提供するフィクスチャ"""
    file_path = os.path.join(temp_dir, "test_document.md")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(test_markdown_content)
    return file_path

@pytest.fixture(scope="session")
def mock_embeddings():
    """テスト用のモック埋め込みモデルを提供するフィクスチャ"""
    return MockEmbeddings()

@pytest.fixture(scope="session")
def test_embedder(mock_embeddings):
    """テスト用の埋め込みモデルを提供するフィクスチャ"""
    # 実際のEmbedderクラスを使用するが、内部の埋め込みモデルをモックに置き換える
    embedder = Embedder(model_name="intfloat/multilingual-e5-small")
    embedder.embeddings = mock_embeddings
    return embedder

@pytest.fixture(scope="session")
def test_index_path(temp_dir):
    """テスト用のインデックスパスを提供するフィクスチャ"""
    return os.path.join(temp_dir, "index.faiss")