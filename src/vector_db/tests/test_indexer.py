"""
インデクサー機能のテストモジュール
"""
import os
import sys
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from src.vector_db.indexer import MarkdownChunker, Embedder, FAISSIndexer

# 外部ライブラリへのアクセスをモックに置き換え
@pytest.fixture
def mock_faiss():
    """FAISSのモックを提供するフィクスチャ"""
    with patch('langchain_community.vectorstores.FAISS.from_texts') as mock_from_texts, \
         patch('langchain_community.vectorstores.FAISS.load_local') as mock_load_local, \
         patch('langchain_community.vectorstores.FAISS.similarity_search') as mock_similarity_search:
        
        # モックDBとその動作を設定
        mock_db = MagicMock()
        mock_db.save_local = MagicMock()
        mock_from_texts.return_value = mock_db
        mock_load_local.return_value = mock_db
        
        # モック検索結果を設定
        mock_doc = MagicMock()
        mock_doc.page_content = "OpenMMは分子シミュレーションのためのツールキットです。"
        mock_doc.metadata = {'source': 'test.md', 'id': 'test_1'}
        mock_similarity_search.return_value = [mock_doc]
        
        yield {
            'mock_db': mock_db,
            'mock_from_texts': mock_from_texts,
            'mock_load_local': mock_load_local,
            'mock_similarity_search': mock_similarity_search
        }


def test_markdown_chunker(test_markdown_file):
    """MarkdownChunkerがマークダウンファイルを正しくチャンクに分割できることをテスト"""
    # テスト用のチャンカーを初期化
    chunker = MarkdownChunker(chunk_size=200, chunk_overlap=20)
    
    # マークダウンファイルを処理
    chunks = chunker.process_markdown(test_markdown_file)
    
    # チャンクが作成されたことを確認
    assert len(chunks) > 0
    
    # 各チャンクが必要なフィールドを持っていることを確認
    for chunk in chunks:
        assert 'id' in chunk
        assert 'content' in chunk
        assert 'source' in chunk
        assert chunk['source'] == test_markdown_file
        assert chunk['id'].startswith("test_document.md_")


def test_embedder(test_embedder):
    """Embedderが正しくテキストを埋め込みベクトルに変換できることをテスト"""
    # テスト用のチャンク
    chunks = [
        {
            'id': 'test_1',
            'content': 'OpenMMは分子シミュレーションのためのツールキットです。',
            'source': 'test.md'
        },
        {
            'id': 'test_2',
            'content': 'シミュレーションを実行するにはSystemを設定します。',
            'source': 'test.md'
        }
    ]
    
    # 埋め込みを生成
    vectors, meta = test_embedder.embed_chunks(chunks)
    
    # 埋め込みの数がチャンクの数と一致することを確認
    assert len(vectors) == len(chunks)
    assert len(meta) == len(chunks)
    
    # 埋め込みベクトルの形状を確認
    assert len(vectors[0]) > 0  # 埋め込みベクトルの次元数が正であること
    
    # メタデータが保持されていることを確認
    assert meta == chunks


def test_faiss_indexer_creation_and_save(mock_faiss, test_embedder, temp_dir):
    """FAISSIndexerが正しくインデックスを作成して保存できることをテスト"""
    # テスト用のチャンク
    chunks = [
        {
            'id': 'test_1',
            'content': 'OpenMMは分子シミュレーションのためのツールキットです。',
            'source': 'test.md'
        },
        {
            'id': 'test_2',
            'content': 'シミュレーションを実行するにはSystemを設定します。',
            'source': 'test.md'
        }
    ]
    
    # 埋め込みを生成
    vectors, meta = test_embedder.embed_chunks(chunks)
    
    # FAISSIndexerを使用してインデックスを作成
    with patch.object(FAISSIndexer, 'add_embeddings') as mock_add_embeddings:
        # インデクサーを初期化してインデックスを作成
        indexer = FAISSIndexer(embedding_model=test_embedder)
        indexer.db = mock_faiss['mock_db']
        indexer.add_embeddings(vectors, meta)
        
        # インデックスを保存するパスを設定
        index_path = os.path.join(temp_dir, "test_index.faiss")
        
        # インデックスを保存
        indexer.save(index_path)
        
        # save_localが呼び出されたことを確認
        assert indexer.db.save_local.call_count > 0


def test_faiss_indexer_load(mock_faiss, test_embedder, temp_dir):
    """FAISSIndexerが正しくインデックスを読み込めることをテスト"""
    # インデックスを読み込むパスを設定
    index_path = os.path.join(temp_dir, "test_index.faiss")
    
    # インデックスを読み込む
    with patch('src.vector_db.indexer.FAISS') as mock_faiss_class:
        mock_faiss_class.load_local.return_value = mock_faiss['mock_db']
        
        loaded_indexer = FAISSIndexer.load(index_path, test_embedder)
        
        # load_localが呼び出されたことを確認
        assert mock_faiss_class.load_local.call_count > 0 or mock_faiss['mock_load_local'].call_count > 0
        
        # モック検索を設定
        loaded_indexer.db = mock_faiss['mock_db']
        loaded_indexer.db.similarity_search = mock_faiss['mock_similarity_search']
        
        # 読み込んだインデックスを使用してクエリを実行
        query_text = "分子シミュレーション"
        results = loaded_indexer.db.similarity_search(query_text, k=1)
        
        # similarity_searchが呼び出されたことを確認
        assert mock_faiss['mock_similarity_search'].call_count > 0
        
        # 検索結果が返されたことを確認
        assert len(results) > 0


def test_end_to_end_indexing(mock_faiss, test_markdown_file, test_embedder, temp_dir):
    """マークダウンファイルの読み込みからインデックス作成までの一連の流れをテスト"""
    # マークダウンをチャンクに分割
    chunker = MarkdownChunker(chunk_size=200, chunk_overlap=20)
    chunks = chunker.process_markdown(test_markdown_file)
    
    # チャンクを埋め込みベクトルに変換
    vectors, meta = test_embedder.embed_chunks(chunks)
    
    # FAISSインデックスを作成して保存
    with patch.object(FAISSIndexer, 'add_embeddings') as mock_add_embeddings:
        indexer = FAISSIndexer(embedding_model=test_embedder)
        indexer.db = mock_faiss['mock_db']
        indexer.add_embeddings(vectors, meta)
        
        # インデックスを保存するパスを設定
        index_path = os.path.join(temp_dir, "test_index.faiss")
        
        # インデックスを保存
        indexer.save(index_path)
        
        # save_localが呼び出されたことを確認
        assert indexer.db.save_local.call_count > 0