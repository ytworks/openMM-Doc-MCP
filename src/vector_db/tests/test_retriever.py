"""
リトリーバー機能のテストモジュール
"""
import os
import sys
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from src.vector_db.indexer import MarkdownChunker, Embedder, FAISSIndexer
from src.vector_db.retriever import FAISSRetriever


# 外部ライブラリへのアクセスをモックに置き換え
@pytest.fixture
def mock_faiss():
    """FAISSのモックを提供するフィクスチャ"""
    with patch('langchain_community.vectorstores.FAISS.from_texts') as mock_from_texts, \
         patch('langchain_community.vectorstores.FAISS.load_local') as mock_load_local, \
         patch('langchain_community.vectorstores.FAISS.similarity_search') as mock_similarity_search, \
         patch('langchain_community.vectorstores.FAISS.similarity_search_by_vector') as mock_similarity_search_by_vector:
        
        # モックDBとその動作を設定
        mock_db = MagicMock()
        mock_db.save_local = MagicMock()
        mock_from_texts.return_value = mock_db
        mock_load_local.return_value = mock_db
        
        # モック検索結果を設定
        mock_doc1 = MagicMock()
        mock_doc1.page_content = "OpenMMは分子シミュレーションのためのツールキットです。"
        mock_doc1.metadata = {'source': 'test.md', 'id': 'test_1'}
        
        mock_doc2 = MagicMock()
        mock_doc2.page_content = "シミュレーションを実行するにはSystemを設定します。"
        mock_doc2.metadata = {'source': 'test.md', 'id': 'test_2'}
        
        mock_similarity_search.return_value = [mock_doc1, mock_doc2]
        mock_similarity_search_by_vector.return_value = [mock_doc1]
        
        yield {
            'mock_db': mock_db,
            'mock_from_texts': mock_from_texts,
            'mock_load_local': mock_load_local,
            'mock_similarity_search': mock_similarity_search,
            'mock_similarity_search_by_vector': mock_similarity_search_by_vector
        }


@pytest.fixture
def create_test_index(mock_faiss, test_markdown_file, temp_dir, test_embedder):
    """テスト用のFAISSインデックスを作成するフィクスチャ"""
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
        index_path = os.path.join(temp_dir, "index")
        
        # インデックスを保存
        indexer.save(index_path)
    
    return {
        'index_path': index_path,
        'embedder': test_embedder,
        'chunks': chunks,
        'mock_faiss': mock_faiss
    }


def test_retriever_initialization(mock_faiss, create_test_index):
    """FAISSRetrieverが正しく初期化できることをテスト"""
    # テスト用のインデックスパスを取得
    index_path = create_test_index['index_path']
    
    # リトリーバーを初期化
    with patch('src.vector_db.retriever.FAISS') as mock_faiss_class:
        mock_faiss_class.load_local.return_value = mock_faiss['mock_db']
        
        retriever = FAISSRetriever(index_path)
        
        # load_localが呼び出されたことを確認
        assert mock_faiss_class.load_local.call_count > 0
        
        # リトリーバーが正しく初期化されたことを確認
        assert retriever.db is not None
        assert retriever.embeddings is not None


def test_search_by_text(mock_faiss, create_test_index):
    """テキストによる検索機能をテスト"""
    # テスト用のインデックスパスを取得
    index_path = create_test_index['index_path']
    
    # リトリーバーを初期化
    with patch('src.vector_db.retriever.FAISS') as mock_faiss_class:
        mock_faiss_class.load_local.return_value = mock_faiss['mock_db']
        
        retriever = FAISSRetriever(index_path)
        retriever.db.similarity_search = mock_faiss['mock_similarity_search']
        
        # テキスト検索を実行
        results = retriever.search_by_text("OpenMM", top_k=2)
        
        # similarity_searchが呼び出されたことを確認
        assert mock_faiss['mock_similarity_search'].call_count > 0
        
        # 検索結果が正しく返されたことを確認
        assert len(results) == 2
        assert isinstance(results, list)
        
        # 各結果が必要なフィールドを含んでいることを確認
        for result in results:
            assert 'content' in result
            assert 'metadata' in result
            assert 'source' in result
            assert 'id' in result


def test_search_by_vector(mock_faiss, create_test_index):
    """ベクトルによる検索機能をテスト"""
    # テスト用のインデックスパスと埋め込みモデルを取得
    index_path = create_test_index['index_path']
    embedder = create_test_index['embedder']
    
    # リトリーバーを初期化
    with patch('src.vector_db.retriever.FAISS') as mock_faiss_class:
        mock_faiss_class.load_local.return_value = mock_faiss['mock_db']
        
        retriever = FAISSRetriever(index_path, embeddings=embedder.embeddings)
        retriever.db.similarity_search_by_vector = mock_faiss['mock_similarity_search_by_vector']
        
        # クエリテキストをベクトルに変換
        query_text = "分子シミュレーション"
        query_vector = embedder.embeddings.embed_query(query_text)
        
        # ベクトル検索を実行
        results = retriever.search_by_vector(query_vector, top_k=2)
        
        # similarity_search_by_vectorが呼び出されたことを確認
        assert mock_faiss['mock_similarity_search_by_vector'].call_count > 0
        
        # 検索結果が正しく返されたことを確認
        assert len(results) > 0
        assert isinstance(results, list)
        
        # 各結果が必要なフィールドを含んでいることを確認
        for result in results:
            assert 'content' in result
            assert 'metadata' in result
            assert 'source' in result
            assert 'id' in result


def test_integration_indexer_and_retriever(mock_faiss, test_markdown_file, temp_dir, test_embedder):
    """インデクサーとリトリーバーの連携をテスト"""
    # マークダウンをチャンクに分割
    chunker = MarkdownChunker(chunk_size=200, chunk_overlap=20)
    chunks = chunker.process_markdown(test_markdown_file)
    
    # チャンクを埋め込みベクトルに変換
    vectors, meta = test_embedder.embed_chunks(chunks)
    
    # FAISSインデックスを作成して保存
    with patch.object(FAISSIndexer, 'add_embeddings') as mock_add_embeddings, \
         patch('src.vector_db.retriever.FAISS') as mock_faiss_class:
        
        indexer = FAISSIndexer(embedding_model=test_embedder)
        indexer.db = mock_faiss['mock_db']
        indexer.add_embeddings(vectors, meta)
        
        # インデックスを保存するパスを設定
        index_path = os.path.join(temp_dir, "test_index.faiss")
        
        # インデックスを保存
        indexer.save(index_path)
        
        # save_localが呼び出されたことを確認
        assert indexer.db.save_local.call_count > 0
        
        # リトリーバー用のモック設定
        mock_faiss_class.load_local.return_value = mock_faiss['mock_db']
        
        # リトリーバーを初期化して検索を実行
        retriever = FAISSRetriever(index_path)
        retriever.db.similarity_search = mock_faiss['mock_similarity_search']
        
        results = retriever.search_by_text("OpenMM", top_k=2)
        
        # similarity_searchが呼び出されたことを確認
        assert mock_faiss['mock_similarity_search'].call_count > 0
        
        # 検索結果が正しく返されたことを確認
        assert len(results) > 0
        
        # 検索結果が期待通りであることを確認
        assert "OpenMM" in results[0]['content']