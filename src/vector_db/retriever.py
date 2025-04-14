"""
FAISSインデックスからの検索機能を提供するモジュール
"""
from typing import List, Dict, Any, Optional

from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document


class FAISSRetriever:
    """FAISSインデックスを使用して類似テキストを検索するクラス"""
    
    def __init__(self, index_path: str, embeddings=None):
        """
        Args:
            index_path: インデックスファイルのパス
            embeddings: 埋め込みモデル（指定がない場合は新しく作成）
        """
        # 埋め込みモデルの設定
        if embeddings is None:
            self.embeddings = HuggingFaceEmbeddings(
                model_name="intfloat/multilingual-e5-large"
            )
        else:
            self.embeddings = embeddings
        
        # インデックスの読み込み
        self.db = FAISS.load_local(index_path, self.embeddings, allow_dangerous_deserialization=True)
    
    def search_by_text(self, query_text: str, top_k: int = 5) -> List[Dict[str, Any]]:
        """テキストを使用して類似テキストを検索
        
        Args:
            query_text: 検索クエリのテキスト
            top_k: 取得する結果の最大数
            
        Returns:
            類似テキストとメタデータのリスト
        """
        # LangChainのsimilarity_searchメソッドを使用
        docs = self.db.similarity_search(query_text, k=top_k)
        
        # 結果をまとめる
        results = []
        for doc in docs:
            # Document オブジェクトから情報を取得
            results.append({
                'content': doc.page_content,
                'metadata': doc.metadata,
                'source': doc.metadata.get('source', ''),
                'id': doc.metadata.get('id', '')
            })
        
        return results
    
    def search_by_vector(self, query_vector: List[float], top_k: int = 5) -> List[Dict[str, Any]]:
        """ベクトルを使用して類似テキストを検索
        
        Args:
            query_vector: 検索クエリのベクトル
            top_k: 取得する結果の最大数
            
        Returns:
            類似テキストとメタデータのリスト
        """
        # LangChainのsimilarity_search_by_vectorメソッドを使用
        docs = self.db.similarity_search_by_vector(query_vector, k=top_k)
        
        # 結果をまとめる
        results = []
        for doc in docs:
            results.append({
                'content': doc.page_content,
                'metadata': doc.metadata,
                'source': doc.metadata.get('source', ''),
                'id': doc.metadata.get('id', '')
            })
        
        return results
