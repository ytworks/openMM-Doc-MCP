"""
マークダウンドキュメントの解析とFAISSインデックスの作成を行うモジュール
"""
import os
from typing import List, Dict, Tuple, Any

from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import FAISS


class MarkdownChunker:
    """マークダウンファイルをチャンクに分割するクラス"""
    
    def __init__(self, chunk_size: int = 500, chunk_overlap: int = 50):
        """
        Args:
            chunk_size: 各チャンクの最大文字数
            chunk_overlap: 連続するチャンク間のオーバーラップ文字数
        """
        self.chunk_size = chunk_size
        self.chunk_overlap = chunk_overlap
        self.text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size,
            chunk_overlap=chunk_overlap,
            separators=["\n## ", "\n### ", "\n#### ", "\n##### ", "\n###### ", "\n", " ", ""]
        )
    
    def read_markdown_file(self, file_path: str) -> str:
        """マークダウンファイルを読み込む
        
        Args:
            file_path: マークダウンファイルのパス
            
        Returns:
            ファイルの内容
        """
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    
    def process_markdown(self, file_path: str) -> List[Dict[str, str]]:
        """マークダウンファイルを処理してチャンクのリストを返す
        
        Args:
            file_path: マークダウンファイルのパス
            
        Returns:
            チャンクのリスト
        """
        content = self.read_markdown_file(file_path)
        
        # LangChainのテキストスプリッターでチャンキング
        texts = self.text_splitter.split_text(content)
        
        # メタデータとIDを追加
        chunks = []
        for i, text in enumerate(texts):
            chunk = {
                'id': f"{os.path.basename(file_path)}_{i}",
                'content': text,
                'source': file_path
            }
            chunks.append(chunk)
        
        return chunks


class Embedder:
    """テキストをベクトル化するクラス"""
    
    def __init__(self, model_name: str = "intfloat/multilingual-e5-large"):
        """
        Args:
            model_name: 使用する埋め込みモデルの名前
        """
        self.embeddings = HuggingFaceEmbeddings(
            model_name=model_name
        )
    
    def embed_chunks(self, chunks: List[Dict[str, str]]) -> Tuple[List[List[float]], List[Dict[str, str]]]:
        """チャンクリストを埋め込みベクトルに変換する
        
        Args:
            chunks: チャンクのリスト
            
        Returns:
            埋め込みベクトルの配列とチャンク情報
        """
        texts = [chunk['content'] for chunk in chunks]
        embeddings = self.embeddings.embed_documents(texts)
        return embeddings, chunks


class FAISSIndexer:
    """FAISS インデックスを作成・保存するクラス"""
    
    def __init__(self, embedding_model: Any = None):
        """
        Args:
            embedding_model: 埋め込みに使用するモデル
        """
        self.embedding_model = embedding_model
        self.db = None
    
    def add_embeddings(self, embeddings: List[List[float]], metadata: List[Dict[str, str]]) -> None:
        """埋め込みベクトルとメタデータをインデックスに追加する
        
        Args:
            embeddings: 埋め込みベクトルの配列
            metadata: 対応するメタデータのリスト
        """
        # 直接FAISS.from_embeddingsを使うわけではないので、
        # contentsとembeddingsを使って独自にDBを作成
        texts = [item['content'] for item in metadata]
        self.db = FAISS.from_texts(texts, self.embedding_model.embeddings, metadatas=metadata)
    
    def save(self, index_path: str, metadata_path: str = None) -> None:
        """インデックスとメタデータを保存する
        
        Args:
            index_path: インデックスファイルのパス
            metadata_path: 未使用（LangChainのFAISSはインデックスとメタデータを一緒に保存）
        """
        # ディレクトリ名を取得
        index_dir = os.path.dirname(index_path)
        # 拡張子を除いたベース名を取得
        base_name = os.path.splitext(os.path.basename(index_path))[0]
        # LangChainのFAISSはディレクトリに保存する
        self.db.save_local(os.path.join(index_dir, base_name))
    
    @classmethod
    def load(cls, index_path: str, embedding_model: Any, metadata_path: str = None):
        """インデックスとメタデータを読み込む
        
        Args:
            index_path: インデックスファイルのパス
            embedding_model: 埋め込みモデル
            metadata_path: 未使用
            
        Returns:
            FAISSIndexer インスタンス
        """
        # ディレクトリ名を取得
        index_dir = os.path.dirname(index_path)
        # 拡張子を除いたベース名を取得
        base_name = os.path.splitext(os.path.basename(index_path))[0]
        
        indexer = cls(embedding_model=embedding_model)
        # LangChainのFAISSを読み込む
        indexer.db = FAISS.load_local(
            os.path.join(index_dir, base_name),
            embedding_model.embeddings,
            allow_dangerous_deserialization=True
        )
        
        return indexer
