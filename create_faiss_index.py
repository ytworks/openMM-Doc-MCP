#!/usr/bin/env python
"""
マークダウンドキュメントからFAISSベクターインデックスを作成するスクリプト
"""
import os
import argparse
import glob
import traceback
from typing import List
import logging

from src.vector_db.indexer import MarkdownChunker, Embedder, FAISSIndexer

# ロギングの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_markdown_files(docs_path: str) -> List[str]:
    """
    指定されたディレクトリ内のマークダウンファイルを取得する
    
    Args:
        docs_path: マークダウンファイルが格納されているディレクトリパス
        
    Returns:
        マークダウンファイルパスのリスト
    """
    md_files = glob.glob(os.path.join(docs_path, "**/*.md"), recursive=True)
    logger.info(f"{len(md_files)}件のマークダウンファイルを検出しました")
    return md_files


def main():
    """
    メイン処理
    """
    # コマンドライン引数の解析
    parser = argparse.ArgumentParser(description="マークダウンドキュメントからFAISSベクターDBを作成")
    parser.add_argument("--docs-path", type=str, default="specs", 
                        help="マークダウンファイルが格納されているディレクトリパス")
    parser.add_argument("--output-dir", type=str, default="data/indices",
                        help="出力ディレクトリ")
    parser.add_argument("--index-name", type=str, default="docs",
                        help="インデックスの名前")
    parser.add_argument("--chunk-size", type=int, default=500,
                        help="チャンクサイズ（文字数）")
    parser.add_argument("--chunk-overlap", type=int, default=50,
                        help="チャンクオーバーラップ（文字数）")
    parser.add_argument("--embedding-model", type=str, 
                        default="intfloat/multilingual-e5-large",
                        help="埋め込みモデル")
    args = parser.parse_args()
    
    # 出力ディレクトリの作成
    os.makedirs(args.output_dir, exist_ok=True)
    
    # インデックスファイルのパス
    index_path = os.path.join(args.output_dir, f"{args.index_name}.index")
    
    # マークダウンファイルの取得
    md_files = get_markdown_files(args.docs_path)
    if not md_files:
        logger.error(f"マークダウンファイルが見つかりませんでした: {args.docs_path}")
        return
    
    # チャンキング
    logger.info("マークダウンのチャンキングを開始します...")
    chunker = MarkdownChunker(chunk_size=args.chunk_size, chunk_overlap=args.chunk_overlap)
    all_chunks = []
    
    for md_file in md_files:
        logger.info(f"処理中: {md_file}")
        chunks = chunker.process_markdown(md_file)
        logger.info(f"{len(chunks)}個のチャンクを抽出しました")
        all_chunks.extend(chunks)
    
    logger.info(f"合計{len(all_chunks)}個のチャンクを抽出しました")
    
    # 埋め込み
    logger.info(f"埋め込みモデルを初期化: {args.embedding_model}")
    embedder = Embedder(model_name=args.embedding_model)
    
    logger.info("テキストの埋め込みを実行中...")
    embeddings, chunks = embedder.embed_chunks(all_chunks)
    
    # FAISSインデックスの作成
    logger.info("FAISSインデックスの作成")
    indexer = FAISSIndexer(embedding_model=embedder)
    
    # 埋め込みとメタデータを追加
    indexer.add_embeddings(embeddings, chunks)
    
    # インデックスを保存
    logger.info(f"インデックスを保存: {index_path}")
    indexer.save(index_path)
    
    logger.info("処理が完了しました")
    logger.info(f"インデックスディレクトリ: {os.path.join(args.output_dir, args.index_name)}")
    logger.info(f"総チャンク数: {len(chunks)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"エラーが発生しました: {e}")
        traceback.print_exc()
