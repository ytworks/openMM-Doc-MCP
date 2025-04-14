#!/usr/bin/env python
"""
分子シミュレーションに関連する情報をベクターDBから検索するサンプルスクリプト
"""
import os
import argparse
import time
from src.vector_db.retriever import FAISSRetriever


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(
        description="分子シミュレーションに関連する情報をベクターDBから検索します"
    )
    parser.add_argument(
        "--query",
        type=str,
        default="分子シミュレーション",
        help="検索クエリ（デフォルト: 分子シミュレーション）"
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="表示する結果の数（デフォルト: 5）"
    )
    parser.add_argument(
        "--index-dir",
        type=str,
        default="data/indices/docs",
        help="インデックスディレクトリのパス（デフォルト: data/indices/docs）"
    )
    args = parser.parse_args()

    # インデックスパスの確認
    index_path = os.path.abspath(args.index_dir)
    if not os.path.exists(index_path):
        print(f"エラー: インデックスディレクトリが存在しません: {index_path}")
        return

    print(f"検索クエリ: '{args.query}'")
    print(f"検索結果数: {args.top_k}")
    print(f"インデックスパス: {index_path}")
    print("検索中...")

    # リトリーバーの初期化
    start_time = time.time()
    print("ベクターDBを読み込み中...")
    retriever = FAISSRetriever(index_path)
    print(f"ベクターDB読み込み完了 ({time.time() - start_time:.2f}秒)")
    
    # 検索実行
    print(f"'{args.query}'で検索を実行中...")
    search_start = time.time()
    results = retriever.search_by_text(args.query, top_k=args.top_k)
    print(f"検索完了 ({time.time() - search_start:.2f}秒)")
    
    # ファイルに結果を保存
    result_file = "search_results.txt"
    with open(result_file, "w", encoding="utf-8") as f:
        f.write(f"検索クエリ: '{args.query}'\n")
        f.write(f"検索結果 ({len(results)}件):\n\n")
        
        for i, result in enumerate(results, 1):
            content = result["content"]
            source = result["metadata"].get("source", "不明")
            
            f.write(f"【結果 {i}】\n")
            f.write(f"ソース: {source}\n")
            f.write(f"内容: {content}\n")
            f.write("-" * 80 + "\n")
    
    # 結果の表示
    print(f"検索結果 ({len(results)}件):\n")
    
    for i, result in enumerate(results, 1):
        content = result["content"]
        source = result["metadata"].get("source", "不明")
        
        print(f"【結果 {i}】")
        print(f"ソース: {source}")
        print(f"内容: {content}")
        print("-" * 80)
    
    print(f"\n結果はファイル '{result_file}' にも保存されました")


if __name__ == "__main__":
    main()