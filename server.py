#!/usr/bin/env python
"""
MCP server for document similarity search using a vector database
"""
import json
import logging
import sys
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

# 必要なライブラリのインポート
from mcp.server.fastmcp import FastMCP
from src.vector_db.retriever import FAISSRetriever

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Initialize MCP server
mcp = FastMCP("OpenMM Document Search Service")

# プロジェクトのルートディレクトリを取得
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Define the path to the FAISS index with absolute path
DEFAULT_INDEX_PATH = os.path.join(ROOT_DIR, "data/indices/docs")
logger.info(f"Setting index path to: {DEFAULT_INDEX_PATH}")

# Initialize the retriever
retriever = None

def get_retriever():
    """Lazy load the retriever to avoid loading the model on import"""
    global retriever
    if retriever is None:
        try:
            retriever = FAISSRetriever(str(DEFAULT_INDEX_PATH))
            logger.info(f"Initialized FAISSRetriever with index from {DEFAULT_INDEX_PATH}")
        except Exception as e:
            logger.error(f"Failed to initialize retriever: {e}")
            raise
    return retriever

@mcp.tool()
def search_documents(
    query: str,
    top_k: int = 5,
    index_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Search for similar documents based on a query string
    
    Args:
        query: The search query text
        top_k: Number of results to return (default 5)
        index_path: Optional custom path to FAISS index
        
    Returns:
        Dict: Dictionary containing search results
    """
    try:
        if not query:
            return {"error": "No query provided"}
        
        # Use custom index path if provided
        if index_path:
            try:
                custom_retriever = FAISSRetriever(index_path)
                results = custom_retriever.search_by_text(query, top_k=top_k)
            except Exception as e:
                return {"error": f"Failed to load custom index from {index_path}: {str(e)}"}
        else:
            # Use the default retriever
            try:
                results = get_retriever().search_by_text(query, top_k=top_k)
            except Exception as e:
                return {"error": f"Search failed: {str(e)}"}
        
        # Format response
        return {
            "query": query,
            "results_count": len(results),
            "results": results,
            "metadata": {
                "top_k": top_k,
                "index_path": index_path or str(DEFAULT_INDEX_PATH)
            }
        }
            
    except Exception as e:
        logger.exception("Error occurred during document search")
        return {"error": f"An error occurred: {str(e)}"}

@mcp.tool()
def get_index_info(index_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Get information about the currently loaded vector database index
    
    Args:
        index_path: Optional custom path to FAISS index
        
    Returns:
        Dict: Dictionary containing index information
    """
    try:
        path_to_use = index_path or str(DEFAULT_INDEX_PATH)
        
        # Check if the index exists
        index_file = Path(path_to_use) / "index.faiss"
        index_pkl = Path(path_to_use) / "index.pkl"
        
        if not index_file.exists() or not index_pkl.exists():
            return {
                "error": f"Index files not found at {path_to_use}",
                "exists": False
            }
        
        # Basic index info
        index_info = {
            "exists": True,
            "path": path_to_use,
            "index_file": str(index_file),
            "index_pkl": str(index_pkl),
            "index_file_size": index_file.stat().st_size,
            "index_pkl_size": index_pkl.stat().st_size
        }
        
        # Try to load the index to get more information
        try:
            if index_path:
                retriever = FAISSRetriever(index_path)
            else:
                retriever = get_retriever()
                
            # Add more details if available
            db = retriever.db
            if hasattr(db, "index"):
                if hasattr(db.index, "ntotal"):
                    index_info["document_count"] = db.index.ntotal
                
            index_info["embedding_model"] = retriever.embeddings.model_name
            
        except Exception as e:
            index_info["load_error"] = str(e)
        
        return index_info
        
    except Exception as e:
        logger.exception("Error occurred while retrieving index info")
        return {"error": f"An error occurred: {str(e)}"}

if __name__ == "__main__":
    try:
        # Initialize the retriever to validate it works
        retriever = get_retriever()
        print(f"Successfully initialized FAISSRetriever with index from {DEFAULT_INDEX_PATH}", file=sys.stderr)
        
        print("Starting MCP server for document similarity search...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"Server startup error: {str(e)}", file=sys.stderr)
        sys.exit(1)