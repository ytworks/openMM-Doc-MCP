# OpenMM Documentation MCP Server

A search server for OpenMM molecular dynamics simulation documentation. This server complies with the Model Context Protocol (MCP) and is optimized for integration with LLMs.

## Overview

This server vectorizes OpenMM documentation and makes it searchable. Users can send queries to find relevant document sections. It uses FAISS vector database to achieve fast search capabilities.

OpenMM is a high-performance simulation library for biological molecular systems, and this server makes it easier to search through its documentation. It's particularly useful for:

- Searching for documentation on specific molecular dynamics methods
- Finding usage instructions for OpenMM functions and classes
- Retrieving information about simulation parameter settings

## Features

- Document Search: Find relevant documents based on natural language queries
- Index Information: Retrieve information about the vector database
- MCP Compliant: Optimized for AI/LLM integration
- Fast Search: Efficient vector search using FAISS
- Flexible Queries: Support for both Japanese and English queries

## Technology Stack

This project uses the following key technologies:

- **FastAPI**: REST API implementation
- **FAISS**: Efficient vector search indexing
- **LangChain**: Document processing and LLM integration
- **HuggingFace Models**: Document embedding generation
- **Pytest**: Test automation

## Setup

### Prerequisites

- Python 3.9 or higher
- `uv` package manager (recommended) or `pip`
- Minimum 8GB RAM (16GB+ recommended for index creation)

### Required Packages

Main dependencies:

```
fastapi==0.110.0
uvicorn==0.28.0
faiss-cpu==1.8.0
langchain==0.1.12
langchain-community==0.0.29
langchain-text-splitters==0.0.1
sentence-transformers==2.5.1
pydantic==2.6.4
pytest==8.3.5
```

All dependencies are listed in the `requirements.txt` file.

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/openMM-Doc-MCP.git
cd openMM-Doc-MCP

# Create and activate virtual environment (using uv)
uv venv

# Install packages with uv (recommended)
uv pip install -r requirements.txt

# Or, if using pip
# python -m venv .venv
# source .venv/bin/activate  # For Linux
# .venv\Scripts\activate     # For Windows
# pip install -r requirements.txt
```

### Environment Configuration

You can set environment variables as needed:

```bash
# Set server port (default is 8080)
export MCP_SERVER_PORT=8888

# Set index directory (optional)
export MCP_INDEX_DIR="/path/to/custom/index"
```

For Windows environment variables:

```
set MCP_SERVER_PORT=8888
set MCP_INDEX_DIR=C:\path\to\custom\index
```

### Creating the Index

To create an index of the documents, run the following command:

```bash
uv run python create_faiss_index.py
```

Optional parameters:

```bash
# Specify a custom document directory
uv run python create_faiss_index.py --docs_dir "/path/to/docs"

# Specify a custom output directory
uv run python create_faiss_index.py --output_dir "/path/to/output"

# Specify a custom embedding model
uv run python create_faiss_index.py --embedding_model "intfloat/multilingual-e5-large"
```

Index creation may take time. The process uses significant system resources (especially memory). By default, index files are created in the `data/indices/docs/` directory.

## Usage

### Starting the Server

```bash
uv run python server.py
```

By default, the server listens at http://localhost:8080.

### Executing Queries

To search directly from the command line:

```bash
uv run python search_molecular_simulation.py "basic principles of molecular dynamics"
```

Or you can send an HTTP request:

```
POST http://localhost:8080/query
Content-Type: application/json

{
  "query": "basic principles of molecular dynamics",
  "top_k": 5
}
```

## Testing

This project includes tests for the vector database and MCP server.

### Running Tests

To run all tests:

```bash
uv run -m pytest
```

To run specific tests:

```bash
# Run only server tests
uv run -m pytest tests/test_server.py

# Run vector database tests
uv run -m pytest src/vector_db/tests/
```

Add the `-v` option for detailed output:

```bash
uv run -m pytest tests/test_server.py -v
```

### Test Structure

- `tests/test_server.py`: MCP server functionality tests
- `src/vector_db/tests/`: Vector database related tests
  - `test_indexer.py`: Indexer tests
  - `test_retriever.py`: Retriever tests

## API Specification

### Endpoints

This server provides the following API endpoints:

1. **Search Endpoint**
   - URL: `/query`
   - Method: POST
   - Request:
     ```json
     {
       "query": "text to search for",
       "top_k": 5,  // optional, default is 5
       "index_path": "/path/to/index"  // optional
     }
     ```
   - Response:
     ```json
     {
       "query": "text to search for",
       "results_count": 5,
       "results": [
         {
           "content": "Document content...",
           "metadata": {
             "source": "filename"
           },
           "score": 0.95
         },
         // ...other results
       ],
       "metadata": {
         "top_k": 5,
         "index_path": "/path/to/index"
       }
     }
     ```

2. **Index Information Endpoint**
   - URL: `/index-info`
   - Method: GET
   - Query Parameters:
     - `index_path`: Path to the index (optional)
   - Response:
     ```json
     {
       "exists": true,
       "path": "/path/to/index",
       "index_file": "/path/to/index/index.faiss",
       "index_pkl": "/path/to/index/index.pkl",
       "index_file_size": 1024,
       "index_pkl_size": 2048,
       "document_count": 100,
       "embedding_model": "model-name"
     }
     ```

3. **Health Check Endpoint**
   - URL: `/health`
   - Method: GET
   - Response:
     ```json
     {
       "status": "ok",
       "version": "1.0.0"
     }
     ```

Detailed API specifications are available in the following files:
- English: [specs/apispec_en.md](specs/apispec_en.md)
- Japanese: [specs/apispec_ja.md](specs/apispec_ja.md)

## Troubleshooting

### Common Issues

1. **Server Won't Start**
   - Check that dependencies are correctly installed
   - Verify that the port is not already in use
   - Check logs for specific error messages

2. **No Search Results Returned**
   - Verify that index files have been created correctly
   - Ensure that the query is not empty
   - Check that the index path is correct

3. **Memory Errors**
   - Index creation requires sufficient memory (16GB+ recommended)
   - If using large models, try switching to smaller ones

### Debugging

To enable detailed logging, set the environment variable:

```bash
export DEBUG=true
uv run python server.py
```

### Test Execution Issues

1. **NumPy Errors**
   - If you encounter errors with `libgfortran.5.dylib`, running tests with `uv run` can resolve this
   - In some environments, you may need to install the library with `conda install -c conda-forge libgfortran`

2. **Tests Timeout**
   - Increase test timeout: `uv run -m pytest --timeout=30`

## Performance Tuning

### Latency Optimization

- Use smaller indices
- Enable GPU acceleration with `faiss-gpu` (if you have compatible GPUs)
- Adjust indexing parameters

### Memory Usage Optimization

For processing large document collections, you can optimize memory usage with these settings:

```python
# Adjust chunk size in create_faiss_index.py
--chunk_size 256  # default is 512
--chunk_overlap 20  # default is 50
```

## Developer Information

### Coding Conventions

Please follow these coding conventions for this project:

- PEP 8 style guide
- Provide appropriate docstrings for functions and classes
- Maintain test coverage

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Create a Pull Request

### Extension Guide

#### Adding New Vector Stores

Currently using FAISS, but you can add other vector stores:

1. Add a new indexer class in `src/vector_db/indexer.py`
2. Add a corresponding retriever class in `src/vector_db/retriever.py`
3. Create tests for the new indexer/retriever
4. Make it selectable via environment variables or command line arguments

## Directory Structure

```
openMM-Doc-MCP/
├── create_faiss_index.py   # Index creation script
├── search_molecular_simulation.py # CLI search
├── server.py               # MCP server implementation
├── data/
│   └── indices/
│       └── docs/
│           ├── index.faiss # FAISS index file
│           └── index.pkl   # Pickle index file
├── specs/
│   ├── apispec_en.md       # API specification (English)
│   └── apispec_ja.md       # API specification (Japanese)
├── src/
│   └── vector_db/          # Vector database related modules
│       ├── indexer.py      # Indexer implementation
│       ├── retriever.py    # Retriever implementation
│       └── tests/          # Vector DB tests
│           ├── conftest.py
│           ├── test_indexer.py
│           └── test_retriever.py
└── tests/
    └── test_server.py      # Server tests
```

## License

This project is provided under the [LICENSE](LICENSE).