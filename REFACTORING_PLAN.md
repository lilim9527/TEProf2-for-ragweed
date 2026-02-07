# TEProf2 Refactoring Plan - Python 3.12+ Migration

## 📁 Proposed Modern Project Structure

```
teprof2/
├── pyproject.toml              # Modern Python packaging (PEP 518)
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml             # Conda environment
│
├── src/
│   └── teprof2/
│       ├── __init__.py
│       ├── __version__.py
│       │
│       ├── core/               # Core data structures and utilities
│       │   ├── __init__.py
│       │   ├── genome_interval.py      # GenomeIntervalHandler class
│       │   ├── config.py               # Configuration management
│       │   ├── exceptions.py           # Custom exceptions
│       │   └── validators.py           # Input validation
│       │
│       ├── io/                 # Input/Output handlers
│       │   ├── __init__.py
│       │   ├── gtf_parser.py           # GTF/GFF3 parsing with pandas
│       │   ├── bam_handler.py          # BAM file handling with pysam
│       │   ├── bed_handler.py          # BED file handling
│       │   └── reference_loader.py     # Reference genome data loading
│       │
│       ├── annotation/         # Annotation logic
│       │   ├── __init__.py
│       │   ├── te_annotator.py         # TE annotation (replaces rmskhg38_annotate*)
│       │   ├── gene_annotator.py       # Gene annotation
│       │   └── intron_annotator.py     # Intron annotation
│       │
│       ├── quantification/     # Expression quantification
│       │   ├── __init__.py
│       │   ├── tpm_calculator.py       # TPM/FPKM calculation
│       │   ├── coverage_analyzer.py    # Coverage analysis
│       │   └── expression_processor.py # Expression processing
│       │
│       ├── analysis/           # Analysis modules
│       │   ├── __init__.py
│       │   ├── promoter_finder.py      # TE promoter analysis
│       │   ├── fusion_detector.py      # TE-gene fusion detection
│       │   └── statistics.py           # Statistical analysis
│       │
│       ├── utils/              # Utility functions
│       │   ├── __init__.py
│       │   ├── interval_ops.py         # Interval operations
│       │   ├── parallel.py             # Parallel processing utilities
│       │   └── logging_config.py       # Logging configuration
│       │
│       └── cli/                # Command-line interface
│           ├── __init__.py
│           ├── main.py                 # Main CLI entry point
│           ├── annotate.py             # Annotation commands
│           └── quantify.py             # Quantification commands
│
├── tests/                      # Unit and integration tests
│   ├── __init__.py
│   ├── test_genome_interval.py
│   ├── test_annotation.py
│   ├── test_quantification.py
│   └── fixtures/               # Test data
│
├── scripts/                    # Standalone scripts
│   ├── convert_legacy_data.py  # Convert old pickle files
│   └── benchmark.py            # Performance benchmarking
│
├── docs/                       # Documentation
│   ├── installation.md
│   ├── usage.md
│   ├── api_reference.md
│   └── migration_guide.md      # Guide for users migrating from v1
│
└── examples/                   # Example workflows
    ├── ambrosia_workflow.py    # Ragweed-specific example
    └── human_workflow.py       # Human genome example
```

## 🎯 Key Architectural Improvements

### 1. **Separation of Concerns**
- **Core**: Data structures independent of I/O
- **IO**: All file operations isolated
- **Annotation/Quantification**: Business logic separated
- **CLI**: User interface layer

### 2. **Dependency Injection**
- Configuration passed explicitly
- No global state
- Testable components

### 3. **Type Safety**
- Full type hints (PEP 484)
- Runtime validation with pydantic
- Static type checking with mypy

### 4. **Error Handling**
- Custom exception hierarchy
- Graceful degradation for missing contigs
- Detailed error messages with context

## 📦 Modern Dependencies

```toml
[project]
name = "teprof2"
version = "2.0.0"
requires-python = ">=3.12"
dependencies = [
    "pysam>=0.22.0",           # BAM/SAM handling
    "pandas>=2.2.0",           # Tabular data
    "polars>=0.20.0",          # Fast dataframe operations
    "pyarrow>=15.0.0",         # Efficient data storage
    "intervaltree>=3.1.0",     # Interval operations
    "pybedtools>=0.9.1",       # BED file operations
    "numpy>=1.26.0",           # Numerical operations
    "scipy>=1.12.0",           # Statistical functions
    "pydantic>=2.6.0",         # Data validation
    "typer>=0.9.0",            # CLI framework
    "rich>=13.7.0",            # Beautiful terminal output
    "loguru>=0.7.2",           # Better logging
    "tqdm>=4.66.0",            # Progress bars
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0.0",
    "pytest-cov>=4.1.0",
    "mypy>=1.8.0",
    "ruff>=0.2.0",             # Fast linter/formatter
    "black>=24.0.0",
]
```

## 🔧 Migration Strategy

### Phase 1: Core Infrastructure (Week 1-2)
1. Set up new project structure
2. Implement GenomeIntervalHandler
3. Create configuration system
4. Set up logging and error handling

### Phase 2: I/O Layer (Week 2-3)
1. GTF/GFF3 parser with pandas
2. BAM handler with pysam
3. BED file handler
4. Reference data loader

### Phase 3: Annotation Logic (Week 3-4)
1. Port TE annotation logic
2. Implement safe dictionary access
3. Dynamic column detection
4. Contig-aware processing

### Phase 4: Quantification (Week 4-5)
1. TPM/FPKM calculator
2. Coverage analyzer
3. Expression processor

### Phase 5: Testing & Optimization (Week 5-6)
1. Unit tests for all modules
2. Integration tests with Ambrosia data
3. Performance optimization
4. Documentation

## 🚀 Performance Optimizations

### 1. **Memory Efficiency**
- Use `__slots__` for frequently instantiated classes
- Stream processing for large files
- Lazy loading of reference data
- Polars for large dataframes (faster than pandas)

### 2. **Parallel Processing**
- `concurrent.futures` for I/O-bound tasks
- `multiprocessing` for CPU-bound tasks
- Batch processing for millions of reads

### 3. **Caching**
- LRU cache for repeated lookups
- Memoization for expensive computations
- Pickle replacement with faster formats (parquet, arrow)

## 🛡️ Robustness Improvements

### 1. **Safe Dictionary Access**
```python
# Old (Python 2.7)
value = dictionary[key]  # KeyError if missing

# New (Python 3.12+)
value = dictionary.get(key, default_value)
# Or with defaultdict
from collections import defaultdict
safe_dict = defaultdict(lambda: default_value)
```

### 2. **Dynamic Column Detection**
```python
# Old: Hardcoded indices
intron_num = line.split("; ")[8].split(" ")[1]  # IndexError prone

# New: Schema-aware parsing
def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """Parse GTF attributes into dictionary."""
    attrs = {}
    for item in attr_string.split(";"):
        if item.strip():
            key, value = item.strip().split(" ", 1)
            attrs[key] = value.strip('"')
    return attrs

# Access safely
attrs = parse_gtf_attributes(line)
intron_num = attrs.get("intron_number", "0")
```

### 3. **Contig-Safe Operations**
```python
# Old: Assumes all contigs have data
te_data = contig_dict[contig_name]  # KeyError for missing contigs

# New: Graceful handling
te_data = contig_dict.get(contig_name, [])
if not te_data:
    logger.warning(f"No TE data for contig {contig_name}, skipping")
    continue
```

## 📊 Ambrosia-Specific Optimizations

### 1. **Fragmented Genome Handling**
- Efficient storage for thousands of small contigs
- Batch processing to reduce overhead
- Memory-mapped files for large datasets

### 2. **TE Promoter Analysis**
- Configurable TSS upstream window (default: 2kb)
- Efficient interval tree for overlap detection
- Parallel processing per contig

### 3. **Non-Standard GTF Support**
- Flexible attribute parsing
- Column mapping configuration
- Validation with detailed error messages

## 🔍 Key Problem Solutions

### Problem 1: KeyError for Missing Contigs
**Root Cause**: Direct dictionary access without checking existence
**Solution**: Use `.get()` with defaults, defaultdict, or explicit checks

### Problem 2: IndexError from Hardcoded Indices
**Root Cause**: Assuming fixed column structure
**Solution**: Parse attributes into dictionaries, use named access

### Problem 3: Memory Issues with Large Genomes
**Root Cause**: Loading entire datasets into memory
**Solution**: Stream processing, chunking, memory-mapped files

### Problem 4: Slow Processing
**Root Cause**: Sequential processing, inefficient data structures
**Solution**: Parallel processing, intervaltree, optimized algorithms
