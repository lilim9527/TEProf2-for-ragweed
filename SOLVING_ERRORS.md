# Solving KeyError and IndexError Problems

## Problem Analysis

The original TEProf2 code frequently crashes with:
1. **KeyError**: `KeyError: 'contigX'` - when accessing dictionaries for contigs without TE data
2. **IndexError**: `IndexError: list index out of range` - when using hardcoded indices on variable-length data

## Root Causes

### 1. KeyError: Direct Dictionary Access

```python
# Original code (Python 2.7) - CRASHES on missing keys
te_data = contig_dict[contig_name]  # KeyError if contig_name not in dict
```

**Why it fails**: Fragmented genomes like Ambrosia have thousands of contigs. Not all contigs have TE annotations, but the code assumes they do.

### 2. IndexError: Hardcoded Indices

```python
# Original code - CRASHES on non-standard GTF
intron_num = line.split("; ")[8].split(" ")[1]  # IndexError if < 9 elements
```

**Why it fails**: Non-model organisms often have non-standard GTF files with different numbers of attributes.

## Solutions

### Solution 1: Safe Dictionary Access

#### Method A: `.get()` with Default Value

```python
# ✅ SAFE: Returns default if key missing
te_data = contig_dict.get(contig_name, [])

if not te_data:
    logger.warning(f"No TE data for contig {contig_name}, skipping")
    continue

# Process te_data safely
```

#### Method B: `defaultdict`

```python
from collections import defaultdict

# ✅ SAFE: Automatically creates empty list for missing keys
contig_dict = defaultdict(list)

# No KeyError possible!
te_data = contig_dict[contig_name]  # Returns [] if missing
```

#### Method C: Explicit Check

```python
# ✅ SAFE: Check before access
if contig_name in contig_dict:
    te_data = contig_dict[contig_name]
    # Process te_data
else:
    logger.warning(f"Contig {contig_name} not found in TE database")
    continue
```

### Solution 2: Dynamic Attribute Parsing

#### Problem: Hardcoded Indices

```python
# ❌ FRAGILE: Assumes fixed structure
def parse_gtf_old(line):
    fields = line.split('\t')
    attributes = fields[8].split('; ')
    intron_num = attributes[8].split(' ')[1]  # CRASHES if < 9 attributes
    return intron_num
```

#### Solution: Parse into Dictionary

```python
# ✅ ROBUST: Handles variable attributes
def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attributes into dictionary.

    Handles:
    - Variable number of attributes
    - Different attribute orders
    - Missing attributes
    - Various GTF formats
    """
    attributes = {}

    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue

        # Split on first whitespace
        parts = item.split(None, 1)
        if len(parts) == 2:
            key, value = parts
            # Remove quotes
            value = value.strip('"').strip("'")
            attributes[key] = value

    return attributes

# Usage
attrs = parse_gtf_attributes(gtf_line[8])

# Safe access with default
intron_num = attrs.get("intron_number", "0")
gene_name = attrs.get("gene_name", "unknown")
```

### Solution 3: Schema Validation

```python
from typing import Optional
from pydantic import BaseModel, validator

class GTFRecord(BaseModel):
    """Validated GTF record."""

    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attributes: dict[str, str]

    @validator('start', 'end')
    def validate_coordinates(cls, v):
        if v < 0:
            raise ValueError(f"Coordinate cannot be negative: {v}")
        return v

    @validator('strand')
    def validate_strand(cls, v):
        if v not in ('+', '-', '.'):
            raise ValueError(f"Invalid strand: {v}")
        return v

    def get_attribute(self, key: str, default: str = "None") -> str:
        """Safely get attribute with default."""
        return self.attributes.get(key, default)

# Usage
try:
    record = GTFRecord(
        seqname=fields[0],
        source=fields[1],
        feature=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=fields[5],
        strand=fields[6],
        frame=fields[7],
        attributes=parse_gtf_attributes(fields[8])
    )

    # Safe access
    intron_num = record.get_attribute("intron_number", "0")

except ValueError as e:
    logger.error(f"Invalid GTF record: {e}")
    continue
```

## Complete Example: Refactored Annotation Function

### Before (Python 2.7 - Fragile)

```python
# ❌ CRASHES on fragmented genomes
def annotate_transcript_old(transcript_id, contig_dict, gencode_dict):
    # KeyError if contig not in dict
    te_list = contig_dict[transcript_id]

    # IndexError if wrong format
    gene_name = te_list[0].split("; ")[8].split(" ")[1]

    # KeyError if gene not in gencode
    gene_info = gencode_dict[gene_name]

    return gene_info
```

### After (Python 3.12+ - Robust)

```python
# ✅ ROBUST: Handles all edge cases
from typing import Optional
import logging

logger = logging.getLogger(__name__)

def annotate_transcript_new(
    transcript_id: str,
    contig_dict: dict[str, list],
    gencode_dict: dict[str, dict]
) -> Optional[dict]:
    """
    Annotate transcript with TE information.

    Handles:
    - Missing contigs (no KeyError)
    - Variable GTF formats (no IndexError)
    - Missing genes (no KeyError)

    Args:
        transcript_id: Transcript identifier
        contig_dict: Dictionary of contig -> TE list
        gencode_dict: Dictionary of gene -> info

    Returns:
        Annotation dictionary or None if annotation fails
    """
    # SAFE: Use .get() with default
    te_list = contig_dict.get(transcript_id, [])

    if not te_list:
        logger.debug(f"No TE data for transcript {transcript_id}")
        return None

    # SAFE: Parse attributes into dictionary
    try:
        attributes = parse_gtf_attributes(te_list[0])
        gene_name = attributes.get("gene_name", "unknown")
    except (IndexError, ValueError) as e:
        logger.warning(f"Failed to parse TE attributes for {transcript_id}: {e}")
        return None

    # SAFE: Check before access
    if gene_name not in gencode_dict:
        logger.debug(f"Gene {gene_name} not found in Gencode")
        return None

    gene_info = gencode_dict[gene_name]

    return {
        "transcript_id": transcript_id,
        "gene_name": gene_name,
        "te_count": len(te_list),
        "gene_info": gene_info
    }
```

## Real-World Example: Ambrosia Genome

### Problem Scenario

```python
# Ambrosia genome: 50,000 contigs
# Only 15,000 contigs have TE annotations
# GTF has variable attribute counts (6-12 attributes per line)

contigs = ["contig_1", "contig_2", ..., "contig_50000"]
te_dict = {
    "contig_1": [...],  # Has TEs
    "contig_5": [...],  # Has TEs
    # contig_2, contig_3, contig_4 have NO TEs!
}
```

### Old Code (Crashes)

```python
# ❌ CRASHES on contig_2, contig_3, contig_4
for contig in contigs:
    te_data = te_dict[contig]  # KeyError!
    process_te_data(te_data)
```

### New Code (Robust)

```python
# ✅ WORKS: Handles missing contigs gracefully
from collections import defaultdict

# Option 1: Use defaultdict
te_dict = defaultdict(list)
# ... populate te_dict ...

for contig in contigs:
    te_data = te_dict[contig]  # Returns [] if missing
    if te_data:
        process_te_data(te_data)
    else:
        logger.debug(f"No TE data for {contig}, skipping")

# Option 2: Use .get()
for contig in contigs:
    te_data = te_dict.get(contig, [])
    if te_data:
        process_te_data(te_data)

# Option 3: Check membership
for contig in contigs:
    if contig in te_dict:
        te_data = te_dict[contig]
        process_te_data(te_data)
```

## Testing Strategy

### Unit Tests for Edge Cases

```python
import pytest

def test_missing_contig():
    """Test handling of missing contig."""
    contig_dict = {"contig_1": ["data"]}

    # Should not raise KeyError
    result = annotate_transcript_new("contig_999", contig_dict, {})
    assert result is None

def test_malformed_gtf():
    """Test handling of malformed GTF attributes."""
    # GTF with only 3 attributes (not standard 9+)
    malformed_attr = "gene_id ENSG001; gene_name TP53"

    attrs = parse_gtf_attributes(malformed_attr)
    assert attrs.get("gene_id") == "ENSG001"
    assert attrs.get("gene_name") == "TP53"
    assert attrs.get("missing_attr", "default") == "default"

def test_empty_te_list():
    """Test handling of empty TE list."""
    contig_dict = {"contig_1": []}

    result = annotate_transcript_new("contig_1", contig_dict, {})
    assert result is None
```

## Best Practices Summary

### ✅ DO

1. **Use `.get()` for dictionary access**
   ```python
   value = dict.get(key, default)
   ```

2. **Use `defaultdict` for accumulation**
   ```python
   from collections import defaultdict
   data = defaultdict(list)
   ```

3. **Parse structured data into dictionaries**
   ```python
   attrs = parse_attributes(string)
   value = attrs.get("key", "default")
   ```

4. **Validate inputs**
   ```python
   if not data:
       logger.warning("No data found")
       return None
   ```

5. **Use type hints**
   ```python
   def func(data: dict[str, list]) -> Optional[dict]:
       ...
   ```

### ❌ DON'T

1. **Don't use direct dictionary access without checking**
   ```python
   value = dict[key]  # Can raise KeyError
   ```

2. **Don't use hardcoded indices**
   ```python
   value = list[8]  # Can raise IndexError
   ```

3. **Don't assume data structure**
   ```python
   # Assumes exactly 9 elements
   parts = line.split("; ")
   value = parts[8]
   ```

4. **Don't ignore errors silently**
   ```python
   try:
       value = dict[key]
   except KeyError:
       pass  # Silent failure!
   ```

## Conclusion

By following these patterns, TEProf2 v2.0 can handle:
- ✅ Fragmented genomes with thousands of contigs
- ✅ Non-standard GTF files
- ✅ Missing annotations
- ✅ Variable data formats

The code is now **robust**, **maintainable**, and **production-ready** for non-model organism analysis.