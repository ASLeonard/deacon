---
description: Plan for adding paired-end index building, frequency filtering, index subtraction, and dual-index competitive binning for long reads with ratio-based haplotype assignment
---

# Plan: Multi-Index Long-Read Binning with Frequency Filtering

This plan adds four major capabilities to deacon: (1) paired-end FASTQ index building, (2) frequency-based minimizer filtering during index construction, (3) index subtraction to retain exclusive minimizers, and (4) dual-index competitive binning for long reads with ratio-based haplotype assignment.

**Key Design Decisions:**
- Filtering now REQUIRES two indexes for haplotype binning - no backward compatibility with single-index mode
- **Long reads only** - no support for paired-end reads in filtering (only single-end)
- Reads assigned to haplotype 0 (ambiguous), 1, or 2 based on hit counts and ratio thresholds
- Output flags: `--hap0`, `--hap1`, `--hap2` for three haplotype bins
- Ratio calculation: `hits1 / hits2` (if hits2=0 and hits1>0, assign to hap1)
- Original threshold flags (`--abs-threshold`, `--rel-threshold`, `--deplete`) are removed

---

## Feature 1: Build Indexes Using Paired-End FASTQ

### Goal
Allow index building from paired-end FASTQ files where minimizers from both mates are merged into a single index.

### Files to Modify

#### `src/index.rs`
- Add `input2_path: Option<PathBuf>` field to `IndexConfig`
- Modify `build()` function to detect paired-end mode when `input2_path.is_some()`
- When paired: use `PairedParallelProcessor` trait instead of `ParallelProcessor`
- Implement `BuildIndexPairedProcessor` struct that:
  - Maintains per-thread `local_minimizers: RapidHashSet<u128>`
  - In `process_pair()`: extracts minimizers from both mates, adds all to local set
  - In `on_batch_complete()`: merges local set into global `MinimizerSet`
- Merge logic: collect minimizers from mate1 and mate2 into same set (union operation)

#### `src/main.rs`
- Add `--input2 <PATH>` argument to `IndexCommands::Build`
- Pass `input2_path` to `IndexConfig`
- Update help text to document paired-end support

### Implementation Details
- Use existing `paraseq` paired-end reader capabilities
- Single output index file containing minimizers from both mates
- No change to serialization format (still stores u128 minimizers)

---

## Feature 2: Count Minimizer Frequency and Filter by Threshold

### Goal
During index building, track how many times each minimizer occurs and only retain minimizers that appear at least N times in the input dataset.

### Files to Modify

#### `src/lib.rs`
- Add type alias: `pub type RapidHashMap128 = HashMap<u128, u32, FixedRapidHasher>`
- Add `pub struct MinimizerFrequencies { counts: RapidHashMap128 }` with methods:
  - `new() -> Self`
  - `increment(&mut self, minimizer: u128)`
  - `filter_by_threshold(&self, min_freq: u32) -> MinimizerSet` - returns set of minimizers with count ≥ threshold
  - `len(&self) -> usize`
  - `histogram(&self) -> Vec<(u32, usize)>` - returns frequency bins for logging

#### `src/index.rs`
- Add `min_frequency: Option<u32>` field to `IndexConfig` (default None = no filtering)
- Modify `BuildIndexProcessor` to track counts:
  - When `min_frequency.is_some()`: use `MinimizerFrequencies` instead of `MinimizerSet`
  - Each minimizer occurrence increments counter
  - After processing all sequences, call `filter_by_threshold()` to get final set
- Add logging:
  - Print frequency histogram before filtering (e.g., "1-10: 50000, 11-100: 20000, 101+: 5000")
  - Print "Retained X / Y minimizers (Z%) after frequency filtering"

#### `src/main.rs`
- Add `--min-frequency <N>` argument to `IndexCommands::Build`
- Pass `min_frequency` to `IndexConfig`
- Help text: "Only retain minimizers occurring at least N times (default: 1, no filtering)"

### Implementation Details
- Count increments are per occurrence (not per sequence)
- Paired-end mode: each minimizer occurrence in either mate counts separately
- u32 counter supports up to 4 billion occurrences per minimizer
- Histogram output: print all frequency values 1-1023 individually, combine all counts ≥1024 into final "1024+" bin

---

## Feature 3: Index Subtraction for Exclusive Minimizers

### Goal
Create a new index containing only minimizers that appear in index A but not in index B (set difference A - B).

### Files to Modify

#### `src/index.rs`
- Add new function: `pub fn subtract(index_a_path: &Path, index_b_path: &Path, output_path: &Path) -> Result<()>`
  - Load both indexes with `load_header_and_count()` to validate k/w match
  - Load both as `MinimizerSet` via `load_minimizers()`
  - Compute set difference: iterate over A, skip if in B, collect remainder
  - Dump result via `dump_minimizers()`
- Add validation:
  - Error if k or w mismatch between indexes
  - Log: "Index A: X minimizers, Index B: Y minimizers, A-B: Z minimizers (W% retained)"

#### `src/main.rs`
- Add new subcommand: `IndexCommands::Subtract { index_a, index_b, output }`
  - Arguments:
    - `index_a`: First index (PathBuf)
    - `index_b`: Index to subtract (PathBuf)
    - `-o/--output`: Output index (PathBuf)
- Wire to `IndexCommands::execute()` to call `index::subtract()`

### Implementation Details
- Works with existing u128 minimizer sets
- Does not modify original indexes (creates new output)
- Can chain operations: `(A - B) - C` by running subtract twice
- Enables creating species-specific indexes: e.g., `human_exclusive = human - mouse`

---

## Feature 4: Dual-Index Competitive Binning with Ratio Thresholds

### Goal
**BREAKING CHANGE**: Filter mode now requires two indexes and assigns long reads to haplotype 1 (hap1), haplotype 2 (hap2), or ambiguous (hap0) based on hit counts and configurable ratio thresholds.

**Long reads only** - no paired-end support in filtering.

### Assignment Logic

For each read:
1. Extract distinct minimizers from sequence
2. Count hits: `hits1 = count(minimizers in index1)`, `hits2 = count(minimizers in index2)`
3. Compute ratio: `ratio = hits1 / hits2`
4. Assign haplotype:
   - **Haplotype 0 (Ambiguous)** if:
     - Both `hits1 == 0` AND `hits2 == 0` (no hits in either index)
     - Both `hits1 > 0` AND `hits2 > 0` AND `ratio < 1/min_ratio_threshold` AND `ratio > min_ratio_threshold` (hits too similar)
   - **Haplotype 1** if:
     - `hits1 > 0` AND `hits2 == 0` (only hits in index1)
     - `ratio >= min_ratio_threshold` (index1 has sufficiently more hits)
   - **Haplotype 2** if:
     - `hits1 == 0` AND `hits2 > 0` (only hits in index2)
     - `ratio <= 1/min_ratio_threshold` (index2 has sufficiently more hits)

### Files to Modify

#### `src/lib.rs`
- Add enum:
```rust
pub enum ReadAssignment {
    Ambiguous,  // status 0
    Index1,     // status 1
    Index2,     // status 2
}
```

#### `src/filter.rs`
**BREAKING CHANGES**: Simplify and redesign for dual-index mode only

- **Remove** from `FilterConfig`:
  - `abs_threshold`, `rel_threshold` (replaced by ratio threshold)
  - `deplete` flag (not applicable to binning)
  - `output_path`, `output2_path` (replaced by status-based outputs)
  
- **Replace** `FilterConfig` with:
```rust
pub struct FilterConfig {
    // Required: two indexes
    pub index1_path: PathBuf,
    pub index2_path: PathBuf,
    
    // Input (long reads only)
    pub input_path: PathBuf,
    
    // Output (haplotype-based)
    pub hap0_path: Option<PathBuf>,  // ambiguous reads (optional)
    pub hap1_path: PathBuf,          // haplotype 1 reads
    pub hap2_path: PathBuf,          // haplotype 2 reads
    
    // Assignment thresholds
    pub min_ratio_threshold: f64,  // default: 2.0 (ratio = hits1/hits2)
    
    // Processing options
    pub prefix_length: usize,
    pub threads: usize,
    pub compression_level: u32,
    pub rename: bool,
    pub quiet: bool,
    pub debug: bool,
}
```

- **Redesign** `FilterProcessor`:
```rust
struct FilterProcessor {
    config: Arc<FilterConfig>,
    minimizers_a: Arc<MinimizerSet>,  // index1
    minimizers_b: Arc<MinimizerSet>,  // index2
    
    // Three output writers
    writer_hap0: Option<Arc<Mutex<BufWriter>>>,  // ambiguous
    writer_hap1: Arc<Mutex<BufWriter>>,          // haplotype 1
    writer_hap2: Arc<Mutex<BufWriter>>,          // haplotype 2
    
    // Local buffers
    buffers: Buffers,
    local_buffer_hap0: Vec<u8>,
    local_buffer_hap1: Vec<u8>,
    local_buffer_hap2: Vec<u8>,
}
```

- **Add method**: `fn assign_read(&self, minimizers: &[u128]) -> ReadAssignment`
```rust
fn assign_read(&self, minimizers: &[u128]) -> ReadAssignment {
    let hits1 = minimizers.iter().filter(|m| self.minimizers_a.contains(m)).count();
    let hits2 = minimizers.iter().filter(|m| self.minimizers_b.contains(m)).count();
    
    // No hits in either index
    if hits1 == 0 && hits2 == 0 {
        return ReadAssignment::Ambiguous;
    }
    
    // Only hits in index1
    if hits1 > 0 && hits2 == 0 {
        return ReadAssignment::Index1;
    }
    
    // Only hits in index2
    if hits1 == 0 && hits2 > 0 {
        return ReadAssignment::Index2;
    }
    
    // Both indexes have hits - compute ratio (hits1 / hits2)
    let ratio = hits1 as f64 / hits2 as f64;
    let min_ratio = self.config.min_ratio_threshold;
    
    // Check if ratio is sufficiently extreme
    if ratio >= min_ratio {
        return ReadAssignment::Index1;  // Index1 dominates
    } else if ratio <= 1.0 / min_ratio {
        return ReadAssignment::Index2;  // Index2 dominates
    } else {
        return ReadAssignment::Ambiguous;  // Too similar (between 1/min_ratio and min_ratio)
    }
}
```

- **Rewrite** `process()` (single-end only):
  - Extract minimizers via `fill_minimizers()`
  - Call `assign_read(minimizers)`
  - Route to appropriate writer (hap0, hap1, or hap2) based on assignment
  - Update local stats
  - Format record and write to local buffer
  - Flush to global writer on batch complete

- **Update** `FilterSummary`:
```rust
pub struct FilterSummary {
    pub total_reads: usize,
    pub hap0_count: usize,  // ambiguous
    pub hap1_count: usize,  // haplotype 1
    pub hap2_count: usize,  // haplotype 2
}
```
  - Print: "Haplotype 0 (ambiguous): X (Y%), Haplotype 1: A (B%), Haplotype 2: C (D%)"

#### `src/main.rs`
**BREAKING CHANGES**: Simplify CLI for dual-index long-read binning only

- **Remove** `Filter` as top-level command, **replace** with simpler structure:
```
deacon filter \
  --index1 <PATH> \           # required
  --index2 <PATH> \           # required
  --input <PATH> \            # required (long reads)
  --hap0 <PATH> \             # optional (ambiguous reads, omit to discard)
  --hap1 <PATH> \             # required (haplotype 1)
  --hap2 <PATH> \             # required (haplotype 2)
  --min-ratio <FLOAT> \       # default: 2.0 (ratio = hits1/hits2)
  --prefix-length <N> \       # default: 0 (whole sequence)
  --threads <N> \
  --compression-level <N> \
  --rename \                  # sequential numbering
  --quiet \
  --debug
```

- **Remove** flags: `--abs-threshold`, `--rel-threshold`, `--deplete`, `--output`, `--input2`, `--output2`, `--tie-strategy`
- **Validation**:
  - Error if `--index1` or `--index2` missing
  - Error if `--hap1` or `--hap2` missing
  - Error if k/w mismatch between indexes (check headers)
  - Warn if `--hap0` omitted (ambiguous reads will be discarded)
  - Error if `--min-ratio < 1.0` (must be ≥1.0, typically 1.5-3.0)

### Implementation Details
- Both indexes loaded into memory (no lazy loading)
- **Long reads only** - no paired-end processing in filter mode
- Distinct minimizers counted per read
- Output format: standard FASTQ/FASTA (no header modification)
- Ratio threshold applied: `hits1/hits2 >= min_ratio` for hap1, `hits1/hits2 <= 1/min_ratio` for hap2

---

## Workflow Example

### Use Case: Diploid Haplotype Binning (Long Reads)

```bash
# Step 1: Build haplotype 1 index with frequency filtering
deacon index build \
  --input haplotype1_R1.fq.gz \
  --input2 haplotype1_R2.fq.gz \
  --min-frequency 10 \
  -k 31 -w 15 \
  -o hap1.idx

# Step 2: Build haplotype 2 index with frequency filtering
deacon index build \
  --input haplotype2_R1.fq.gz \
  --input2 haplotype2_R2.fq.gz \
  --min-frequency 10 \
  -k 31 -w 15 \
  -o hap2.idx

# Step 3: Create exclusive indexes (remove shared minimizers)
deacon index subtract hap1.idx hap2.idx -o hap1_exclusive.idx
deacon index subtract hap2.idx hap1.idx -o hap2_exclusive.idx

# Step 4: Bin long reads by haplotype with ratio threshold
deacon filter \
  --index1 hap1_exclusive.idx \
  --index2 hap2_exclusive.idx \
  --input long_reads.fq.gz \
  --hap0 ambiguous.fq.gz \
  --hap1 haplotype1_reads.fq.gz \
  --hap2 haplotype2_reads.fq.gz \
  --min-ratio 2.0 \
  --threads 16
```

**Interpretation**:
- `--min-ratio 2.0`: ratio `hits1/hits2` must be ≥2.0 for hap1 or ≤0.5 for hap2
- Reads with ratio between 0.5 and 2.0 → hap0 (ambiguous)
- Reads with no hits → hap0 (ambiguous)
- Long reads only (no paired-end support in filtering)

---

## Testing Strategy

### New Tests to Add

#### `tests/index_paired_build.rs`
- Test paired-end index building
- Verify minimizers from both mates are merged
- Check single-end vs paired-end produces different indexes

#### `tests/index_frequency_filter.rs`
- Build index with `--min-frequency 5`
- Verify rare minimizers removed
- Check histogram output
- Test edge case: frequency threshold higher than max count (empty index)

#### `tests/index_subtract.rs`
- Create two indexes with overlapping minimizers
- Run subtract operation
- Verify result contains only A-B minimizers
- Test k/w mismatch error handling

#### `tests/filter_dual_index.rs`
- **Haplotype 1**: read with hits1=10, hits2=2, ratio=5.0 ≥ 2.0
- **Haplotype 2**: read with hits1=3, hits2=12, ratio=0.25 ≤ 0.5
- **Haplotype 0 (no hits)**: read with hits1=0, hits2=0
- **Haplotype 0 (ratio too close)**: read with hits1=10, hits2=8, ratio=1.25 (between 0.5 and 2.0)
- **Edge case**: read with hits1=5, hits2=0 → haplotype 1 (only index1 has hits)
- **Edge case**: read with hits1=0, hits2=7 → haplotype 2 (only index2 has hits)

#### `tests/filter_output_routing.rs`
- Verify hap0/hap1/hap2 reads written to correct outputs
- Test omitting `--hap0` (ambiguous reads discarded)
- Long reads only (no paired-end test needed)

---

## Data Structure Summary

### New Types in `src/lib.rs`

```rust
// Frequency tracking during index build
pub type RapidHashMap128 = HashMap<u128, u32, FixedRapidHasher>;

pub struct MinimizerFrequencies {
    counts: RapidHashMap128,
}

// Read assignment for dual-index filtering
pub enum ReadAssignment {
    Ambiguous,  // haplotype 0
    Index1,     // haplotype 1
    Index2,     // haplotype 2
}
```

---

## Implementation Order

1. **Feature 2** (frequency filtering) - foundational, needed for quality indexes
2. **Feature 1** (paired-end building) - extends Feature 2 to paired data
3. **Feature 3** (subtract) - enables creating exclusive indexes
4. **Feature 4** (dual-index binning) - breaking changes, implement last

---

## Design Decisions (Resolved)

1. **Paired-end support**: ✅ **Long reads only** - no paired-end input support in filtering, only single-end long reads

2. **Default `--min-ratio`**: ✅ **2.0** - ratio threshold for haplotype assignment

3. **Histogram output**: ✅ **Print all individual values 1-1023, combine ≥1024 into "1024+"** - complete frequency distribution without binning

4. **Ambiguous reads**: ✅ **Warn once at start** if `--hap0` omitted (ambiguous reads will be discarded)

5. **Ratio calculation**: ✅ **hits1/hits2** (not max/min) - if hits2=0 and hits1>0, automatically assign to haplotype 1

6. **Output naming**: ✅ **--hap0, --hap1, --hap2** - three output files for haplotype bins
