Project plan: minimizer counts + two-index compare mode

## STEP 0: Simplify codebase - Remove u64 variant, keep only u128

### Rationale
The current codebase maintains dual code paths for u64 and u128 minimizers, with branching logic throughout. This adds ~500+ lines of duplicated code and complexity. By standardizing on u128 for all k values (even small ones), we significantly simplify the codebase at the cost of using 16 bytes per minimizer instead of 8 bytes for k<=32.

### Files to modify

#### 1) src/minimizers.rs
- **Remove** `decode_u64()` function entirely
- **Keep** `decode_u128()` as the single decode function (rename to `decode` for simplicity)
- **Remove** k<=32 check in `compute_minimizers()` - always use `Buffers::new_u128()`
- **Remove** `U64` match arm in `fill_minimizers()` - only keep `U128` variant
- **Update tests**: Remove u64-specific tests, update remaining to use u128

#### 2) src/lib.rs
- **Remove** `U64(RapidHashSet<u64>)` variant from `MinimizerSet` enum
- **Remove** `U64(Vec<u64>)` variant from `MinimizerVec` enum
- **Simplify** all methods on `MinimizerSet`: remove match arms for U64 variant
- **Remove** `decode_u64` from public exports
- **Remove** `is_u64()` method entirely

#### 3) src/index.rs
- **Remove** k<=32 branching in `load_minimizers()` - always load as u128
- **Remove** k<=32 branching in `dump_minimizers()` - always dump as u128
- **Remove** k<=32 branching in `build()` - always use u128 processor
- **Remove** `local_minimizers_u64` and `global_minimizers_u64` fields
- **Simplify** all match statements to only handle U128 variant

#### 4) src/filter.rs
- **Remove** `decode_u64` from imports
- **Remove** `Buffers::new_u64()` - make `new_u128()` the default
- **Simplify** all match statements to only handle U128 variant

#### 5) Update validation constraints
- k can now go up to 61 (u128 supports up to 64 bp with 2-bit encoding)
- Keep k+w <= 96 and k+w even constraints

### Expected changes summary
- **Lines removed**: ~500-600 (all U64 code paths, match arms, tests)
- **Memory impact**: 2x memory usage for k<=32 (acceptable tradeoff for simplicity)
- **Performance impact**: Negligible (u128 operations are well-optimized)

---

## STEP 1: Add minimizer count tracking

Summary

Add an optional minimizer-counts artifact produced by indexing (a fast HashMap of minimizer -> u32 counts). At filter-time, load the counts companion file and optionally drop rarely-occurring minimizers before building the in-memory MinimizerSet used for read matching. Add a new two-index compare filtering mode so a read (or read-pair) is assigned to the index with more distinct minimizer hits (tie-breaking configurable).

Files to change (file-by-file plan)

1) src/lib.rs
- Add type alias:
  - `pub type RapidHashMap128 = HashMap<u128, u32, FixedRapidHasher>;`
- Add `pub enum MinimizerCounts { U128(RapidHashMap128) }` (only u128 variant after Step 0).
Purpose: centralize count types for use across index & filter.

2) src/index.rs
- Add `dump_minimizer_counts(path: &Path, counts: &MinimizerCounts) -> Result<()>` and `load_minimizer_counts(path: &Path) -> Result<MinimizerCounts>`.
- Optionally add a small `--emit-counts` flag to the index `build` routine that writes the companion counts file `<index>.counts.bin`.
Purpose: persist/load minimizer occurrence counts without changing existing index binary format.

3) src/minimizers.rs
- Add helper `pub fn filter_minimizers_by_count(counts: &MinimizerCounts, min_count: u32) -> MinimizerSet` that builds a MinimizerSet containing only keys with count >= min_count.
Purpose: convert counts map into the MinimizerSet expected by `FilterProcessor`.

4) src/filter.rs
- Extend `FilterConfig` to accept second index path (`index2_path: Option<PathBuf>`) and compare options (`compare_tie` enum, `compare_require_thresholds: bool`, `compare_min_count: Option<u32>`).
- `FilterProcessor` gains an optional `minimizers_b: Option<MinimizerSet>` and logic to load it.
- Add method `should_choose_index(&mut self, seq/mates) -> (Option<usize>, hits1, hits2)` that computes distinct hits for each index, sums across mates, applies tie-breaker and threshold checks, and returns assigned index (1 or 2) or None.
- Update core processing pipeline: when `index2` provided, use `should_choose_index` to decide whether to keep/deplete/route reads and how to update summary.
Purpose: implement two-index compare mode and integrate counts-filtered MinimizerSets.

5) src/main.rs
- Add CLI flags to `filter` subcommand: `--index2 <PATH>`, `--compare-tie {prefer1,prefer2,none}`, `--compare-require-thresholds`, `--compare-min-count <N>`, and pass these into `FilterConfig`.
Purpose: expose the new mode to users.

6) tests/
- Add tests:
  - `tests/index_counts.rs`: building an index with `--emit-counts` produces a counts file and `load_minimizer_counts` returns expected mapping.
  - `tests/filter_two_index.rs`: single read and paired read scenarios where index1 wins, index2 wins, and ties (test tie strategy choices).
  - `tests/filter_min_count.rs`: ensure `--compare-min-count` removes rare minimizers before matching and changes assignment accordingly.
Purpose: cover new code paths and compatibility.

Serialization / file format

- Implement companion counts file as compact binary using serde + bincode. File layout suggestion:
  - Header: format_version (u8), kmer_length (u8), window_size (u8)
  - Number of entries (usize)
  - For each entry: minimizer (u64 or u128 in little-endian) then count (u32 little-endian)
- Filename: if index file is `db.idx` then counts file `db.idx.counts.bin` (or `db.counts`); document choice in README.
- Rationale: bincode is fast and compact; companion file avoids breaking existing index format and tools.

CLI choices & behavior details

- New flags (under `filter`):
  - `--index2 <PATH>` : second index path to compare against.
  - `--compare-min-count <N>` : drop minimizers with global counts < N from both indexes before matching.
  - `--compare-tie {prefer1,prefer2,none}` : default `prefer1`.
  - `--compare-require-thresholds` : when set, winner must additionally meet the existing `abs`/`rel` thresholds to be considered a match.
- Per-read algorithm:
  - Extract distinct minimizers per mate; apply prefix-length.
  - For each index, count distinct minimizers found (deduplicated per-read).
  - Sum across mates for paired reads.
  - Compare hit1 and hit2; winner = argmax; ties resolved according to `compare-tie` (or treated as no-match if `none`).
  - If `compare-require-thresholds` is set, ensure winner satisfies the standard thresholds (abs and relative) using combined total minimizer count.
- Output routing: minimal change path — keep existing single output; annotate summary with assigned index. Optionally later add `--output1`/`--output2` to route reads to separate files.

Backwards compatibility & rollout

- If counts file is missing, filter falls back to current behavior.
- Companion file approach preserves compatibility; embedding counts into index would require version bump and compatibility code — avoid for now.

Edge cases & tests (high-priority)

- Missing counts file: ensure fallback behaves identically.
- Index header mismatch (k/w): error and refuse to compare.
- Tie behavior tests for all `compare-tie` choices.
- u128 indexing: ensure large-k indexes use u128 keys for counts and filtering.

Next steps

- Implement the minimal API additions in `src/lib.rs` and `src/index.rs` to support counts IO.
- Implement `filter_minimizers_by_count` in `src/minimizers.rs` and integrate loading in `src/filter.rs`.
- Add CLI flags in `src/main.rs` and wire to `FilterConfig`.
- Add tests and run `cargo test`.

If you want, I can now implement these edits and run the test suite, starting with the minimal companion-file IO and filter integration. Otherwise I can adjust the plan (embed counts vs companion file; output routing) before coding.

## Notes

Do not worry about backward compatability
Make histogram of minimizer frequency
Only save keys of the hashmap, i.e. a hashset, do not write the counts