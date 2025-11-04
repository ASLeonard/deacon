use crate::index::load_minimizers_cached;
use crate::minimizers::KmerHasher;
use crate::{FilterConfig, MinimizerSet, ReadAssignment};
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use packed_seq::{PackedNSeqVec, SeqVec, u32x8};
use paraseq::Record;
use paraseq::fastx::Reader;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::io::{self, BufWriter, Write};
use std::sync::Arc;
use std::time::Instant;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // 8MB output buffer
const DEFAULT_BUFFER_SIZE: usize = 64 * 1024;

type BoxedWriter = Box<dyn Write + Send>;

/// Check input file paths exist
fn check_input_paths(config: &FilterConfig) -> Result<()> {
    if !config.index1_path.exists() {
        return Err(anyhow::anyhow!(
            "Index 1 file does not exist: {}",
            config.index1_path.display()
        ));
    }
    
    if !config.index2_path.exists() {
        return Err(anyhow::anyhow!(
            "Index 2 file does not exist: {}",
            config.index2_path.display()
        ));
    }

    // Skip stdin case
    if config.input_path != "-" && !std::path::Path::new(config.input_path).exists() {
        return Err(anyhow::anyhow!(
            "Input file does not exist: {}",
            config.input_path
        ));
    }

    Ok(())
}

/// Check if file metadata len < 5 (catches empty uncompressed files only)
fn is_empty_file(path: &str) -> Result<bool> {
    if path == "-" {
        return Ok(false); // Can't check stdin
    }
    let metadata = std::fs::metadata(path)
        .map_err(|e| anyhow::anyhow!("Failed to read file metadata {}: {}", path, e))?;
    Ok(metadata.len() < 5)
}

/// Error thrown when paraseq reads empty compressed file
fn is_empty_input_error(err: &anyhow::Error) -> bool {
    err.to_string().contains("failed to fill whole buffer")
}

/// Create a paraseq reader from optional path (stdin if None or "-")
fn create_paraseq_reader(path: Option<&str>) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    match path {
        None | Some("-") => {
            let stdin_reader = Box::new(std::io::stdin()) as Box<dyn std::io::Read + Send>;
            Reader::new(stdin_reader)
                .map_err(|e| anyhow::anyhow!("Failed to create stdin reader: {}", e))
        }
        Some(p) => {
            Reader::from_path(p).map_err(|e| anyhow::anyhow!("Failed to open file {}: {}", p, e))
        }
    }
}

/// Format a single record into a buffer (FASTA/FASTQ format)
/// `seq` is the newline-stripped sequence corresponding to the record from `record.seq()`.
fn format_record_to_buffer<R: Record>(
    record: &R,
    seq: &[u8],
    counter: u64,
    rename: bool,
    buffer: &mut Vec<u8>,
) -> Result<()> {
    let is_fasta = record.qual().is_none();

    // Write header line (> or @ prefix)
    buffer.push(if is_fasta { b'>' } else { b'@' });

    // Write sequence identifier
    if rename {
        buffer.extend_from_slice(counter.to_string().as_bytes());
    } else {
        buffer.extend_from_slice(record.id());
    }
    buffer.push(b'\n');

    // Write sequence
    buffer.extend_from_slice(seq);
    buffer.push(b'\n');

    // FASTQ: write + and quality
    if !is_fasta {
        buffer.extend_from_slice(b"+\n");
        buffer.extend_from_slice(record.qual().unwrap());
        buffer.push(b'\n');
    }

    Ok(())
}

/// Open a boxed writer for output path (stdout if None)
fn open_output_writer(path: Option<&std::path::Path>, compression_level: u8) -> Result<BoxedWriter> {
    match path {
        None => Ok(Box::new(BufWriter::with_capacity(
            OUTPUT_BUFFER_SIZE,
            io::stdout(),
        ))),
        Some(p) => {
            let file = OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(p)
                .with_context(|| format!("Failed to create output file: {}", p.display()))?;

            let writer: BoxedWriter = if let Some(ext) = p.extension() {
                match ext.to_str() {
                    Some("gz") => {
                        let encoder = flate2::write::GzEncoder::new(
                            file,
                            flate2::Compression::new(compression_level as u32),
                        );
                        Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, encoder))
                    }
                    Some("zst") | Some("zstd") => {
                        let encoder = zstd::Encoder::new(file, compression_level as i32)?
                            .auto_finish();
                        Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, encoder))
                    }
                    _ => Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file)),
                }
            } else {
                Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file))
            };

            Ok(writer)
        }
    }
}

/// Summary struct for JSON output
#[derive(Serialize, Deserialize, Debug)]
pub struct FilterSummary {
    version: String,
    index1: String,
    index2: String,
    input: String,
    hap0_output: Option<String>,
    hap1_output: String,
    hap2_output: String,
    k: u8,
    w: u8,
    min_ratio_threshold: f64,
    prefix_length: usize,
    rename: bool,
    seqs_in: u64,
    hap0_count: u64,
    hap1_count: u64,
    hap2_count: u64,
    bp_in: u64,
    bp_hap0: u64,
    bp_hap1: u64,
    bp_hap2: u64,
    elapsed_seconds: f64,
}

/// Processing statistics tracked per thread and globally
#[derive(Clone, Default, Debug)]
pub(crate) struct ProcessingStats {
    pub total_seqs: u64,
    pub hap0_count: u64,
    pub hap1_count: u64,
    pub hap2_count: u64,
    pub total_bp: u64,
    pub bp_hap0: u64,
    pub bp_hap1: u64,
    pub bp_hap2: u64,
    pub output_hap0_counter: u64,
    pub output_hap1_counter: u64,
    pub output_hap2_counter: u64,
    pub last_reported: u64,
}

impl ProcessingStats {
    fn merge_from(&mut self, other: &Self) {
        self.total_seqs += other.total_seqs;
        self.hap0_count += other.hap0_count;
        self.hap1_count += other.hap1_count;
        self.hap2_count += other.hap2_count;
        self.total_bp += other.total_bp;
        self.bp_hap0 += other.bp_hap0;
        self.bp_hap1 += other.bp_hap1;
        self.bp_hap2 += other.bp_hap2;
    }
}

/// Buffers for minimizer computation
#[derive(Clone)]
pub(crate) struct Buffers {
    pub packed_nseq: PackedNSeqVec,
    pub positions: Vec<u32>,
    pub minimizers: crate::MinimizerVec,
    pub cache: (simd_minimizers::Cache, Vec<u32x8>, Vec<u32x8>),
}

impl Buffers {
    pub fn new_u64() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U64(Vec::new()),
            cache: Default::default(),
        }
    }

    pub fn new_u128() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U128(Vec::new()),
            cache: Default::default(),
        }
    }
}

/// Main processor for dual-index competitive binning
#[derive(Clone)]
struct FilterProcessor {
    // Two minimizer indexes
    minimizers_a: &'static MinimizerSet,
    minimizers_b: &'static MinimizerSet,
    
    kmer_length: u8,
    window_size: u8,
    min_ratio_threshold: f64,
    prefix_length: usize,
    rename: bool,
    debug: bool,

    hasher: KmerHasher,

    // Local buffers (per thread)
    local_buffer_hap0: Vec<u8>,
    local_buffer_hap1: Vec<u8>,
    local_buffer_hap2: Vec<u8>,
    local_stats: ProcessingStats,
    buffers: Buffers,

    // Global shared state
    global_writer_hap0: Option<Arc<Mutex<BoxedWriter>>>,
    global_writer_hap1: Arc<Mutex<BoxedWriter>>,
    global_writer_hap2: Arc<Mutex<BoxedWriter>>,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    filtering_start_time: Instant,
}

impl FilterProcessor {
    /// Flush local buffers to global writers
    fn flush_buffers(&mut self) -> Result<()> {
        // Flush hap0
        if let Some(ref writer) = self.global_writer_hap0 {
            if !self.local_buffer_hap0.is_empty() {
                let mut w = writer.lock();
                w.write_all(&self.local_buffer_hap0)?;
                self.local_buffer_hap0.clear();
            }
        }

        // Flush hap1
        if !self.local_buffer_hap1.is_empty() {
            let mut w = self.global_writer_hap1.lock();
            w.write_all(&self.local_buffer_hap1)?;
            self.local_buffer_hap1.clear();
        }

        // Flush hap2
        if !self.local_buffer_hap2.is_empty() {
            let mut w = self.global_writer_hap2.lock();
            w.write_all(&self.local_buffer_hap2)?;
            self.local_buffer_hap2.clear();
        }

        Ok(())
    }

    /// Update progress spinner if present
    fn update_spinner(&mut self) {
        if let Some(ref spinner) = self.spinner {
            if self.local_stats.total_seqs >= self.local_stats.last_reported + 10000 {
                let mut stats = self.global_stats.lock();
                stats.merge_from(&self.local_stats);
                
                spinner.lock().set_message(format!(
                    "{} reads | hap0: {} | hap1: {} | hap2: {}",
                    stats.total_seqs,
                    stats.hap0_count,
                    stats.hap1_count,
                    stats.hap2_count
                ));

                self.local_stats = ProcessingStats::default();
                self.local_stats.last_reported = stats.total_seqs;
            }
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for FilterProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        // Convert sequence to packed format
        let seq_bytes = record.seq();
        let seq_len = seq_bytes.len();
        
        // Skip if shorter than kmer length
        if seq_len < self.kmer_length as usize {
            self.local_stats.total_seqs += 1;
            self.local_stats.total_bp += seq_len as u64;
            self.local_stats.hap0_count += 1;
            self.local_stats.bp_hap0 += seq_len as u64;
            return Ok(());
        }

        // Clear packed sequence buffers
        self.buffers.packed_nseq.seq.clear();
        self.buffers.packed_nseq.ambiguous.clear();
        
        // Push sequence into packed format
        self.buffers.packed_nseq.ambiguous.push_ascii(&seq_bytes);
        
        // Apply prefix length if specified
        let effective_len = if self.prefix_length > 0 {
            self.prefix_length.min(seq_len)
        } else {
            seq_len
        };

        let seq_bp = effective_len;

        // Compute minimizers and assign
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;
        let m = simd_minimizers::canonical_minimizers(k, w)
            .hasher(&self.hasher)
            .run_skip_ambiguous_windows_with_buf(
                self.buffers.packed_nseq.as_slice(),
                &mut self.buffers.positions,
                &mut self.buffers.cache,
            );

        // Populate minimizer buffer and compute assignment
        let assignment = match (self.minimizers_a, self.minimizers_b, &mut self.buffers.minimizers) {
            (MinimizerSet::U64(set_a), MinimizerSet::U64(set_b), crate::MinimizerVec::U64(mins)) => {
                mins.clear();
                mins.extend(m.pos_and_values_u64().map(|(_pos, val)| val));
                
                // Compute assignment
                let hits1 = mins.iter().filter(|&&m| set_a.contains(&m)).count();
                let hits2 = mins.iter().filter(|&&m| set_b.contains(&m)).count();
                
                if hits1 == 0 && hits2 == 0 {
                    ReadAssignment::Ambiguous
                } else if hits2 == 0 {
                    ReadAssignment::Index1
                } else if hits1 == 0 {
                    ReadAssignment::Index2
                } else {
                    let ratio = hits1 as f64 / hits2 as f64;
                    if ratio >= self.min_ratio_threshold {
                        ReadAssignment::Index1
                    } else if ratio <= 1.0 / self.min_ratio_threshold {
                        ReadAssignment::Index2
                    } else {
                        ReadAssignment::Ambiguous
                    }
                }
            }
            (MinimizerSet::U128(set_a), MinimizerSet::U128(set_b), crate::MinimizerVec::U128(mins)) => {
                mins.clear();
                mins.extend(m.pos_and_values_u128().map(|(_pos, val)| val));
                
                // Compute assignment
                let hits1 = mins.iter().filter(|&&m| set_a.contains(&m)).count();
                let hits2 = mins.iter().filter(|&&m| set_b.contains(&m)).count();
                
                if hits1 == 0 && hits2 == 0 {
                    ReadAssignment::Ambiguous
                } else if hits2 == 0 {
                    ReadAssignment::Index1
                } else if hits1 == 0 {
                    ReadAssignment::Index2
                } else {
                    let ratio = hits1 as f64 / hits2 as f64;
                    if ratio >= self.min_ratio_threshold {
                        ReadAssignment::Index1
                    } else if ratio <= 1.0 / self.min_ratio_threshold {
                        ReadAssignment::Index2
                    } else {
                        ReadAssignment::Ambiguous
                    }
                }
            }
            _ => {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Minimizer type mismatch between indexes").into());
            }
        };

        // Update stats and route to appropriate output
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_bp as u64;

        match assignment {
            ReadAssignment::Ambiguous => {
                self.local_stats.hap0_count += 1;
                self.local_stats.bp_hap0 += seq_bp as u64;
                
                if let Some(_) = self.global_writer_hap0 {
                    self.local_stats.output_hap0_counter += 1;
                    format_record_to_buffer(
                        &record,
                        &seq_bytes,
                        self.local_stats.output_hap0_counter,
                        self.rename,
                        &mut self.local_buffer_hap0,
                    )?;
                }
            }
            ReadAssignment::Index1 => {
                self.local_stats.hap1_count += 1;
                self.local_stats.bp_hap1 += seq_bp as u64;
                self.local_stats.output_hap1_counter += 1;
                
                format_record_to_buffer(
                    &record,
                    &seq_bytes,
                    self.local_stats.output_hap1_counter,
                    self.rename,
                    &mut self.local_buffer_hap1,
                )?;
            }
            ReadAssignment::Index2 => {
                self.local_stats.hap2_count += 1;
                self.local_stats.bp_hap2 += seq_bp as u64;
                self.local_stats.output_hap2_counter += 1;
                
                format_record_to_buffer(
                    &record,
                    &seq_bytes,
                    self.local_stats.output_hap2_counter,
                    self.rename,
                    &mut self.local_buffer_hap2,
                )?;
            }
        }

        // Flush buffers if they get large
        if self.local_buffer_hap0.len() > OUTPUT_BUFFER_SIZE
            || self.local_buffer_hap1.len() > OUTPUT_BUFFER_SIZE
            || self.local_buffer_hap2.len() > OUTPUT_BUFFER_SIZE
        {
            self.flush_buffers()?;
        }

        self.update_spinner();

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Flush local buffers
        self.flush_buffers()?;

        // Merge local stats into global
        {
            let mut stats = self.global_stats.lock();
            stats.merge_from(&self.local_stats);
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

/// Main entry point for dual-index competitive binning
pub fn run(config: &FilterConfig) -> Result<()> {
    let start_time = Instant::now();

    // Validate input paths
    check_input_paths(config)?;

    // Check for empty input files
    if is_empty_file(config.input_path)? {
        eprintln!("Warning: Input file appears to be empty");
        return Ok(());
    }

    // Load both minimizer indexes
    let (minimizers_a, header_a) = load_minimizers_cached(config.index1_path)?;
    let (minimizers_b, header_b) = load_minimizers_cached(config.index2_path)?;

    // Verify indexes have matching parameters
    if header_a.kmer_length != header_b.kmer_length || header_a.window_size != header_b.window_size {
        return Err(anyhow::anyhow!(
            "Index parameters do not match: index1 has k={}, w={}, index2 has k={}, w={}",
            header_a.kmer_length, header_a.window_size, header_b.kmer_length, header_b.window_size
        ));
    }

    let kmer_length = header_a.kmer_length;
    let window_size = header_a.window_size;

    // Open output writers
    let writer_hap0 = if let Some(path) = config.hap0_path {
        Some(Arc::new(Mutex::new(open_output_writer(
            Some(path),
            config.compression_level,
        )?)))
    } else {
        None
    };

    let writer_hap1 = Arc::new(Mutex::new(open_output_writer(
        Some(config.hap1_path),
        config.compression_level,
    )?));

    let writer_hap2 = Arc::new(Mutex::new(open_output_writer(
        Some(config.hap2_path),
        config.compression_level,
    )?));

    // Create progress spinner if not quiet
    let spinner = if !config.quiet {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        pb.set_draw_target(ProgressDrawTarget::stderr());
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    // Initialize global stats
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));

    // Determine if indexes use u64 or u128
    let use_u128 = match minimizers_a {
        MinimizerSet::U64(_) => false,
        MinimizerSet::U128(_) => true,
    };

    // Create reader
    let reader = create_paraseq_reader(Some(config.input_path))?;

    // Create processor
    let mut processor = FilterProcessor {
        minimizers_a,
        minimizers_b,
        kmer_length,
        window_size,
        min_ratio_threshold: config.min_ratio_threshold,
        prefix_length: config.prefix_length,
        rename: config.rename,
        debug: config.debug,
        hasher: KmerHasher::new(kmer_length as usize),
        local_buffer_hap0: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
        local_buffer_hap1: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
        local_buffer_hap2: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
        local_stats: ProcessingStats::default(),
        buffers: if use_u128 {
            Buffers::new_u128()
        } else {
            Buffers::new_u64()
        },
        global_writer_hap0: writer_hap0.clone(),
        global_writer_hap1: writer_hap1.clone(),
        global_writer_hap2: writer_hap2.clone(),
        global_stats: global_stats.clone(),
        spinner: spinner.clone(),
        filtering_start_time: start_time,
    };

    // Determine number of threads
    let num_threads = if config.threads == 0 {
        rayon::current_num_threads()
    } else {
        config.threads
    };

    // Process in parallel
    let result = reader.process_parallel(&mut processor, num_threads);

    // Process result
    result.map_err(|e| anyhow::anyhow!("Processing error: {}", e))?;

    // Finalize writers
    if let Some(w) = writer_hap0 {
        Arc::try_unwrap(w)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap hap0 writer"))?
            .into_inner()
            .flush()?;
    }

    Arc::try_unwrap(writer_hap1)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap hap1 writer"))?
        .into_inner()
        .flush()?;

    Arc::try_unwrap(writer_hap2)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap hap2 writer"))?
        .into_inner()
        .flush()?;

    // Clear spinner
    if let Some(sp) = spinner {
        Arc::try_unwrap(sp)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap spinner"))?
            .into_inner()
            .finish_and_clear();
    }

    let elapsed = start_time.elapsed();
    let final_stats = Arc::try_unwrap(global_stats)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap global stats"))?
        .into_inner();

    if !config.quiet {
        eprintln!("\nFiltering complete!");
        eprintln!("  Total reads:  {}", final_stats.total_seqs);
        eprintln!("  Ambiguous:    {} ({:.2}%)", 
            final_stats.hap0_count,
            100.0 * final_stats.hap0_count as f64 / final_stats.total_seqs.max(1) as f64);
        eprintln!("  Haplotype 1:  {} ({:.2}%)", 
            final_stats.hap1_count,
            100.0 * final_stats.hap1_count as f64 / final_stats.total_seqs.max(1) as f64);
        eprintln!("  Haplotype 2:  {} ({:.2}%)", 
            final_stats.hap2_count,
            100.0 * final_stats.hap2_count as f64 / final_stats.total_seqs.max(1) as f64);
        eprintln!("  Time elapsed: {:.2}s", elapsed.as_secs_f64());
    }

    Ok(())
}
