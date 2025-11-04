use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use deacon::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, FilterConfig, IndexConfig, diff_index, dump_index,
    index_info, intersect_index, subtract_index, union_index,
};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::PathBuf;

#[derive(Parser, Serialize, Deserialize)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    #[arg(long)]
    /// Execute command using an existing server process
    use_server: bool,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum Commands {
    /// Build, compose and inspect minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Retain or deplete sequence records with sufficient minimizer hits to an indexed query
    Filter {
        /// Path to first minimizer index file (haplotype 1)
        #[arg(long = "index1", required = true)]
        index1: PathBuf,

        /// Path to second minimizer index file (haplotype 2)
        #[arg(long = "index2", required = true)]
        index2: PathBuf,

        /// Path to input fastx file (long reads only, or - for stdin)
        #[arg(short = 'i', long = "input", default_value = "-")]
        input: String,

        /// Path to output file for ambiguous reads (haplotype 0) - optional
        #[arg(long = "hap0")]
        hap0: Option<PathBuf>,

        /// Path to output file for haplotype 1 reads
        #[arg(long = "hap1", required = true)]
        hap1: PathBuf,

        /// Path to output file for haplotype 2 reads
        #[arg(long = "hap2", required = true)]
        hap2: PathBuf,

        /// Minimum ratio threshold for haplotype assignment (hits1/hits2 >= min_ratio for hap1)
        #[arg(long = "min-ratio", default_value_t = 2.0)]
        min_ratio: f64,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'R', long = "rename", default_value_t = false)]
        rename: bool,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Output compression level (1-9 for gz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 3)]
        compression_level: u8,

        /// Output sequences with minimizer hits to stderr
        #[arg(long = "debug", default_value_t = false)]
        debug: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
    /// Start/stop a server process for reduced latency filtering
    Server {
        #[command(subcommand)]
        command: ServerCommands,
    },
    /// Show citation information
    Cite,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum ServerCommands {
    /// Start the server
    Start {
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,
    },
    /// Stop the running server
    Stop,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum IndexCommands {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (supports gz, zst and xz compression)
        input: PathBuf,

        /// Path to optional second paired fastx file (for paired-end reads)
        #[arg(long = "input2")]
        input2: Option<PathBuf>,

        /// K-mer length used for indexing (k+w-1 must be <= 96 and odd)
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: u8,

        /// Minimizer window size used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Suppress sequence header output
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,

        /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
        #[arg(short = 'e', long = "entropy-threshold", default_value = "0.0")]
        entropy_threshold: f32,

        /// Minimum frequency threshold for minimizer filtering (only retain minimizers occurring at least N times)
        #[arg(short = 'm', long = "min-frequency")]
        min_frequency: Option<u32>,
    },
    /// Combine multiple minimizer indexes (A ∪ B…)
    Union {
        /// Path(s) to one or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Intersect multiple minimizer indexes (A ∩ B…)
    Intersect {
        /// Path(s) to two or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Subtract minimizers in index B from index A (A - B), creating exclusive minimizers
    Subtract {
        /// Path to first index file (A)
        #[arg(required = true)]
        index_a: PathBuf,

        /// Path to second index file to subtract (B)
        #[arg(required = true)]
        index_b: PathBuf,

        /// Path to output file
        #[arg(short = 'o', long = "output", required = true)]
        output: PathBuf,
    },
    /// Subtract minimizers from FASTX or index (legacy, use Subtract for index-index operations)
    Diff {
        /// Path to first index file
        #[arg(required = true)]
        first: PathBuf,

        /// Path to second index file or FASTX file (or - for stdin when using FASTX)
        #[arg(required = true)]
        second: PathBuf,

        /// K-mer length (required if second argument is FASTX file, 1-32)
        #[arg(short = 'k', long = "kmer-length", value_parser = clap::value_parser!(u8).range(1..=32))]
        kmer_length: Option<u8>,

        /// Window size (required if second argument is FASTX file)
        #[arg(short = 'w', long = "window-size")]
        window_size: Option<u8>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Dump minimizer index to fasta
    Dump {
        /// Path to index file
        index: PathBuf,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Show index information
    Info {
        /// Path to index file
        index: PathBuf,
    },
}

#[derive(Serialize, Deserialize)]
enum Message {
    /// client -> server
    Command(Commands),
    /// server -> client
    Done,
}

fn print_citation() {
    println!("Bede Constantinides, John Lees, Derrick W Crook.");
    println!("\"Deacon: fast sequence filtering and contaminant depletion\"");
    println!("bioRxiv 2025.06.09.658732");
    println!("https://doi.org/10.1101/2025.06.09.658732");
}

fn main() -> Result<()> {
    // Check we have either AVX2 or NEON
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    let cli = Cli::parse();

    if let Commands::Server { command } = &cli.command {
        if !cli.use_server {
            // Running server commands directly (not via --use-server)
            match command {
                ServerCommands::Start { threads } => {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(*threads)
                        .build_global()
                        .context("Failed to initialize thread pool")?;

                    // Remove existing socket if present
                    let _ = std::fs::remove_file("deacon_server_socket");
                    let listener = UnixListener::bind("deacon_server_socket")?;
                    for stream in listener.incoming() {
                        let mut stream = stream.unwrap();
                        let mut message = vec![];
                        let mut buf = vec![0; 10000];
                        loop {
                            let len = stream.read(&mut buf)?;
                            let buf = &buf[..len];
                            message.extend_from_slice(buf);
                            if buf.contains(&0) {
                                assert_eq!(buf.last(), Some(&0));
                                message.pop();
                                break;
                            }
                        }
                        let message: Message = serde_json::from_slice(&message).unwrap();
                        match message {
                            Message::Command(Commands::Server {
                                command: ServerCommands::Stop,
                            }) => {
                                serde_json::to_writer(stream, &Message::Done)?;
                                let _ = std::fs::remove_file("deacon_server_socket");
                                break;
                            }
                            Message::Command(commands) => {
                                process_command(&commands)?;
                                serde_json::to_writer(stream, &Message::Done)?;
                            }
                            Message::Done => {
                                unreachable!("Server should not receive `Done` messages.")
                            }
                        }
                    }

                    return Ok(());
                }
                ServerCommands::Stop => {
                    panic!("Use `deacon --use-server server stop` to stop the server.")
                }
            }
        }
        // If --use-server is set, fall through to client code below
    }

    if cli.use_server {
        let mut stream = UnixStream::connect("deacon_server_socket")?;
        serde_json::to_writer(&stream, &Message::Command(cli.command))?;
        stream.write(b"\0")?;
        stream.flush()?;
        let message: Message = serde_json::from_reader(stream).unwrap();
        match message {
            Message::Done => {}
            _ => unreachable!("The client only expects to receive `Done` messages."),
        }

        return Ok(());
    }

    process_command(&cli.command)?;

    Ok(())
}

fn process_command(command: &Commands) -> Result<(), anyhow::Error> {
    match &command {
        Commands::Cite => {
            print_citation();
        }
        Commands::Server { .. } => {
            unreachable!("Server commands are handled before this function is called")
        }
        Commands::Index { command } => match command {
            IndexCommands::Build {
                input,
                input2,
                kmer_length,
                window_size,
                output,
                threads,
                quiet,
                entropy_threshold,
                min_frequency,
            } => {
                let config = IndexConfig {
                    input_path: input.clone(),
                    input2_path: input2.clone(),
                    kmer_length: *kmer_length,
                    window_size: *window_size,
                    output_path: output.clone(),
                    threads: *threads,
                    quiet: *quiet,
                    entropy_threshold: *entropy_threshold,
                    min_frequency: *min_frequency,
                };
                config
                    .execute()
                    .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            IndexCommands::Union { inputs, output } => {
                union_index(inputs, output.as_deref())
                    .context("Failed to run index union command")?;
            }
            IndexCommands::Intersect { inputs, output } => {
                intersect_index(inputs, output.as_deref())
                    .context("Failed to run index intersect command")?;
            }
            IndexCommands::Subtract {
                index_a,
                index_b,
                output,
            } => {
                subtract_index(index_a, index_b, output)
                    .context("Failed to run index subtract command")?;
            }
            IndexCommands::Diff {
                first,
                second,
                kmer_length,
                window_size,
                output,
                threads,
            } => {
                diff_index(
                    first,
                    second,
                    *kmer_length,
                    *window_size,
                    *threads,
                    output.as_deref(),
                )
                .context("Failed to run index diff command")?;
            }
            IndexCommands::Dump { index, output } => {
                dump_index(index, output.as_deref()).context("Failed to run index dump command")?;
            }
        },
        Commands::Filter {
            index1,
            index2,
            input,
            hap0,
            hap1,
            hap2,
            min_ratio,
            prefix_length,
            rename,
            threads,
            compression_level,
            debug,
            quiet,
        } => {
            let config = FilterConfig {
                index1_path: index1,
                index2_path: index2,
                input_path: input,
                hap0_path: hap0.as_deref(),
                hap1_path: hap1,
                hap2_path: hap2,
                min_ratio_threshold: *min_ratio,
                prefix_length: *prefix_length,
                rename: *rename,
                threads: *threads,
                compression_level: *compression_level,
                debug: *debug,
                quiet: *quiet,
            };
            config.execute().context("Failed to run filter command")?;
        }
    }

    Ok(())
}
