use bytemuck::cast_slice_mut;
use parse_int::parse;
use rand::{thread_rng, Rng};
use std::fs::File;
use std::io::{BufRead, Error, ErrorKind, Read, Result, Stdin};
use std::num::Wrapping as W;
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;
use structopt::StructOpt;

/// XSH-RR output transformation from PCG.
fn compute_xsh_rr(state: W<u64>) -> u32 {
    let xorshifted = ((state >> 18) ^ state) >> 27;
    let rotation = (state >> 59).0 as u32;

    (xorshifted.0 as u32).rotate_right(rotation)
}

/// Inverts XSH-RR given one particular rotation guess.
fn invert_xsh_rr(rotation: u32, output: u32) -> W<u64> {
    let mut state = W(rotation as u64) << 59;

    let recovered = W(output.rotate_left(rotation) as u64);

    state |= (recovered >> 19) << 46;

    state |= (((recovered >> 1) ^ (state >> 46)) & W(0x3ffff)) << 28;

    state |= ((recovered ^ (state >> 45)) & W(1)) << 27;

    state
}

/// The multiplicative constant from PCG-XSH-RR
const C: W<u64> = W(6_364_136_223_846_793_005);

fn parse_output(line: &str) -> Result<u32> {
    parse::<u32>(line).map_err(|err| Error::new(ErrorKind::Other, err))
}

fn read_output(stdin: &Stdin) -> Result<u32> {
    let mut line = String::new();
    stdin.read_line(&mut line)?;
    parse_output(&line)
}

struct CandidateState {
    epsilon_1: W<u64>,
    increment: W<u64>,

    state: W<u64>,
}

struct LookupTable {
    table: Box<[u64]>,
}

impl LookupTable {
    pub fn open(path: &Path) -> Result<Self> {
        let mut raw_table_file = File::open(path)?;

        let mut table = vec![0; 0x800_0000].into_boxed_slice();
        raw_table_file.read_exact(cast_slice_mut(&mut table))?;

        Ok(Self { table })
    }

    /// Returns the (zeta, beta) solution for N when one exists.
    pub fn query(&self, n: W<u64>) -> Option<(W<u64>, W<u64>)> {
        if let Some((zeta, beta)) = self.check(n) {
            return Some((zeta, beta));
        }

        if let Some((zeta, beta)) = self.check(-n) {
            return Some((-zeta, -beta));
        }

        if let Some((zeta, beta)) = self.check(n - W(1)) {
            return Some((zeta, beta - W(1 << 27)));
        }

        if let Some((zeta, beta)) = self.check(-n - W(1)) {
            return Some((-zeta, W(1 << 27) - beta));
        }

        None
    }

    fn check(&self, n: W<u64>) -> Option<(W<u64>, W<u64>)> {
        if let Some(zeta) = self.binary_search(n) {
            Some((zeta, zeta * C - (n << 27)))
        } else {
            None
        }
    }

    fn binary_search(&self, n: W<u64>) -> Option<W<u64>> {
        let index = match self.table.binary_search(&((n << 27).0)) {
            Ok(index) | Err(index) => index,
        };

        // IMPORTANT: we must mask out the highest 27 bits of N here as adjusting
        // the value of N to locate symmetric solutions might result in N > 2^37.

        if self.table[index] >> 27 == n.0 & 0x1f_ffff_ffff {
            Some(W(self.table[index] & 0x7ff_ffff))
        } else {
            None
        }
    }
}

#[derive(Debug)]
struct PartialStateTuple {
    s0_high_bits: W<u64>,
    s1_high_bits: W<u64>,

    zeta: W<u64>,
    beta: W<u64>,
}

#[derive(Debug, Default)]
struct OutputPredictions {
    frequencies: Vec<(u32, i32)>,
}

impl OutputPredictions {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn record_prediction(&mut self, prediction: u32) {
        for (stored, frequency) in &mut self.frequencies {
            if *stored == prediction {
                *frequency += 1;
                return;
            }
        }

        self.frequencies.push((prediction, 1));
    }

    pub fn clear_predictions(&mut self) {
        self.frequencies.clear();
    }

    pub fn display_predictions(&mut self) {
        // Sort predicted outputs by decreasing probability
        self.frequencies.sort_unstable_by_key(|&(_, f)| -f);

        let mut count = 0.0;

        for &(_, frequency) in &self.frequencies {
            count += frequency as f64;
        }

        println!();

        for &(output, frequency) in &self.frequencies {
            println!(
                "    0x{:08X} ({:.2}% probability)",
                output,
                frequency as f64 / count * 100.0
            );
        }

        println!();
    }

    pub fn unique_prediction(&self) -> Option<u32> {
        if self.frequencies.len() == 1 {
            Some(self.frequencies[0].0)
        } else {
            None
        }
    }
}

fn recover_partial_state(
    table: &LookupTable,
    outputs: &[u32; 4],
) -> Result<(PartialStateTuple, usize)> {
    // In theory, we need to check all 2^27 possible internal states to see if any of them
    // generate output #4 correctly. However, in practice there are usually only 4 or less
    // different outputs which will be produced, so we can just try a few random states.

    // If we fail, just try increasingly many random states and eventually try them all if
    // we are really unlucky. That means we are fast but never produce false negatives.

    const STRATEGIES: &[Strategy] = &[
        Strategy::Random { iterations: 16 },
        Strategy::Random { iterations: 256 },
        Strategy::ExhaustiveSearch,
    ];

    for strategy in STRATEGIES {
        let mut guesses = 0;

        for s0_rotation in 0..32 {
            for s1_rotation in 0..32 {
                for s2_rotation in 0..32 {
                    // Uncover the 37 highest bits of all three states

                    let s0_high_bits = invert_xsh_rr(s0_rotation, outputs[0]);
                    let s1_high_bits = invert_xsh_rr(s1_rotation, outputs[1]);
                    let s2_high_bits = invert_xsh_rr(s2_rotation, outputs[2]);

                    // Compute the star expressions and the LHS

                    let s_star_0 = s0_high_bits >> 27;
                    let s_star_1 = s1_high_bits >> 27;
                    let s_star_2 = s2_high_bits >> 27;

                    let n = (s_star_1 - s_star_2 + C * (s_star_1 - s_star_0)) & W(0x1f_ffff_ffff);

                    // NOTE: we don't actually need to calculate beta since it is predetermined
                    // by the value of zeta, but we compute it anyway to make the code clearer.

                    if let Some((zeta, beta)) = table.query(n) {
                        guesses += 1;

                        if verify_state(s0_high_bits, s1_high_bits, zeta, outputs[3], *strategy) {
                            return Ok((
                                PartialStateTuple {
                                    s0_high_bits,
                                    s1_high_bits,
                                    zeta,
                                    beta,
                                },
                                guesses,
                            ));
                        }
                    }
                }
            }
        }
    }

    Err(Error::new(
        ErrorKind::Other,
        "output sequence not produced by PCG-XSH-RR",
    ))
}

#[derive(Clone, Copy)]
enum Strategy {
    Random { iterations: u32 },
    ExhaustiveSearch,
}

fn verify_state(
    s0_high_bits: W<u64>,
    s1_high_bits: W<u64>,
    zeta: W<u64>,
    x3: u32,
    strategy: Strategy,
) -> bool {
    let mut rng = thread_rng();

    let (iterations, random_search) = match strategy {
        Strategy::Random { iterations } => (iterations, true),
        Strategy::ExhaustiveSearch => (1 << 27, false),
    };

    for i in 0..iterations {
        let epsilon_1 = if random_search {
            W(rng.gen::<u64>() & 0x7ff_ffff)
        } else {
            W(i as u64)
        };

        let mut state = s0_high_bits | (zeta + epsilon_1);
        let increment = ((s1_high_bits | epsilon_1) - state * C) | W(1);

        for _ in 0..3 {
            state = state * C + increment;
        }

        if compute_xsh_rr(state) == x3 {
            return true;
        }
    }

    false
}

fn display_recovered_state(state: W<u64>, increment: W<u64>) {
    println!("\n    pcg32_random_t state = {{");
    println!("        .state = 0x{:016X}", state);
    println!("        .inc   = 0x{:016X}", increment);
    println!("    }};\n");
}

fn run(args: Opt) -> Result<()> {
    println!("{}", ASCII_HEADER);
    let stdin = std::io::stdin();

    println!("[-] Starting clock.");
    let start_time = Instant::now();

    let table = LookupTable::open(&args.table).map_err(|err| {
        println!("[!] Failed to load precomputed table!");
        err
    })?;

    println!("[+] Loaded precomputed table.");

    println!("[-] Reading 4 outputs to recover target (zeta, beta) parameters.");

    let (state, guesses) = recover_partial_state(
        &table,
        &[
            read_output(&stdin)?,
            read_output(&stdin)?,
            read_output(&stdin)?,
            read_output(&stdin)?,
        ],
    )?;

    // We only needed the table for the early filtering phase (everything
    // beyond that is just standard LCG recovery) so reclaim some memory.

    drop(table);

    println!(
        "[+] Target parameters recovered after {:.2} seconds and {} guess(es).",
        start_time.elapsed().as_secs_f64(),
        guesses
    );

    println!(
        "\n     (zeta, beta) = ({}, {})\n",
        state.zeta.0 as i32, state.beta.0 as i32
    );

    println!("[-] Enumerating internal states...");
    let mut predictions = OutputPredictions::new();

    let mut candidates = vec![];

    for i in 0..(1u64 << 27) {
        let epsilon_1 = W(i);

        let s0_high_bits = state.s0_high_bits;
        let s1_high_bits = state.s1_high_bits;

        let mut state = s0_high_bits | (state.zeta + epsilon_1);
        let increment = ((s1_high_bits | epsilon_1) - state * C) | W(1);

        for _ in 0..4 {
            state = state * C + increment;
        }

        predictions.record_prediction(compute_xsh_rr(state));

        candidates.push(CandidateState {
            epsilon_1,
            increment,
            state,
        })
    }

    println!("[+] Enumeration complete.");

    if !args.recovery {
        println!("[+] Output #5 will be one of the following:");

        predictions.display_predictions();
    }

    for (outputs, line) in stdin.lock().lines().enumerate() {
        let output = parse_output(&line?)?;

        if !args.recovery {
            println!(
                "[-] Reading output #{} (with value 0x{:08X})",
                outputs + 5,
                output
            );
        }

        if candidates.len() == 1 {
            println!(
                "[+] Generator internal state fully recovered after {:.2} seconds:",
                start_time.elapsed().as_secs_f64()
            );

            display_recovered_state(
                state.s0_high_bits | (state.zeta + candidates[0].epsilon_1),
                candidates[0].increment,
            );

            return Ok(());
        }

        if let Some(prediction) = predictions.unique_prediction() {
            if output != prediction {
                return Err(Error::new(
                    ErrorKind::Other,
                    "output sequence not produced by PCG-XSH-RR",
                ));
            }
        } else {
            candidates.retain(|candidate_state| compute_xsh_rr(candidate_state.state) == output);

            if candidates.is_empty() {
                return Err(Error::new(
                    ErrorKind::Other,
                    "output sequence not produced by PCG-XSH-RR",
                ));
            }

            if args.recovery {
                println!(
                    "[+] Pruned to {} states after {} outputs and {:.2} seconds.",
                    candidates.len(),
                    outputs + 5,
                    start_time.elapsed().as_secs_f64()
                );
            }
        }

        predictions.clear_predictions();

        for candidate_state in &mut candidates {
            candidate_state.state = candidate_state.state * C + candidate_state.increment;

            predictions.record_prediction(compute_xsh_rr(candidate_state.state));
        }

        if !args.recovery {
            println!("[+] Output #{} will be one of the following:", outputs + 6);

            predictions.display_predictions();
        }
    }

    println!("[-] No more outputs available, below is a candidate internal state:");

    display_recovered_state(candidates[0].state, candidates[0].increment);

    println!("[-] WARNING: this may eventually diverge from the target generator.");

    Ok(())
}

fn main() {
    if let Err(err) = run(Opt::from_args()) {
        eprintln!("\nfatal error: {}", err);
        process::exit(1); // report failure
    }
}

#[derive(StructOpt)]
#[structopt(about)]
struct Opt {
    #[structopt(long = "recovery")]
    recovery: bool,

    #[structopt(parse(from_os_str))]
    table: PathBuf,
}

const ASCII_HEADER: &str = r#"
  ____   ___  ___    ____  ____  ____   __   __ _  ____  ____
 (  _ \ / __)/ __)  (  _ \(  _ \(  __) / _\ (  / )(  __)(  _ \
  ) __/( (__( (_ \   ) _ ( )   / ) _) /    \ )  (  ) _)  )   /
 (__)   \___)\___/  (____/(__\_)(____)\_/\_/(__\_)(____)(__\_)

         PCG-XSH-RR Output Prediction & State Recovery
"#;
