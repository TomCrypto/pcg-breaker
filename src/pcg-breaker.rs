use bytemuck::{bytes_of_mut, cast_slice_mut};
use itertools::iproduct;
use parse_int::parse;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind, Read, Result, Stdin};
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
const A: W<u64> = W(6_364_136_223_846_793_005);

/// The multiplicative inverse of A (modulo 2^64)
const A_INV: W<u64> = W(13_877_824_140_714_322_085);

fn parse_output(line: &str) -> Result<u32> {
    parse::<u32>(line).map_err(|err| Error::new(ErrorKind::Other, err))
}

fn read_output_lines(stdin: &mut BufReader<Stdin>) -> Result<u32> {
    let mut line = String::new();
    stdin.read_line(&mut line)?;
    parse_output(&line)
}

fn read_output_bytes(stdin: &mut BufReader<Stdin>) -> Result<u32> {
    let mut output = 0u32;
    stdin.read_exact(bytes_of_mut(&mut output))?;
    Ok(output)
}

struct LookupTable {
    table: Box<[u64]>,
}

impl LookupTable {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut raw_table_file = File::open(path)?;

        let mut table = vec![0; 0x800_0000].into_boxed_slice();
        raw_table_file.read_exact(cast_slice_mut(&mut table))?;

        Ok(Self { table })
    }

    /// Queries the lookup table for its value given N.
    pub fn query(&self, n: W<u64>) -> Option<W<u64>> {
        if let Some(beta) = self.scan_table(n) {
            return Some(beta);
        }

        if let Some(beta) = self.scan_table(-n) {
            return Some(-beta);
        }

        if n == W(1) {
            return None;
        }

        if let Some(beta) = self.scan_table(W(1) + n) {
            return Some(beta - W(0x800_0000));
        }

        if n == -W(1) {
            return None;
        }

        if let Some(beta) = self.scan_table(W(1) - n) {
            return Some(W(0x800_0000) - beta);
        }

        None
    }

    fn scan_table(&self, n: W<u64>) -> Option<W<u64>> {
        let n2 = n.0 & 0x1f_ffff_ffff;

        let estimate = (n2 >> 10) as isize;

        let lo = (estimate - 2).max(0) as usize;
        let hi = (estimate + 10).min(0x7ff_ffff) as usize;

        for entry in &self.table[lo..=hi] {
            if entry >> 27 == n2 {
                return Some(W(entry & 0x7ff_ffff));
            }
        }

        None
    }
}

struct Predictor {
    table: LookupTable,
    last_output: u32,
    triple: Triple,
}

impl Predictor {
    /// Initializes the output predictor with four initial outputs.
    fn new(table: LookupTable, outputs: [u32; 4]) -> Result<Self> {
        for (s0_rot, s1_rot, s2_rot) in iproduct!(0..32, 0..32, 0..32) {
            let s0_star = invert_xsh_rr(s0_rot, outputs[0]) >> 27;
            let s1_star = invert_xsh_rr(s1_rot, outputs[1]) >> 27;
            let s2_star = invert_xsh_rr(s2_rot, outputs[2]) >> 27;

            let n = (A * (s1_star - s0_star) + (s1_star - s2_star)) & W(0x1f_ffff_ffff);

            if let Some(beta) = table.query(n) {
                let epsilon_min: i32 = (beta.0 as i32).max(0);
                let epsilon_max: i32 = (beta.0 as i32 + 134_217_728).min(134_217_728);

                if Self::test_state(s1_star, s2_star, epsilon_min, epsilon_max, beta, outputs[3]) {
                    return Ok(Self {
                        table,
                        last_output: outputs[3],
                        triple: Triple {
                            sj_star: s1_star,
                            sk_star: s2_star,

                            beta,

                            epsilon_min,
                            epsilon_max,
                        },
                    });
                }
            }
        }

        Err(Error::new(
            ErrorKind::Other,
            "output sequence not produced by PCG-XSH-RR",
        ))
    }

    /// Submits the next output produced by the target PCG generator.
    pub fn submit_next_output(&mut self, output: u32) -> Result<()> {
        for sk_rot in 0..32 {
            let si_star = self.triple.sj_star;
            let sj_star = self.triple.sk_star;
            let sk_star = invert_xsh_rr(sk_rot, self.last_output) >> 27;

            let n = (A * (sj_star - si_star) + (sj_star - sk_star)) & W(0x1f_ffff_ffff);

            if let Some(beta) = self.table.query(n) {
                let epsilon_min = (self.triple.epsilon_min + beta.0 as i32).max(0);
                let epsilon_max = (self.triple.epsilon_max + beta.0 as i32).min(134_217_728);

                if Self::test_state(sj_star, sk_star, epsilon_min, epsilon_max, beta, output) {
                    self.triple = Triple {
                        sj_star,
                        sk_star,

                        beta,

                        epsilon_min,
                        epsilon_max,
                    };

                    self.last_output = output;
                    return Ok(());
                }
            }
        }

        Err(Error::new(
            ErrorKind::Other,
            "output sequence not produced by PCG-XSH-RR",
        ))
    }

    /// Returns either one or two future PCG outputs.
    pub fn predict_future_output(&self) -> [u32; 2] {
        let epsilon2 = self.triple.epsilon_min;
        let epsilon1 = W(epsilon2 as u64) - self.triple.beta;

        let sj = (self.triple.sj_star << 27) + epsilon1;
        let sk = (self.triple.sk_star << 27) + W(epsilon2 as u64);

        let increment = (sk - A * sj) | W(1);

        let state = A * sk + increment;

        let output1 = compute_xsh_rr(A * state + increment);

        let epsilon2 = self.triple.epsilon_max - 1;
        let epsilon1 = W(epsilon2 as u64) - self.triple.beta;

        let sj = (self.triple.sj_star << 27) + epsilon1;
        let sk = (self.triple.sk_star << 27) + W(epsilon2 as u64);

        let increment = (sk - A * sj) | W(1);

        let state = A * sk + increment;

        let output2 = compute_xsh_rr(A * state + increment);

        [output1, output2]
    }

    /// Returns the number of candidate states left.
    pub fn remaining_candidate_count(&self) -> usize {
        (self.triple.epsilon_max - self.triple.epsilon_min) as usize
    }

    /// Returns the set of all remaining candidate states.
    pub fn remaining_candidates(&self) -> Vec<FullState> {
        let mut states = Vec::with_capacity(self.remaining_candidate_count());

        for epsilon_k in self.triple.epsilon_min..self.triple.epsilon_max {
            let sj = (self.triple.sj_star << 27) + W(epsilon_k as u64) - self.triple.beta;
            let sk = (self.triple.sk_star << 27) + W(epsilon_k as u64);

            let inc = (sk - A * sj) | W(1);

            states.push(FullState {
                state: A * sk + inc,
                inc,
            });
        }

        states
    }

    fn test_state(
        sj_star: W<u64>,
        sk_star: W<u64>,
        min: i32,
        max: i32,
        beta: W<u64>,
        output: u32,
    ) -> bool {
        let epsilon1 = W(min as u64) - beta;
        let sj = (sj_star << 27) + epsilon1;
        let sk = (sk_star << 27) + W(min as u64);

        let increment = (sk - A * sj) | W(1);
        let min_state = A * sk + increment;

        if output == compute_xsh_rr(min_state) {
            return true;
        }

        let epsilon1 = W((max - 1) as u64) - beta;
        let sj = (sj_star << 27) + epsilon1;
        let sk = (sk_star << 27) + W((max - 1) as u64);

        let increment = (sk - A * sj) | W(1);
        let max_state = A * sk + increment;

        if output == compute_xsh_rr(max_state) {
            return true;
        }

        false
    }
}

#[derive(Debug)]
struct FullState {
    state: W<u64>,
    inc: W<u64>,
}

struct Triple {
    sj_star: W<u64>,
    sk_star: W<u64>,

    beta: W<u64>,

    epsilon_min: i32,
    epsilon_max: i32,
}

fn display_predictions(count: usize, outputs: [u32; 2]) {
    if outputs[0] == outputs[1] {
        println!("\n[+] Output #{} will be 0x{:08X}\n", count, outputs[0]);
    } else {
        println!(
            "\n[+] Output #{} will be 0x{:08X} OR 0x{:08X}\n",
            count, outputs[0], outputs[1]
        );
    }
}

fn run(args: Opt) -> Result<()> {
    println!("{}", ASCII_HEADER);

    println!("[-] Starting clock.");
    let start_time = Instant::now();

    let table = LookupTable::open(&args.table).map_err(|err| {
        println!("[!] Failed to load precomputed table!");
        err
    })?;

    println!("[+] Loaded precomputed table.");

    println!("[-] Reading 4 outputs to initialize the predictor.");

    let mut stdin = BufReader::new(std::io::stdin());

    let read_output = if args.binary {
        read_output_bytes
    } else {
        read_output_lines
    };

    let mut predictor = Predictor::new(
        table,
        [
            read_output(&mut stdin)?,
            read_output(&mut stdin)?,
            read_output(&mut stdin)?,
            read_output(&mut stdin)?,
        ],
    )?;

    println!(
        "[+] Predictor initialized after {:.2} seconds.",
        start_time.elapsed().as_secs_f64(),
    );

    if !args.recovery {
        display_predictions(5, predictor.predict_future_output());
    }

    let mut remaining_candidates = vec![];
    let mut outputs = 4;

    while let Ok(output) = read_output(&mut stdin) {
        outputs += 1;

        if !args.recovery {
            println!(
                "[-] Reading output #{} (with value 0x{:08X})",
                outputs,
                output
            );
        }

        if remaining_candidates.is_empty() {
            const THRESHOLD: usize = 1000;

            let count = predictor.remaining_candidate_count();

            predictor.submit_next_output(output)?;

            if args.recovery && predictor.remaining_candidate_count() != count {
                println!(
                    "[+] Pruned to {} states after {} outputs and {:.2} seconds.",
                    predictor.remaining_candidate_count(),
                    outputs,
                    start_time.elapsed().as_secs_f64()
                );
            }

            if args.recovery && predictor.remaining_candidate_count() <= THRESHOLD {
                remaining_candidates = predictor.remaining_candidates();
            } else if !args.recovery {
                display_predictions(outputs + 1, predictor.predict_future_output());
            }
        } else {
            for state in &mut remaining_candidates {
                state.state = A * state.state + state.inc;
            }

            let count = remaining_candidates.len();

            remaining_candidates.retain(|state| compute_xsh_rr(state.state) == output);

            if remaining_candidates.len() != count {
                println!(
                    "[+] Pruned to {} states after {} outputs and {:.2} seconds.",
                    remaining_candidates.len(),
                    outputs,
                    start_time.elapsed().as_secs_f64()
                );
            }

            if remaining_candidates.is_empty() {
                return Err(Error::new(
                    ErrorKind::Other,
                    "output sequence not produced by PCG-XSH-RR",
                ));
            } else if let [recovered] = remaining_candidates.as_mut_slice() {
                println!("[-] State recovery complete, rewinding state...");

                for _ in 0..outputs - 1 {
                    recovered.state = A_INV * (recovered.state - recovered.inc);
                }

                println!(
                    "[+] Generator internal state fully recovered after {:.2} seconds:",
                    start_time.elapsed().as_secs_f64()
                );

                println!("\n    pcg32_random_t state = {{");
                println!("        .state = 0x{:016X}", recovered.state);
                println!("        .inc   = 0x{:016X}", recovered.inc);
                println!("    }};\n");

                return Ok(());
            }
        }
    }

    if args.recovery {
        println!("[-] Not enough outputs available to complete state recovery.");
    }

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

    #[structopt(long = "binary")]
    binary: bool,

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
