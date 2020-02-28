# PCG-Breaker

This program is able to observe the outputs of a [PCG-XSH-RR][1] pseudorandom number generator with a completely unknown state (both `state` and `inc` are unknown) and quickly predict future outputs with near 100% accuracy. Optionally, given enough outputs, it can perform state recovery of the generator under observation. Note that state recovery will require several dozen million consecutive outputs; this is not a flaw and _any_ state recovery algorithm will need that many outputs to completely pin down the PCG generator's internal state, as there is significant overlap in output sequences for different internal states. This is a well-known characteristic of LCGs with unknown increment.

The program uses a heavily specialized deterministic algorithm for efficiently tracking all possible internal states of the generator and discarding impossible states as additional outputs are observed. The algorithm will, almost all of the time, predict the next output given a sequence of outputs with perfect accuracy; in the rare cases it doesn't, it will predict exactly two possible outputs with typically very low Hamming distance to each other. The algorithm has the following performance characteristics:

 - Initializing the algorithm with four consecutive outputs requires up to 2^17 lookups into a precomputed table
 - Submitting additional outputs to the algorithm requires up to 128 table lookups into a precomputed table
 - Predicting the next output(s) at any given point requires only two PCG-XSH-RR operations

This demonstrates that the PCG-XSH-RR variant is not at all challenging to predict. Note that the algorithm as-is doesn't currently extend to any other variant of the PCG family of generators; the XSH-RR variant is particularly vulnerable due to the large amount of internal state it leaks.

## Usage

First precompute the lookup table using the `gen-table` binary:

```
$ cargo run --release --bin gen-table
```

This will write a 1GiB file of precomputed data into the current working directory called `table.bin`, which will later be used by the `pcg-breaker` binary.

Then execute the `pcg-breaker` binary with the path to the table, piping/typing in outputs of a PCG generator into its standard input. By default it will accept one ASCII number on each line, either in decimal or 0x-prefixed hexadecimal, but raw native-endian outputs can be accepted with the `--binary` flag. After being given four outputs, the program will begin predicting the next output in the sequence produced by the generator. A typical execution looks like this:

```text
[-] Starting clock.
[+] Loaded precomputed table.
[-] Reading 4 outputs to initialize the predictor.
[+] Predictor initialized after 0.50 seconds.

[+] Output #5 will be 0x5FAAB311 OR 0x5FAABD11

[-] Reading output #5 (with value 0x5FAABD11)

[+] Output #6 will be 0x3D7B6D05 OR 0x3D1B6D05

[-] Reading output #6 (with value 0x3D7B6D05)

[+] Output #7 will be 0xDFE18B58

[-] Reading output #7 (with value 0xDFE18B58)

[+] Output #8 will be 0x964867B9

[-] Reading output #8 (with value 0x964867B9)

[+] Output #9 will be 0xB1DE26E9
```

If you are only interested in recovering the internal state of a generator rather than predicting future outputs, you can pass the `--recovery` flag to the program. Note that recovering the state requires an unrealistically large number of outputs, but it can be used for testing by e.g. piping the output of the provided `pcg.c` program. After having recovered the generator's internal state, the program will display, for example:

```text
[+] Generator internal state fully recovered after 29.57 seconds:

    pcg32_random_t state = {
        .state = 0xBD094A5E7A8A7587
        .inc   = 0x24E8930796B7B111
    };
```

Note that the `state` displayed will be the state used to produce the _very first output_ given to PCG-breaker; you can advance it yourself if needed. Also note that since the lowest bit of the increment in the PCG state is always masked to 1 (and is therefore irrelevant) the program will conventionally report the recovered increment `inc` with its lowest bit masked to 1 as well.

## Performance

The algorithm is currently single-threaded, but is able to begin predicting outputs after only four consecutive outputs and less than one second of computation on commodity hardware. Full state recovery will occur in a minute or two if enough outputs are available. It is possible to parallelize virtually every part of this program though it runs adequately fast already.

Currently there is no feature to skip unknown outputs from a generator, so all outputs must be consecutive.

[1]: https://www.pcg-random.org/download.html
