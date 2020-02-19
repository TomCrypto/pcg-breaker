use bytemuck::cast_slice;
use std::fs::write;
use std::num::Wrapping as W;

/// The multiplicative constant from PCG-XSH-RR
const C: W<u64> = W(6_364_136_223_846_793_005);

fn main() {
    let mut table: Vec<u64> = Vec::with_capacity(0x800_0000);

    for zeta in 0..0x800_0000 {
        table.push(((((W(zeta) * C) >> 27) << 27) | W(zeta)).0);
    }

    table.sort_unstable();

    write("table.bin", cast_slice(&table)).unwrap();
}
