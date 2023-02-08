// Autogenerated by aaindex2rust.py, do not change manually
// POLAR_GRANTHAM amino acid featurisation
use phf::phf_map;

use crate::encodings::get_value;

static POLAR_GRANTHAM_MAP: phf::Map<char, f64> = phf_map! {
    'A' => 8.10,
    'L' => 4.90,
    'R' => 10.50,
    'K' => 11.30,
    'N' => 11.60,
    'M' => 5.70,
    'D' => 13.00,
    'F' => 5.20,
    'C' => 5.50,
    'P' => 8.00,
    'Q' => 10.50,
    'S' => 9.20,
    'E' => 12.30,
    'T' => 8.60,
    'G' => 9.00,
    'W' => 5.40,
    'H' => 10.40,
    'Y' => 6.20,
    'I' => 5.20,
    'V' => 5.90,

};
const POLAR_GRANTHAM_MEAN: f64 = 8.325;
const POLAR_GRANTHAM_STDEV: f64 = 2.62237964452136;

pub fn get(c: char) -> f64 {
    get_value(
        &POLAR_GRANTHAM_MAP,
        c,
        POLAR_GRANTHAM_MEAN,
        POLAR_GRANTHAM_STDEV,
        true,
    )
}
