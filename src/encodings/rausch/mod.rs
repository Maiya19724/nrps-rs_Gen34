// License: GNU Affero General Public License v3 or later
// A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

mod alpha_helix;
mod beta_sheet;
mod beta_turn;
mod hydrogenbond;
mod hydrophobicity_neu1;
mod hydrophobicity_neu2;
mod hydrophobicity_neu3;
mod isoelectric;
mod polar_grantham;
mod polar_radzicka;
mod polar_zimmerman;
mod volume;

pub fn encode(sequence: &str) -> Vec<f64> {
    let capacity = sequence.len() * 12;
    let encoded: Vec<f64> = Vec::with_capacity(capacity);
    sequence
        .chars()
        .map(encode_one)
        .fold(encoded, |mut acc, mut part| {
            acc.append(&mut part);
            acc
        })
}

// NRPSPredictor 2 uses {4,5,6,7,11,10,9,12,3,2,1,8} as the feature order
pub fn encode_one(c: char) -> Vec<f64> {
    vec![
        hydrogenbond::get(c),
        hydrophobicity_neu1::get(c),
        hydrophobicity_neu2::get(c),
        hydrophobicity_neu3::get(c),
        polar_zimmerman::get(c),
        polar_radzicka::get(c),
        polar_grantham::get(c),
        volume::get(c),
        beta_turn::get(c),
        beta_sheet::get(c),
        alpha_helix::get(c),
        isoelectric::get(c),
    ]
}

pub fn legacy_encode(sequence: &str) -> Vec<f64> {
    let capacity = sequence.len() * 12;
    let mut encoded: Vec<f64> = Vec::with_capacity(capacity);

    let mut array: Vec<Vec<f64>> = Vec::with_capacity(12);

    for c in sequence.chars() {
        array.push(encode_one(c));
    }

    for i in 0_usize..12 {
        for a in array.iter().take(sequence.len()) {
            encoded.push(a[i]);
        }
    }

    encoded
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use phf::phf_map;

    static DATA: phf::Map<char, [f64; 12]> = phf_map! {
        'A' => [-0.838548, 0.004378, -1.165652, 0.843010, -0.636447, -0.311135, -0.085800, -1.378285, -0.703066, -0.552985, 1.533010, -0.015368, ],
        'C' => [-0.838548, -0.900312, -1.873538, -1.271292, -0.567156, 1.304263, -1.077266, -1.046168, -0.088058, 0.449910, -1.095007, -0.566286, ],
        'D' => [0.147979, 1.332228, -0.363381, -0.078609, 1.690419, -1.152961, 1.782732, -0.696571, 1.477417, -1.360873, 0.036500, -1.888490, ],
        'E' => [0.147979, 1.157127, -0.457766, -0.783376, 1.699783, -1.118833, 1.515799, -0.074789, -0.116013, -1.834463, 1.861512, -1.627529, ],
        'F' => [-0.838548, -1.527758, 0.863621, 0.300882, -0.620061, 1.201878, -1.191666, 1.166279, -0.926706, 0.979216, 0.474503, -0.316923, ],
        'G' => [-0.838548, 0.383763, -1.495999, 1.005649, -0.636447, -0.709296, 0.257400, -2.005061, 1.589237, -0.775851, -1.569510, -0.032765, ],
        'H' => [0.147979, 0.135703, -0.127419, -1.054441, 1.779374, 0.314547, 0.791266, 0.369698, -0.116013, -0.441553, 0.000000, 0.906696, ],
        'I' => [-0.838548, -1.542349, -0.127419, 0.029817, -0.630361, 1.247382, -1.191666, 0.467086, -1.457849, 1.592097, 0.292002, -0.003769, ],
        'K' => [1.134506, 1.376003, 1.524314, 0.084030, 1.681055, -1.585251, 1.134466, 0.551988, 0.554905, -0.803709, 0.584004, 2.153511, ],
        'L' => [-0.838548, -1.294289, 0.250120, -0.458099, -0.630361, 1.133622, -1.306066, 0.469583, -1.373984, 0.756350, 0.766505, -0.026966, ],
        'M' => [-0.838548, -1.075413, -0.033035, -0.295460, -0.569497, 1.201878, -1.000999, 0.544497, -1.094435, 0.059895, 1.642511, -0.166145, ],
        'N' => [1.134506, 0.938250, -0.268997, 0.409307, -0.478201, -0.788928, 1.248866, -0.511784, 1.309688, -0.385836, -1.204508, -0.357517, ],
        'P' => [-0.838548, 0.573456, 1.099583, 1.710417, -0.562474, 1.008486, -0.123933, -0.556733, 1.589237, -1.333015, -1.569510, 0.158606, ],
        'Q' => [1.134506, 0.952842, -0.080227, 0.138243, -0.471179, -1.073329, 0.829399, 0.105004, -0.088058, 0.199186, 0.401503, -0.218338, ],
        'R' => [3.107561, 1.084168, 0.910813, -2.735040, 1.798101, -1.198465, 0.829399, 1.218719, 0.051717, -0.274403, -0.073000, 2.745023, ],
        'S' => [0.147979, 0.617232, -0.693728, 0.734585, -0.558261, -0.811680, 0.333666, -1.243440, 1.225823, -0.775851, -0.839506, -0.200940, ],
        'T' => [0.147979, 0.471314, -0.457766, 1.059862, -0.558729, -0.550031, 0.104867, -0.591692, -0.032148, 0.449910, -0.620504, -0.212539, ],
        'V' => [-0.838548, -1.177555, -0.882498, -0.349673, -0.630361, 0.997110, -0.924733, -0.154697, -1.122390, 1.870679, 0.219001, -0.038564, ],
        'W' => [0.147979, -0.914904, 1.477122, 1.330926, -0.538129, 0.758213, -1.115399, 2.072733, -1.094435, 0.951358, 0.292002, -0.079158, ],
        'X' => [0.850000, 0.057000, -0.003000, 0.094500, 13.594000, 0.213500, 8.325000, 145.195000, 0.991500, 1.028500, 1.000000, 6.026500, ],
        'Y' => [0.147979, -0.593885, 1.901853, -0.620738, -0.561070, 0.132531, -0.810333, 1.293632, 0.415131, 1.229940, -1.131507, -0.212539, ],
        '-' => [0.850000, 0.057000, -0.003000, 0.094500, 13.594000, 0.213500, 8.325000, 145.195000, 0.991500, 1.028500, 1.000000, 6.026500, ],
    };

    static LEGACY_CONCAT_DATA: phf::Map<&str, [f64; 24]> = phf_map! {
        "AC" => [-0.838548, -0.838548, 0.004378, -0.900312, -1.165652, -1.873538, 0.843010, -1.271292, -0.636447, -0.567156, -0.311135, 1.304263, -0.085800, -1.077266, -1.378285, -1.046168, -0.703066, -0.088058, -0.552985, 0.449910, 1.533010, -1.095007, -0.015368, -0.566286, ]
    };

    #[test]
    fn test_rausch_encoder() {
        for (c, expected) in DATA.entries() {
            let query = c.to_string();
            let got = encode(&query);
            dbg!(&got);
            for (i, value) in got.iter().enumerate() {
                assert_approx_eq!(value.clone(), expected[i]);
            }
        }
    }

    #[test]
    fn test_rausch_legacy_concatenation() {
        for (c, expected) in LEGACY_CONCAT_DATA.entries() {
            let query = c.to_string();
            let got = legacy_encode(&query);
            dbg!(&got);
            for (i, value) in got.iter().enumerate() {
                assert_approx_eq!(value.clone(), expected[i]);
            }
        }
    }
}
