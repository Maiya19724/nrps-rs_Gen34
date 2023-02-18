// License: GNU Affero General Public License v3 or later
// A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
pub mod predictions;
pub mod stachelhaus;

use std::fs::File;
use std::path::PathBuf;

use walkdir::WalkDir;

use crate::config::Config;
use crate::errors::NrpsError;
use crate::svm::models::SVMlightModel;
use predictions::{ADomain, Prediction, PredictionCategory};

#[derive(Debug)]
pub struct Predictor {
    pub models: Vec<SVMlightModel>,
}

impl Predictor {
    pub fn predict(&self, domains: &mut Vec<ADomain>) -> Result<(), NrpsError> {
        for model in self.models.iter() {
            for domain in domains.iter_mut() {
                let score = model.predict_seq(&domain.aa34)?;
                if score > 0.0 {
                    let pred = Prediction {
                        name: model.name.to_string(),
                        score,
                    };
                    domain.add(model.category, pred);
                }
            }
        }
        Ok(())
    }
}

pub fn load_models(config: &Config) -> Result<Vec<SVMlightModel>, NrpsError> {
    let mut models = Vec::with_capacity(1000);

    for category_dir_res in WalkDir::new(&config.model_dir())
        .min_depth(1)
        .max_depth(1)
        .sort_by_file_name()
    {
        let category_dir = category_dir_res?;
        let category = match category_dir.file_name().to_str().unwrap() {
            "NRPS3_THREE_CLUSTER" => PredictionCategory::ThreeCluster,
            "NRPS3_LARGE_CLUSTER" => PredictionCategory::LargeCluster,
            "NRPS3_SMALL_CLUSTER" => PredictionCategory::SmallCluster,
            "NRPS3_SINGLE_CLUSTER" => PredictionCategory::Single,
            "NRPS2_THREE_CLUSTER" => PredictionCategory::LegacyThreeCluster,
            "NRPS2_THREE_CLUSTER_FUNGAL" => PredictionCategory::LegacyThreeClusterFungal,
            "NRPS2_LARGE_CLUSTER" => PredictionCategory::LegacyLargeCluster,
            "NRPS2_SMALL_CLUSTER" => PredictionCategory::LegacySmallCluster,
            "NRPS2_SINGLE_CLUSTER" => PredictionCategory::LegacySingle,
            _ => continue,
        };

        if !config.categories().contains(&category) {
            continue;
        }

        for model_file_res in WalkDir::new(category_dir.path())
            .min_depth(1)
            .max_depth(1)
            .sort_by_file_name()
        {
            let model_file = model_file_res?.path().to_path_buf();
            if let Some(ext) = model_file.extension() {
                if ext != "mdl" {
                    continue;
                }
            } else {
                continue;
            }
            let name = extract_name(&model_file);
            let handle = File::open(&model_file)?;
            models.push(SVMlightModel::from_handle(handle, name, category)?);
        }
    }

    Ok(models)
}

fn extract_name(filename: &PathBuf) -> String {
    let square_brackets: &[_] = &['[', ']'];
    filename
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .trim_matches(square_brackets)
        .to_string()
}
