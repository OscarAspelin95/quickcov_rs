use bio_utils_rs::errors::BioError;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("Invalid file extension: {0}")]
    InvalidFileExtension(String),

    #[error("Failed to read file: {0}")]
    FastaReadError(String),

    #[error("File does not exist: {0}")]
    FileNotFoundError(String),

    #[error("IO error: {0}")]
    IOError(#[from] std::io::Error),

    #[error(transparent)]
    BioError(#[from] BioError),

    #[error("Progress bar template error: {0}")]
    ProgressBarError(#[from] indicatif::style::TemplateError),
}
