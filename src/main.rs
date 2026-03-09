mod args;
mod engine;
mod errors;

use args::Args;
use clap::Parser;
use engine::run;

use crate::errors::AppError;

fn main() -> Result<(), AppError> {
    let args = Args::parse();

    run(args)?;
    Ok(())
}
