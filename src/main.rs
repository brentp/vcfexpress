#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::{Parser, Subcommand};
use std::str::FromStr;

use rust_htslib::bcf::Read;

use vcfexpress::{variant::HeaderMap, vcfexpress::VCFExpress, script_engine::{ScriptLanguage, ScriptConfig, create_engine}};

/// Args take the arguments for clap.
/// Accept the path to VCF or BCF and the lua expressions
#[derive(Parser)]
#[command(version, about, author)]
#[command(arg_required_else_help(true))]
#[command(propagate_version = true)]
#[command(help_template = "
{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Filter a VCF/BCF and optionally print by template expression.
    /// If no template is given the output will be VCF/BCF
    #[command(arg_required_else_help(true))]
    #[command(help_template = "
{name} {version}
{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
")]
    Filter {
        /// Path to input VCF or BCF file
        path: String,

        /// Scripting language to use (lua or javascript). Default: lua
        #[arg(short = 'L', long, default_value = "lua")]
        language: String,

        /// boolean expression(s) to filter the VCF or BCF file. Syntax depends on selected language.
        #[arg(short, long)]
        expression: Vec<String>,

        /// expression(s) to set existing INFO field(s) (new ones can be added in prelude)
        /// e.g. --set-expression "AFmax=math.max(variant.info('AF'), variant.info('AFx'))"
        #[arg(short = 's', long)]
        set_expression: Vec<String>,

        /// template expression for output. Syntax depends on selected language.
        /// Lua: use Luau interpolation '{variant.chrom}'
        /// JavaScript: use template literals '${variant.chrom}'
        #[arg(short, long)]
        template: Option<String>,

        /// File(s) containing code to run once before any variants are processed.
        /// `header` is available here to access or modify the header.
        /// Use --lua-prelude for backward compatibility.
        #[arg(short = 'p', long, alias = "lua-prelude")]
        prelude: Vec<String>,

        /// Optional output file. Default is stdout.
        #[arg(short, long)]
        output: Option<String>,

        /// Run scripting code in sandbox mode.
        #[arg(short = 'b', long)]
        sandbox: bool,
    },
}

fn filter_main(
    path: String,
    language: String,
    expressions: Vec<String>,
    set_expression: Vec<String>,
    template: Option<String>,
    prelude: Vec<String>,
    output: Option<String>,
    sandbox: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // Parse the language
    let script_language = ScriptLanguage::from_str(&language)
        .map_err(|e| format!("Invalid language '{}': {}", language, e))?;

    // Create script configuration
    let config = ScriptConfig {
        sandbox,
        language: script_language,
    };

    // Create the appropriate scripting engine
    let engine = create_engine(&config)?;

    // Use the new engine-based constructor
    let mut vcfexpr = VCFExpress::new_with_engine(
        engine,
        path,
        expressions,
        set_expression,
        template,
        prelude,
        output,
    )?;

    let mut reader = vcfexpr.reader();
    let mut writer = vcfexpr.writer();
    let start_time = std::time::Instant::now();

    let header_map = HeaderMap::new();
    let header = reader.header().clone();
    let mut written = 0;
    let mut total = 0;
    for record in reader.records() {
        let mut record = record?;
        writer.translate(&mut record);
        let mut sob = vcfexpr.evaluate(record, &header, header_map.clone())?;
        written += writer.write(&mut sob)?;
        total += 1;
    }
    log::info!(
        "{} of {} records written ({:.2}%) in {:.1}s. evaluated {:.0} variants/second",
        written,
        total,
        written as f64 / total as f64 * 100.0,
        start_time.elapsed().as_millis() as f64 / 1000.0,
        1000.0 * (total as f64 / start_time.elapsed().as_millis() as f64)
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    match args.command {
        Some(Commands::Filter {
            path,
            language,
            expression,
            set_expression,
            template,
            prelude,
            output,
            sandbox,
        }) => {
            filter_main(
                path,
                language,
                expression,
                set_expression,
                template,
                prelude,
                output,
                sandbox,
            )?;
        }
        None => {
            println!("No command provided");
        }
    }
    Ok(())
}
