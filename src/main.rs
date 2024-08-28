#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::{Parser, Subcommand};

use mlua::Lua;
use rust_htslib::bcf::Read;

use vcfexpress::{variant::HeaderMap, vcfexpress::VCFExpress};

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

        /// boolean Lua expression(s) to filter the VCF or BCF file
        #[arg(short, long)]
        expression: Vec<String>,

        /// expression(s) to set existing INFO field(s) (new ones can be added in prelude)
        /// e.g. --set-expression "AFmax=math.max(variant:info('AF'), variant:info('AFx'))"
        #[arg(short = 's', long)]
        set_expression: Vec<String>,

        /// template expression in luau: https://luau-lang.org/syntax#string-interpolation. e.g. '{variant.chrom}:{variant.pos}'
        #[arg(short, long)]
        template: Option<String>,

        /// File(s) containing lua(u) code to run once before any variants are processed.
        /// `header` is available here to access or modify the header.
        #[arg(short = 'p', long)]
        lua_prelude: Vec<String>,

        /// Optional output file. Default is stdout.
        #[arg(short, long)]
        output: Option<String>,

        /// Run lua code in https://luau.org/sandbox.
        #[arg(short = 'b', long)]
        sandbox: bool,
    },
}

fn filter_main(
    path: String,
    expressions: Vec<String>,
    set_expression: Vec<String>,
    template: Option<String>,
    lua_prelude: Vec<String>,
    output: Option<String>,
    sandbox: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let lua = Lua::new();

    let mut vcfexpr = VCFExpress::new(
        &lua,
        path,
        expressions,
        set_expression,
        template,
        lua_prelude,
        output,
        sandbox,
    )?;

    let mut reader = vcfexpr.reader();
    let mut writer = vcfexpr.writer();

    let header_map = HeaderMap::new();

    for record in reader.records() {
        let mut record = record?;
        writer.translate(&mut record);
        let mut sob = vcfexpr.evaluate(record, header_map.clone())?;
        writer.write(&mut sob)?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    match args.command {
        Some(Commands::Filter {
            path,
            expression,
            set_expression,
            template,
            lua_prelude,
            output,
            sandbox,
        }) => {
            filter_main(
                path,
                expression,
                set_expression,
                template,
                lua_prelude,
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
