use clap::{Parser, Subcommand};

use mlua::Lua;
use rust_htslib::bcf::Read;

use vcfexpr::vcfexpr::VCFExpr;

/// Args take the arguments for clap.
/// Accept the path to VCF or BCF and the lua expressions
#[derive(Parser)]
#[command(version, about, author)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Filter a VCF/BCF and optionally print by template expression.
    /// If no template is given the output will be VCF/BCF
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

        /// File(s) containing lua code to load. May contain functions that will be called by the expressions.
        #[arg(short, long)]
        lua: Vec<String>,

        /// File containing lua code to run once before any variants are processed.
        #[arg(short = 'p', long)]
        lua_prelude: Option<String>,

        /// Optional output file. Default is stdout.
        #[arg(short, long)]
        output: Option<String>,
    },
}

fn filter_main(
    path: String,
    expressions: Vec<String>,
    set_expression: Vec<String>,
    template: Option<String>,
    lua_code: Vec<String>,
    lua_prelude: Option<String>,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let lua = Lua::new();

    let mut vcfexpr = VCFExpr::new(
        &lua,
        path,
        expressions,
        set_expression,
        template,
        lua_prelude,
        output,
    )?;
    for path in lua_code {
        vcfexpr.add_lua_code(&path)?;
    }
    let mut reader = vcfexpr.reader();
    let mut writer = vcfexpr.writer();

    for record in reader.records() {
        let mut record = record?;
        writer.translate(&mut record);
        let mut sob = vcfexpr.evaluate(record)?;
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
            lua: lua_code,
            lua_prelude,
            output,
        }) => {
            filter_main(
                path,
                expression,
                set_expression,
                template,
                lua_code,
                lua_prelude,
                output,
            )?;
        }
        None => {
            println!("No command provided");
        }
    }
    Ok(())
}
