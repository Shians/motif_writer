use std::collections::HashMap;
use polars::prelude::*;
use bio::io::fasta;
use std::fs::File;
use std::io::Write;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let pos_file = &args[1];
    let fasta_file = &args[2];
    let output_file = &args[3];

    let reader = fasta::Reader::from_file(fasta_file).unwrap();
    
    // store fasta sequences as hashmap
    let mut seqs: HashMap<String, String> = HashMap::new();
    for result in reader.records() {
        let record = result.expect("fasta record must be valued");
        let header = record.id().to_string();
        let seq = std::str::from_utf8(record.seq()).unwrap().to_string();
        
        seqs.insert(header, seq);
    }

    // column names for modkit column subset file
    let col_names = [
        "chr",
        "pos",
        "strand",
        "coverage",
        "methy_frac"
    ];

    // read in dataframe
    let mut df = CsvReader::from_path(pos_file)
        .unwrap()
        .with_separator(b'\t')
        .has_header(false)
        .finish()
        .unwrap();

    df.set_column_names(&col_names).expect("column name setting failed.");

    let mut outfile = File::create(output_file).unwrap();
    // iterate over rows in dataframe
    for i in 0..df.height() {
        let chr = df["chr"].str().unwrap().get(i).unwrap();
        let pos = df["pos"].i64().unwrap().get(i).unwrap() as usize;
        let strand = df["strand"].str().unwrap().get(i).unwrap();
        let coverage = df["coverage"].i64().unwrap().get(i).unwrap();
        let methy_frac = df["methy_frac"].f64().unwrap().get(i).unwrap();

        if (pos < 5) | (pos > seqs[chr].len() - 5) {
            continue;
        }

        // get 5bp context surrounding pos
        let chr_seq = &seqs[chr];
        let start = pos - 5;
        let end = pos + 4;
        let context = &chr_seq[start..end];

        writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}", chr, pos, strand, coverage, methy_frac, context).expect("write failed");
    }
}
