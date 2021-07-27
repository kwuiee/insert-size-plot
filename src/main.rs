extern crate byteorder;
#[macro_use]
extern crate clap;
extern crate flate2;
extern crate plotlib;
extern crate plotters;
extern crate serde;
extern crate serde_json;

use std::fs::File;
use std::io::ErrorKind::{InvalidData, UnexpectedEof};
use std::io::{BufRead, BufReader, Error, ErrorKind, Read, Result};

use byteorder::{LittleEndian, ReadBytesExt};
use clap::{App, AppSettings};
use flate2::read::MultiGzDecoder;
use plotters::prelude::*;
use serde::ser::{Serialize, SerializeStruct, Serializer};

use plotlib::page::Page;
use plotlib::repr::Plot;
use plotlib::style::{LineJoin, LineStyle};
use plotlib::view::ContinuousView;

/// Read is paired, first in pair, properly mapped.
const P_FLAG: u16 = 0x1 + 0x2 + 0x40;
/// Read is secondary or supplementary.
const N_FLAG: u16 = 0x100 + 0x800;

fn opterr() -> std::io::Error {
    Error::new(InvalidData, "Option error.")
}

struct BamReader<T: BufRead> {
    reader: T,
}

impl BamReader<BufReader<MultiGzDecoder<File>>> {
    /// Read a bam file from path.
    fn from_path(v: &str) -> Result<Self> {
        let mut file = BufReader::with_capacity(16 * 1024, MultiGzDecoder::new(File::open(v)?));

        // Magic header.
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if magic != [b'B', b'A', b'M', 1] {
            return Err(Error::new(InvalidData, "Wrong BAM magic."));
        };

        // Header text.
        let l_text = file.read_i32::<LittleEndian>()?;
        let mut _text = vec![0u8; l_text as usize];
        file.read_exact(&mut _text)?;

        // Reference and length.
        let n_ref: u32 = file.read_u32::<LittleEndian>()?;
        for _ in 0..n_ref {
            let block_size: usize =
                file.read_u32::<LittleEndian>()? as usize + std::mem::size_of::<u32>();
            let mut _ref_entry = vec![0u8; block_size];
            file.read_exact(&mut _ref_entry)?;
        }

        Ok(Self { reader: file })
    }
}

impl<T: BufRead> BamReader<T> {
    /// Read a record (one line of bam).
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        let mut rem_size = match self.reader.read_u32::<LittleEndian>() {
            Ok(value) => value as usize,
            Err(e) => {
                if e.kind() == UnexpectedEof {
                    return Ok(false);
                } else {
                    return Err(e);
                }
            }
        };
        let mut _sink8 = [0u8; 1];
        let mut _sink16 = [0u8; 2];
        let mut _sink32 = [0u8; 4];
        // Ref id.
        record.set_ref_id(self.reader.read_i32::<LittleEndian>()?);
        // Ref position.
        self.reader.read_exact(&mut _sink32)?;
        // Query name length.
        let _l_name = self.reader.read_u8()? as usize;
        // Mapq.
        self.reader.read_exact(&mut _sink8)?;
        // Bin.
        self.reader.read_exact(&mut _sink16)?;
        // Number of operations in CIGAR.
        let _l_cigar = self.reader.read_u16::<LittleEndian>()? as usize;
        // Flag.
        record.set_flag(self.reader.read_u16::<LittleEndian>()?);
        // Sequence length.
        let _l_seq = self.reader.read_u32::<LittleEndian>()? as usize;
        // Mate ref id.
        record.set_mate_ref_id(self.reader.read_i32::<LittleEndian>()?);
        // Mate pos
        self.reader.read_exact(&mut _sink32)?;
        // Template length.
        record.set_tlen(self.reader.read_i32::<LittleEndian>()?);
        // Query name
        self.reader.read_exact(&mut vec![0u8; _l_name])?;
        // Cigar.
        self.reader
            .read_u32_into::<LittleEndian>(&mut vec![0u32; _l_cigar])?;
        // Sequence.
        self.reader.read_exact(&mut vec![0u8; (_l_seq + 1) / 2])?;
        // Quality.
        self.reader.read_exact(&mut vec![0u8; _l_seq])?;
        rem_size -= 32 + _l_name + _l_cigar * 4 + (_l_seq + 1) / 2 + _l_seq;
        // Optinal fields.
        self.reader.read_exact(&mut vec![0u8; rem_size])?;
        Ok(true)
    }
}

/// Compact read record.
#[derive(Default)]
struct Record {
    ref_id: i32,
    mate_ref_id: i32,
    tlen: i32,
    flag: u16,
}

impl Record {
    fn flag(&self) -> &u16 {
        &self.flag
    }

    fn set_flag(&mut self, v: u16) {
        self.flag = v
    }

    fn tlen(&self) -> &i32 {
        &self.tlen
    }

    fn set_tlen(&mut self, v: i32) {
        self.tlen = v
    }

    fn ref_id(&self) -> &i32 {
        &self.ref_id
    }

    fn set_ref_id(&mut self, v: i32) {
        self.ref_id = v
    }

    fn mate_ref_id(&self) -> &i32 {
        &self.mate_ref_id
    }

    fn set_mate_ref_id(&mut self, v: i32) {
        self.mate_ref_id = v
    }
}

#[derive(Default)]
struct Summary {
    // Pair count all.
    all_count: u32,
    // Insert size mean in all.
    all_mean: f64,
    // Pair count.
    count: u32,
    // Insert size mean.
    mean: f64,
    // Insert size standard deviation.
    std: f64,
    // First quantile.
    q1: usize,
    // Second quantile.
    q2: usize,
    // Third quantile.
    q3: usize,
}

impl Serialize for Summary {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Color", 8)?;
        state.serialize_field("Total count", &self.all_count)?;
        state.serialize_field(
            "Total mean insert size",
            &format!("{:.2}", self.all_mean).parse::<f64>().unwrap(),
        )?;
        state.serialize_field("Qualified read count", &self.count)?;
        state.serialize_field(
            "Qualified mean insize size",
            &format!("{:.2}", self.mean).parse::<f64>().unwrap(),
        )?;
        state.serialize_field(
            "Qualified insert size SD",
            &format!("{:.2}", self.std).parse::<f64>().unwrap(),
        )?;
        state.serialize_field("Qualified Q1", &self.q1)?;
        state.serialize_field("Qualified Q2", &self.q2)?;
        state.serialize_field("Qualified Q3", &self.q3)?;
        state.end()
    }
}

/// Get a proper upper limit value for figure axis.
fn round_max(mut v: f64) -> f64 {
    let mut digits = 0i32;
    if v >= 10f64 {
        while v >= 10f64 {
            v /= 10f64;
            digits += 1;
        }
    } else if v < 1f64 {
        while v < 1f64 {
            v *= 10f64;
            digits -= 1;
        }
    }
    (v.ceil() + 0.1f64) * 10f64.powi(digits)
}

fn cli(bam: &str, pic: &str, upper: &usize, format: &PicFormat) -> Result<()> {
    let mut data = vec![0u32; *upper + 1];
    let mut record = Record::default();
    let mut reader = BamReader::from_path(bam)?;
    let mut sum = Summary::default();

    while reader.read_into(&mut record)? {
        if record.flag() & P_FLAG != P_FLAG
            || record.flag() & N_FLAG != 0
            || record.ref_id() != record.mate_ref_id()
        {
            continue;
        };
        let tlen = record.tlen().abs() as usize;
        sum.all_mean += tlen as f64;
        sum.all_count += 1;
        if &tlen > upper {
            continue;
        };
        data[tlen] += 1;
        sum.mean += tlen as f64;
        sum.count += 1;
    }
    sum.all_mean /= sum.all_count as f64;
    sum.mean /= sum.count as f64;

    // Calculate quantiles and std.
    let mut quantiles = {
        let tmp = sum.count as f64;
        vec![
            ((tmp * 0.25f64) as u32, 0usize),
            ((tmp * 0.5f64) as u32, 0usize),
            ((tmp * 0.75f64) as u32, 0usize),
        ]
    };

    let mut ri = 0usize;
    let mut flag = true;
    let mut accum: u32 = 0;

    let mut height_max: u32 = 0;
    data.iter().enumerate().for_each(|(k, v)| {
        accum += v;
        if flag {
            let (index, value) = &mut quantiles[ri];
            if accum > *index {
                *value = k;
                ri += 1;
            };
            flag = ri < quantiles.len();
        };
        sum.std += (k as f64 - sum.mean).powi(2);
        height_max = u32::max(*v, height_max);
    });
    sum.std = (sum.std / (sum.count as f64)).powf(0.5f64);
    unsafe {
        sum.q1 = quantiles.get_unchecked(0).1;
        sum.q2 = quantiles.get_unchecked(1).1;
        sum.q3 = quantiles.get_unchecked(2).1;
    }
    let height_max: f64 = round_max((height_max as f64) / (sum.count as f64));

    // Plot line.
    match format {
        PicFormat::Svg => {
            let line = Plot::new(
                data.into_iter()
                    .enumerate()
                    .map(|(i, j)| (i as f64, (j as f64) / (sum.count as f64)))
                    .collect(),
            )
            .line_style(
                LineStyle::new()
                    .colour("#FF0000")
                    .linejoin(LineJoin::Round)
                    .width(1.0),
            );
            let view = ContinuousView::new()
                .add(line)
                .x_label("插入片段大小(bp)")
                .y_label("比例");
            Page::single(&view)
                .save(pic)
                .map_err(|_| Error::new(InvalidData, format!("Failed to write {}", pic)))?;
        }
        PicFormat::Png => {
            let root = BitMapBackend::new(pic, (700, 610)).into_drawing_area();
            root.fill(&WHITE)
                .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;

            let mut chart = ChartBuilder::on(&root)
                .x_label_area_size(35)
                .y_label_area_size(40)
                .margin(5)
                .build_cartesian_2d(
                    (0f64..((upper + 1) as f64))
                        .step(1.0)
                        .use_round()
                        .into_segmented(),
                    0f64..height_max,
                )
                .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;

            chart
                .configure_mesh()
                .disable_mesh()
                .bold_line_style(&WHITE.mix(0.3))
                .x_desc("插入片段大小(bp)")
                .y_desc("比例")
                .axis_desc_style((FontFamily::Name("WenQuanYi Zen Hei"), 20))
                .draw()
                .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;

            chart
                .draw_series(LineSeries::new(
                    data.into_iter().enumerate().map(|(i, j)| {
                        (
                            SegmentValue::Exact(i as f64),
                            (j as f64) / (sum.count as f64),
                        )
                    }),
                    RED.stroke_width(2),
                ))
                .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
        }
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&sum).map_err(|e| Error::new(InvalidData, e))?
    );
    Ok(())
}

enum PicFormat {
    Svg,
    Png,
}

impl PicFormat {
    fn from_str(v: &str) -> Result<Self> {
        if v.ends_with(".svg") || v.ends_with(".SVG") {
            Ok(Self::Svg)
        } else if v.ends_with(".png") || v.ends_with(".PNG") {
            Ok(Self::Png)
        } else {
            Err(Error::new(ErrorKind::InvalidData, "No such option."))
        }
    }
}

fn main() -> Result<()> {
    let opts = App::new(crate_name!())
        .author(crate_authors!())
        .about(crate_description!())
        .version(crate_version!())
        .setting(AppSettings::ArgRequiredElseHelp)
        .args_from_usage(
            "
            <pic> -o=[FILE] 'Output pic file path, support `.svg` and `.png` suffix.'
            [upper] -m=[NUMBER] 'Maximum insert size to record, default 500. !Bigger number costs more memory!.'
            <bam> 'Input bam file.'
            ",
        )
        .get_matches();
    let bam: &str = opts.value_of("bam").ok_or_else(opterr)?;
    let pic: &str = opts.value_of("pic").ok_or_else(opterr)?;
    let format: PicFormat = PicFormat::from_str(pic)?;
    let upper: usize = opts
        .value_of("upper")
        .unwrap_or("500")
        .parse()
        .map_err(|_| opterr())?;
    cli(bam, pic, &upper, &format)
}
