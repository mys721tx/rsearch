//! OTU tables, as vsearch's: the abundance of each sample in each OTU, the
//! samples taken from query labels and the OTUs (and their taxonomy) from
//! target labels, written in vsearch's classic, mothur or BIOM 1.0 format.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use crate::label::find_attribute_with;

/// the value of the first of `names` in `header`, as vsearch's
/// `header_find_first_attribute`: an annotation starting the header or
/// following a `;`, with a value up to the next `;`.
fn first_value<'a>(header: &'a str, names: &[&str]) -> Option<&'a str> {
    names
        .iter()
        .filter_map(|name| {
            find_attribute_with(header.as_bytes(), name.as_bytes(), |c| c != b';', true)
                .map(|(start, end)| (start, &header[start + name.len()..end]))
        })
        .min_by_key(|&(start, _)| start)
        .map(|(_, value)| value)
}

/// the sample of a query label: its `sample=` or `barcodelabel=`, whichever
/// comes first, or its leading letters, digits and underscores.
///
/// # Example
/// ```rust
/// # use rsearch::otutable::sample;
/// assert_eq!(sample("s1.r1;size=5"), "s1");
/// assert_eq!(sample("r2;barcodelabel=bb;sample=cc"), "bb");
/// ```
pub fn sample(label: &str) -> &str {
    first_value(label, &["sample=", "barcodelabel="]).unwrap_or_else(|| {
        let end = label
            .bytes()
            .position(|c| !c.is_ascii_alphanumeric() && c != b'_')
            .unwrap_or(label.len());
        &label[..end]
    })
}

/// the OTU of a target label: its `otu=`, or the label up to its first `;`.
///
/// # Example
/// ```rust
/// # use rsearch::otutable::otu;
/// assert_eq!(otu("x;otu=C;tax=d:B"), "C");
/// assert_eq!(otu("t4;x"), "t4");
/// ```
pub fn otu(label: &str) -> &str {
    first_value(label, &["otu="]).unwrap_or_else(|| label.split(';').next().unwrap())
}

/// An OTU table being filled.
#[derive(Default)]
pub struct OtuTable {
    samples: BTreeSet<String>,
    otus: BTreeSet<String>,
    taxonomy: BTreeMap<String, String>,
    counts: BTreeMap<(String, String), u64>,
}

impl OtuTable {
    /// counts `abundance` of the sample of query label `query` in the OTU
    /// of target label `target`. Either may be missing. Only its sample or
    /// OTU, and the OTU's taxonomy, is then recorded.
    pub fn add(&mut self, query: Option<&str>, target: Option<&str>, abundance: u64) {
        let sample = query.map(sample);
        let otu = target.map(otu);
        if let (Some(target), Some(otu)) = (target, otu)
            && let Some(tax) = first_value(target, &["tax="])
        {
            self.taxonomy.insert(otu.to_string(), tax.to_string());
        }
        if let Some(sample) = sample {
            self.samples.insert(sample.to_string());
        }
        if let Some(otu) = otu {
            self.otus.insert(otu.to_string());
        }
        if let (Some(sample), Some(otu)) = (sample, otu)
            && abundance != 0
        {
            *self
                .counts
                .entry((otu.to_string(), sample.to_string()))
                .or_default() += abundance;
        }
    }

    fn count(&self, otu: &str, sample: &str) -> u64 {
        self.counts
            .get(&(otu.to_string(), sample.to_string()))
            .copied()
            .unwrap_or(0)
    }

    /// vsearch's `--otutabout`: an OTU a line, a sample a column, and the
    /// taxonomy last if any OTU has one.
    pub fn classic(&self) -> String {
        let mut out = String::from("#OTU ID");
        for sample in &self.samples {
            let _ = write!(out, "\t{sample}");
        }
        if !self.taxonomy.is_empty() {
            out.push_str("\ttaxonomy");
        }
        out.push('\n');
        for otu in &self.otus {
            out.push_str(otu);
            for sample in &self.samples {
                let _ = write!(out, "\t{}", self.count(otu, sample));
            }
            if !self.taxonomy.is_empty() {
                out.push('\t');
                out.push_str(self.taxonomy.get(otu).map_or("", String::as_str));
            }
            out.push('\n');
        }
        out
    }

    /// vsearch's `--mothur_shared_out`: a sample a line, labelled
    /// `vsearch`, an OTU a column.
    pub fn mothur(&self) -> String {
        let mut out = String::from("label\tGroup\tnumOtus");
        for otu in &self.otus {
            let _ = write!(out, "\t{otu}");
        }
        out.push('\n');
        for sample in &self.samples {
            let _ = write!(out, "vsearch\t{sample}\t{}", self.otus.len());
            for otu in &self.otus {
                let _ = write!(out, "\t{}", self.count(otu, sample));
            }
            out.push('\n');
        }
        out
    }

    /// vsearch's `--biomout`, a sparse BIOM 1.0 table, with its `id`, the
    /// program that `generated_by` it, and its `date`.
    pub fn biom(&self, id: &str, generated_by: &str, date: &str) -> String {
        let mut out = format!(
            "{{\n\t\"id\":\"{id}\",\n\t\"format\": \"Biological Observation Matrix 1.0\",\n\t\
             \"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n\t\
             \"type\": \"OTU table\",\n\t\"generated_by\": \"{generated_by}\",\n\t\"date\": \"{date}\",\n\t\
             \"matrix_type\": \"sparse\",\n\t\"matrix_element_type\": \"int\",\n\t\"shape\": [{},{}],\n",
            self.otus.len(),
            self.samples.len()
        );
        out.push_str("\t\"rows\":[");
        for (i, otu) in self.otus.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "\n\t\t{{\"id\":\"{otu}\", \"metadata\":");
            if self.taxonomy.is_empty() {
                out.push_str("null");
            } else {
                let _ = write!(
                    out,
                    "{{\"taxonomy\":\"{}\"}}",
                    self.taxonomy.get(otu).map_or("", String::as_str)
                );
            }
            out.push('}');
        }
        out.push_str("\n\t],\n\t\"columns\":[");
        for (i, sample) in self.samples.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "\n\t\t{{\"id\":\"{sample}\", \"metadata\":null}}");
        }
        out.push_str("\n\t],\n\t\"data\": [");
        let otu_no: BTreeMap<&str, usize> = self
            .otus
            .iter()
            .enumerate()
            .map(|(i, o)| (o.as_str(), i))
            .collect();
        let sample_no: BTreeMap<&str, usize> = self
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        for (i, ((otu, sample), count)) in self.counts.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "\n\t\t[{},{},{count}]",
                otu_no[otu.as_str()],
                sample_no[sample.as_str()]
            );
        }
        out.push_str("\n\t]\n}\n");
        out
    }
}
