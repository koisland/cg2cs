use std::error::Error;

use crate::cg_str_to_cg_ops;

/// From noodles
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// An alignment match (`M`).
    Match,
    /// An insertion into the reference (`I`).
    Insertion,
    /// A deletion from the reference (`D`).
    Deletion,
    /// A skipped region from the reference (`N`).
    Skip,
    /// A soft clip (`S`).
    SoftClip,
    /// A hard clip (`H`).
    HardClip,
    /// Padding (`P`).
    Pad,
    /// A sequence match (`=`).
    SequenceMatch,
    /// A sequence mismatch (`X`).
    SequenceMismatch,
}

impl From<Kind> for char {
    fn from(value: Kind) -> Self {
        match value {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CigarOp {
    pub(crate) kind: Kind,
    pub(crate) len: usize,
}

impl CigarOp {
    pub(crate) fn new(kind: Kind, len: usize) -> Self {
        Self { kind, len }
    }

    pub fn kind(&self) -> Kind {
        self.kind
    }

    pub fn len(&self) -> usize {
        self.len
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub(crate) enum CGToken {
    Kind(Kind),
    Number,
}

impl TryFrom<u8> for CGToken {
    type Error = Box<dyn Error>;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Ok(match value {
            b'M' => CGToken::Kind(Kind::Match),
            b'I' => CGToken::Kind(Kind::Insertion),
            b'D' => CGToken::Kind(Kind::Deletion),
            b'N' => CGToken::Kind(Kind::Skip),
            b'S' => CGToken::Kind(Kind::SoftClip),
            b'H' => CGToken::Kind(Kind::HardClip),
            b'P' => CGToken::Kind(Kind::Pad),
            b'=' => CGToken::Kind(Kind::SequenceMatch),
            b'X' => CGToken::Kind(Kind::SequenceMismatch),
            b'0'..=b'9' => CGToken::Number,
            _ => Err(format!("Invalid character {}", value as char))?,
        })
    }
}

/// Cigar string.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Cigar {
    pub(crate) repr: String,
    pub(crate) ops: Vec<CigarOp>,
}

impl Cigar {
    /// Create cigar from string.
    /// ```
    /// use cg2cs::Cigar;
    /// 
    /// assert!(Cigar::new("2=1X2D3=").is_ok());
    /// ```
    pub fn new(cg: &str) -> Result<Self, Box<dyn Error>> {
        cg_str_to_cg_ops(cg).map(|ops| Cigar::from(ops))
    }

    /// Get cigar string representation.
    ///
    /// ```
    /// use cg2cs::Cigar;
    /// let cg_str = "2=1X2D3=";
    /// let cg = Cigar::new(cg_str).unwrap();
    /// assert_eq!(cg_str, cg.repr())
    /// ```
    pub fn repr(&self) -> &str {
        &self.repr
    }

    /// Get cigar operations.
    ///
    /// ```
    /// use cg2cs::{cg_str_to_cg_ops, Cigar};
    /// let ops = cg_str_to_cg_ops("2=1X2D3=").unwrap();
    /// let cg = Cigar::from(ops.to_vec());
    /// assert_eq!(&ops, cg.ops())
    /// ```
    pub fn ops(&self) -> &[CigarOp] {
        &self.ops
    }

    /// Get cigar as tag.
    ///
    /// ```
    /// use cg2cs::Cigar;
    /// let cg = Cigar::new("2=1X2D3=").unwrap();
    /// assert_eq!(cg.tag(), "cg:Z:2=1X2D3=")
    /// ```
    pub fn tag(&self) -> String {
        format!("cg:Z:{}", self.repr)
    }
}

impl From<Vec<CigarOp>> for Cigar {
    fn from(ops: Vec<CigarOp>) -> Self {
        let mut repr = String::new();
        for op in ops.iter() {
            let op_char: char = op.kind.into();
            repr.push_str(&format!("{}{op_char}", op.len));
        }
        Self { repr, ops }
    }
}
