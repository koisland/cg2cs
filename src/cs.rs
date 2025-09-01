use std::error::Error;

use crate::{CigarOp, cg::Kind, cs_str_to_cs_ops};

/// `cs` tag operations
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum CSKind {
    /// Identical sequence (short and long form)
    Match,
    /// Substitution: ref to query
    Mismatch,
    /// Insertion to the reference
    Insertion,
    /// Deletion from the reference
    Deletion,
    /// Intron length and splice signal
    /// * Not implemented.
    Intron,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CSOp {
    pub(crate) kind: CSKind,
    pub(crate) len: usize,
    pub(crate) seq: Option<String>,
}

impl CSOp {
    pub(crate) fn new(kind: CSKind, len: usize, seq: Option<&str>) -> Self {
        Self {
            kind,
            len,
            seq: seq.map(|s| s.to_owned()),
        }
    }

    /// `cs` operation kind
    pub fn kind(&self) -> CSKind {
        self.kind
    }

    /// `cs` operation length
    pub fn length(&self) -> usize {
        self.len
    }

    /// `cs` operation associated sequence elements.
    /// * From [`minimap2`](https://lh3.github.io/minimap2/minimap2.html) docs.
    ///
    /// |Op|Regex|Description|
    /// |-|-|-|
    /// |=|\[ACGTN\]+|Identical sequence (long form)|
    /// |*|\[acgtn\]\[acgtn\]|Substitution: ref to query|
    /// |+|\[acgtn\]+|Insertion to the reference|
    /// |-|\[acgtn\]+|Deletion from the reference|
    pub fn seq(&self) -> Option<&str> {
        self.seq.as_deref()
    }
}

impl From<CSOp> for String {
    fn from(op: CSOp) -> Self {
        match op.kind {
            CSKind::Match => {
                if let Some(seq) = &op.seq {
                    format!("={seq}")
                } else {
                    format!(":{}", op.len)
                }
            }
            CSKind::Mismatch => {
                let seq = op.seq.expect("Missing sequence for mismatch.");
                format!("*{seq}")
            }
            CSKind::Deletion => {
                let seq = op.seq.expect("Missing sequence for mismatch.");
                format!("-{}", seq)
            }
            CSKind::Insertion => {
                let seq = op.seq.expect("Missing sequence for mismatch.");
                format!("+{}", seq)
            }
            CSKind::Intron => {
                unimplemented!("{}", format!("CSOp to String non implemented: {:?}", op))
            }
        }
    }
}

impl From<CSOp> for CigarOp {
    fn from(op: CSOp) -> Self {
        let cigar_kind = match op.kind {
            CSKind::Match => Kind::SequenceMatch,
            CSKind::Mismatch => Kind::SequenceMismatch,
            CSKind::Insertion => Kind::Insertion,
            CSKind::Deletion => Kind::Deletion,
            CSKind::Intron => unimplemented!("CSOp to CigarOp for Intron not implemented."),
        };
        CigarOp::new(cigar_kind, op.len)
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub(crate) enum CSToken {
    Identical,
    IdenticalLong,
    Substitution,
    Insertion,
    Deletion,
    Intron,
    Base,
    Number,
}

impl TryFrom<u8> for CSToken {
    type Error = Box<dyn Error>;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Ok(match value {
            b'=' => CSToken::IdenticalLong,
            b':' => CSToken::Identical,
            b'*' => CSToken::Substitution,
            b'+' => CSToken::Insertion,
            b'-' => CSToken::Deletion,
            b'~' => CSToken::Intron,
            b'0'..=b'9' => CSToken::Number,
            b'a' | b't' | b'g' | b'c' | b'n' | b'A' | b'T' | b'G' | b'C' | b'N' => CSToken::Base,
            _ => Err(format!("Invalid character {}", value as char))?,
        })
    }
}

impl TryFrom<CSToken> for CSKind {
    type Error = Box<dyn Error>;

    fn try_from(value: CSToken) -> Result<Self, Self::Error> {
        Ok(match &value {
            CSToken::Identical | CSToken::IdenticalLong => CSKind::Match,
            CSToken::Substitution => CSKind::Mismatch,
            CSToken::Insertion => CSKind::Insertion,
            CSToken::Deletion => CSKind::Deletion,
            CSToken::Intron => CSKind::Intron,
            _r => Err(format!("Invalid token to op {value:?}"))?,
        })
    }
}

impl TryFrom<CSToken> for Kind {
    type Error = Box<dyn Error>;

    fn try_from(value: CSToken) -> Result<Self, Self::Error> {
        Ok(match value {
            CSToken::Identical => Kind::SequenceMatch,
            CSToken::Substitution => Kind::SequenceMismatch,
            CSToken::Insertion => Kind::Insertion,
            CSToken::Deletion => Kind::Deletion,
            _ => Err(format!(
                "Cannot convert cs token {value:?} to cg operation."
            ))?,
        })
    }
}

/// Difference string (`cs`) tag.
/// * See https://lh3.github.io/minimap2/minimap2.html
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CS {
    repr: String,
    ops: Vec<CSOp>,
}

impl CS {
    /// Create `cs` tag from string
    /// ```
    /// use cg2cs::CS;
    /// assert!(CS::new(":2*cg-aa:3").is_ok());
    /// ```
    pub fn new(cs: &str) -> Result<Self, Box<dyn Error>> {
        cs_str_to_cs_ops(cs).map(CS::from)
    }

    /// Get `cs` string representation.
    /// ```
    /// use cg2cs::{cs_str_to_cs_ops, CS};
    /// let cs = CS::new(":2*cg-aa:3").unwrap();
    /// assert_eq!(
    ///     cs.repr(),
    ///     ":2*cg-aa:3"
    /// )
    /// ```
    pub fn repr(&self) -> &str {
        &self.repr
    }

    /// Get `cs` operations.
    /// ```
    /// use cg2cs::{cg_to_cs, CSKind, CS};
    /// let cs = CS::new(":2*cg-aa:3").unwrap();
    /// assert_eq!(cs.ops().len(), 4)
    /// ```
    pub fn ops(&self) -> &[CSOp] {
        &self.ops
    }
    /// Get `cs` as tag.
    /// ```
    /// use cg2cs::CS;
    /// let cs = CS::new(":2*cg-aa:3").unwrap();
    /// assert_eq!(cs.tag(), "cs:Z::2*cg-aa:3")
    /// ```
    pub fn tag(&self) -> String {
        format!("cs:Z:{}", self.repr)
    }
}

impl From<Vec<CSOp>> for CS {
    fn from(ops: Vec<CSOp>) -> Self {
        let mut repr = String::new();
        for op in ops.iter() {
            repr.push_str(&Into::<String>::into(op.clone()));
        }
        Self { repr, ops }
    }
}
