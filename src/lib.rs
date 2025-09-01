use itertools::Itertools;
use std::error::Error;

mod cg;
mod cs;

use cg::CGToken;
use cs::CSToken;

pub use cg::{Cigar, CigarOp, Kind};
pub use cs::{CS, CSKind, CSOp};

/// Convert a cs tag into a Cigar.
///
/// # Args
/// * `cs`
///     * `cs` tag without tag prefix.
/// # Example
/// ```
/// use cg2cs::cs_to_cg;
///
/// let res = cs_to_cg(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
/// assert_eq!(
///     res.repr(),
///     "10=5=5I5D1X1=".to_string()
/// )
/// ```
pub fn cs_to_cg(cs: &str) -> Result<Cigar, Box<dyn Error>> {
    let cg_ops: Vec<CigarOp> = cs_str_to_cs_ops(cs)?
        .into_iter()
        .map(|cs_op| CigarOp::from(cs_op))
        .collect();
    Ok(Cigar::from(cg_ops))
}

/// Convert `cs` string to `cs` operations.
///
/// # Args
/// * `cs`
///     * `cs` string without tag prefix.
///
/// # Returns
/// * Cigar operations as [`CigarOp`].
///
/// # Examples
/// ```
/// use cg2cs::{cs_str_to_cs_ops, CSKind};
///
/// let ops = cs_str_to_cs_ops(":2*cg-aa:3").unwrap();
/// let ops: Vec<(CSKind, usize, Option<String>)> = ops
///     .iter()
///     .map(|op| (op.kind(), op.len(), op.seq().map(|s| s.to_owned())))
///     .collect();
/// let exp = vec![
///     (CSKind::Match, 2, None),
///     (CSKind::Mismatch, 1, Some(String::from("cg"))),
///     (CSKind::Deletion, 2, Some(String::from("aa"))),
///     (CSKind::Match, 3, None),
/// ];
/// assert_eq!(ops, exp)
/// ```
pub fn cs_str_to_cs_ops(cs: &str) -> Result<Vec<CSOp>, Box<dyn Error>> {
    let cs = cs.as_bytes();
    let mut new_ops = vec![];
    let mut curr_op: Option<CSToken> = None;
    // TODO: If intron added, convert to queue.
    for (tk, elems) in &cs
        .iter()
        .chunk_by(|c| TryInto::<CSToken>::try_into(**c).unwrap())
    {
        let elems: Vec<&u8> = elems.collect();

        match (curr_op.as_ref(), &tk) {
            // Consume first token.
            (None, CSToken::Identical)
            | (None, CSToken::IdenticalLong)
            | (None, CSToken::Substitution)
            | (None, CSToken::Insertion)
            | (None, CSToken::Deletion)
            | (None, CSToken::Intron) => curr_op = Some(tk.clone()),
            (None, _) => Err(format!("Invalid starting token: {tk:?}"))?,
            // Then fill in by type.

            // Identical
            (Some(CSToken::Identical), CSToken::Number) => {
                new_ops.push(CSOp {
                    kind: CSKind::Match,
                    len: elems.into_iter().map(|e| char::from(*e)).join("").parse()?,
                    seq: None,
                });
                curr_op.take();
            }
            // Identical long
            (Some(CSToken::IdenticalLong), CSToken::Base) => {
                new_ops.push(CSOp {
                    kind: CSKind::Match,
                    len: elems.len(),
                    seq: Some(elems.into_iter().map(|e| char::from(*e)).join("")),
                });
                curr_op.take();
            }
            // Substitution/mismatch
            (Some(CSToken::Substitution), CSToken::Base) => {
                new_ops.push(CSOp {
                    kind: CSKind::Mismatch,
                    // Only one base affected.
                    len: 1,
                    seq: Some(elems.into_iter().map(|e| char::from(*e)).join("")),
                });
                curr_op.take();
            }
            // Insertion/deletion
            (Some(op), CSToken::Base) => {
                let op_kind: CSKind = op.clone().try_into()?;
                new_ops.push(CSOp {
                    kind: op_kind,
                    len: elems.len(),
                    seq: Some(elems.into_iter().map(|e| char::from(*e)).join("")),
                });
                curr_op.take();
            }
            (Some(CSToken::Intron), _) => {
                unimplemented!("Not implemented.")
            }
            (Some(op), CSToken::Identical)
            | (Some(op), CSToken::IdenticalLong)
            | (Some(op), CSToken::Substitution)
            | (Some(op), CSToken::Insertion)
            | (Some(op), CSToken::Deletion)
            | (Some(op), CSToken::Number)
            | (Some(op), CSToken::Intron) => Err(format!(
                "Invalid matching token: {op:?}: {}",
                elems.into_iter().map(|e| *e as char).join("")
            ))?,
        }
    }
    Ok(new_ops)
}

/// Convert cigar string to cigar operations.
///
/// # Args
/// * `cg`
///     * Cigar string.
///
/// # Returns
/// * Cigar operations as [`CigarOp`].
///
/// # Examples
/// ```
/// use cg2cs::{cg_str_to_cg_ops, Kind};
///
/// let res = cg_str_to_cg_ops("2=1X2D3=").unwrap();
/// let res = res.iter().map(|op| (op.kind(), op.len())).collect::<Vec<(Kind, usize)>>();
/// let exp: Vec<(Kind, usize)> = vec![
///     (Kind::SequenceMatch, 2),
///     (Kind::SequenceMismatch, 1),
///     (Kind::Deletion, 2),
///     (Kind::SequenceMatch, 3)
/// ];
/// assert_eq!(res, exp)
/// ```
pub fn cg_str_to_cg_ops(cg: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    let mut ops = vec![];
    let mut prev_num: Option<usize> = None;
    for (kind, elems) in &cg
        .chars()
        .chunk_by(|c| TryInto::<CGToken>::try_into(*c as u8).unwrap())
    {
        match (prev_num.as_mut(), &kind) {
            (None, CGToken::Number) => {
                prev_num = Some(elems.into_iter().collect::<String>().parse()?)
            }
            (None, CGToken::Kind(kind)) => {
                return Err(format!("Invalid starting token ({kind:?})"))?;
            }
            (Some(number), CGToken::Kind(kind)) => {
                ops.push(CigarOp::new(*kind, *number));
                prev_num.take();
            }
            (Some(number), CGToken::Number) => {
                let curr_num: usize = elems.into_iter().collect::<String>().parse()?;
                return Err(format!(
                    "Invalid number followed by number ({number:?}, {curr_num:?})"
                ))?;
            }
        }
    }

    Ok(ops)
}


/// Convert a cigar string to a difference (`cs`) string.
/// * Adapted from [`minimap2::write_cs_ds_core`](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/format.c#L171)
/// * Intron op and [`ds`](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/NEWS.md?plain=1#L99) tag not supported.
///
/// # Args
/// * `cg`
///     * Cigar string.
/// * `tseq`
///     * Template sequence.
/// * `qseq`
///     * Query sequence.
/// * `no_iden`
///     * Use `:` intead of `=`. See [here](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/minimap.h#L406).
///
/// # Returns
/// * [`CS`]
///
/// # Example
/// ```
/// use cg2cs::{cg_to_cs, CS, CSOp, CSKind};
///
/// let res = cg_to_cs("2=1X2D3=", "ATCAATTT", "ATGTTT", true).unwrap();
/// assert_eq!(
///     res.repr(),
///     ":2*cg-aa:3"
/// )
/// ```
pub fn cg_to_cs(cg: &str, tseq: &str, qseq: &str, no_iden: bool) -> Result<CS, Box<dyn Error>> {
    let (mut q_off, mut t_off) = (0, 0);
    // let (mut q_len, mut t_len) = (0, 0);
    let cg_ops = cg_str_to_cg_ops(cg)?;
    let qseq = qseq.as_bytes();
    let tseq = tseq.as_bytes();

    // // Only required for intron and ds tag
    // for cg_op in cg_ops.iter() {
    //     let (kind, len) = (cg_op.kind, cg_op.len);
    //     if matches!(
    //         kind,
    //         Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
    //     ) {
    //         q_len += len;
    //         t_len += len;
    //     } else if kind == Kind::Insertion {
    //         q_len += len;
    //     } else if matches!(kind, Kind::Deletion | Kind::Skip) {
    //         t_len += len;
    //     }
    // }
    let mut ops = vec![];
    let mut tmp_seq = String::new();

    for cg_op in cg_ops.iter() {
        let (kind, len) = (cg_op.kind, cg_op.len);
        tmp_seq.clear();
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for j in 0..len {
                    let Some(q_chr) = qseq.get(q_off + j).map(|q| char::from(*q)) else {
                        Err(format!(
                            "Cannot extract query sequence element for {cg_op:?} at position {}",
                            q_off + j
                        ))?
                    };
                    let Some(t_chr) = tseq.get(t_off + j).map(|t| char::from(*t)) else {
                        Err(format!(
                            "Cannot extract template sequence element for {cg_op:?} at position {}",
                            t_off + j
                        ))?
                    };
                    if q_chr != t_chr {
                        if !tmp_seq.is_empty() {
                            if !no_iden {
                                // Original set index back. tmp[l_tmp] = 0
                                ops.push(CSOp::new(CSKind::Match, tmp_seq.len(), Some(&tmp_seq)));
                            } else {
                                ops.push(CSOp::new(CSKind::Match, tmp_seq.len(), None));
                            }
                            // Don't set position and instead clear. l_tmp = 0
                            tmp_seq.clear();
                        }
                        // Add mismatches
                        ops.push(CSOp {
                            kind: CSKind::Mismatch,
                            len: 1,
                            seq: Some(format!(
                                "{}{}",
                                t_chr.to_ascii_lowercase(),
                                q_chr.to_ascii_lowercase(),
                            )),
                        })
                    } else {
                        tmp_seq.push(q_chr);
                    }
                }
                if !tmp_seq.is_empty() {
                    if !no_iden {
                        ops.push(CSOp::new(CSKind::Match, tmp_seq.len(), Some(&tmp_seq)));
                    } else {
                        ops.push(CSOp::new(CSKind::Match, tmp_seq.len(), None));
                    }
                    tmp_seq.clear();
                }
                q_off += len;
                t_off += len;
            }
            Kind::Insertion => {
                for j in 0..len {
                    let Some(q_chr) = qseq.get(q_off + j).map(|q| char::from(*q)) else {
                        Err(format!(
                            "Cannot extract query sequence element for {cg_op:?} at position {}",
                            q_off + j
                        ))?
                    };
                    tmp_seq.push(q_chr.to_ascii_lowercase());
                }
                ops.push(CSOp::new(CSKind::Insertion, tmp_seq.len(), Some(&tmp_seq)));
                q_off += len;
            }
            Kind::Deletion => {
                for j in 0..len {
                    let Some(t_chr) = tseq.get(t_off + j).map(|t| char::from(*t)) else {
                        Err(format!(
                            "Cannot extract template sequence element for {cg_op:?} at position {}",
                            t_off + j
                        ))?
                    };
                    tmp_seq.push(t_chr.to_ascii_lowercase());
                }
                ops.push(CSOp::new(CSKind::Deletion, tmp_seq.len(), Some(&tmp_seq)));
                t_off += len;
            }
            _ => unimplemented!("Not implemented: {cg_op:?}"),
        }
    }

    Ok(CS::from(ops))
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn test_cg_to_ops() {
        let cg = "10M4D100I1102=";
        let ops = cg_str_to_cg_ops(cg).unwrap();
        assert_eq!(
            ops,
            [
                CigarOp {
                    kind: Kind::Match,
                    len: 10,
                },
                CigarOp {
                    kind: Kind::Deletion,
                    len: 4,
                },
                CigarOp {
                    kind: Kind::Insertion,
                    len: 100,
                },
                CigarOp {
                    kind: Kind::SequenceMatch,
                    len: 1102,
                },
            ]
        )
    }

    #[test]
    fn test_cg_to_cs_no_iden() {
        let cg = "6=1X43=1X4=1X33=1X17=1X17=1D4=1X10=1X4=1X3=1X20=1X30=1X1=1X22=1X41=1X20=";
        let exp_cs_string = ":6*ca:43*ga:4*tc:33*ga:17*tc:17-a:4*ac:10*ga:4*ag:3*gc:20*ta:30*ct:1*ct:22*ct:41*ga:20";
        let tseq = "GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCTGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTAGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAA";
        let qseq = "GCCGGGAGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGACGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGACTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAACTTAGCCGGGCATGGTGGCGCGCGCCTGTAGTCCCAGCTACACGGGAGGCTGAGGCAGGAGAATGGCGTGAATCTGGGAGGCGGAGCTTGCAGTGAGTCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAAACTCCGTCTCAAAAAAAAAA";
        let cs = cg_to_cs(cg, tseq, qseq, true).unwrap();

        assert_eq!(exp_cs_string, cs.repr())
    }

    #[test]
    fn test_cg_to_cs_iden() {
        let cg = "6=1X43=1X4=1X33=1X17=1X17=1D4=1X10=1X4=1X3=1X20=1X30=1X1=1X22=1X41=1X20=";
        let exp_cs_string = "=GCCGGG*ca=GCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAG*ga=CGGG*tc=GGATCACGAGGTCAGGAGATCGAGACCATCCTG*ga=CTAACACGGTGAAACCC*tc=GTCTCTACTAAAAATAC-a=AAAA*ac=TTAGCCGGGC*ga=TGGT*ag=GCG*gc=GCGCCTGTAGTCCCAGCTAC*ta=CGGGAGGCTGAGGCAGGAGAATGGCGTGAA*ct=C*ct=GGGAGGCGGAGCTTGCAGTGAG*ct=CGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGA*ga=ACTCCGTCTCAAAAAAAAAA";
        let tseq = "GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCTGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTAGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAA";
        let qseq = "GCCGGGAGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGACGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGACTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAACTTAGCCGGGCATGGTGGCGCGCGCCTGTAGTCCCAGCTACACGGGAGGCTGAGGCAGGAGAATGGCGTGAATCTGGGAGGCGGAGCTTGCAGTGAGTCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAAACTCCGTCTCAAAAAAAAAA";
        let cs = cg_to_cs(cg, tseq, qseq, false).unwrap();

        assert_eq!(exp_cs_string, cs.repr())
    }

    #[test]
    fn test_cs_to_cg() {
        let res = cs_to_cg(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
        assert_eq!(
            res,
            Cigar {
                repr: "10=5=5I5D1X1=".to_string(),
                ops: [
                    CigarOp {
                        kind: Kind::SequenceMatch,
                        len: 10
                    },
                    CigarOp {
                        kind: Kind::SequenceMatch,
                        len: 5
                    },
                    CigarOp {
                        kind: Kind::Insertion,
                        len: 5
                    },
                    CigarOp {
                        kind: Kind::Deletion,
                        len: 5
                    },
                    CigarOp {
                        kind: Kind::SequenceMismatch,
                        len: 1
                    },
                    CigarOp {
                        kind: Kind::SequenceMatch,
                        len: 1
                    }
                ]
                .to_vec()
            }
        )
    }
}
