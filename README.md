# `cg2cs`
Library for converting cigar strings to cs tags and vice-versa.

## Usage
```bash
cargo add --git https://github.com/koisland/cg2cs.git
```

### cs to cigar
```rust
use cg2cs::cs_to_cg;

let res = cs_to_cg(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
assert_eq!(
    res.repr(),
    "10=5=5I5D1X1="
)
```

### cigar to cs
```rust
use cg2cs::cg_to_cs;

// cg, template, query, use : or =
let res = cg_to_cs("2=1X2D3=", "ATCAATTT", "ATGTTT", true).unwrap();
assert_eq!(
    res.repr(),
    ":2*cg-aa:3"
)
```

## References
* https://lh3.github.io/minimap2/minimap2.html
* https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/format.c#L171
