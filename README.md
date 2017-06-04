# Ludomath
A Rust library containing various math & number functions,
primarily suited for 2D graphics & game programming.
As such, it has an emphasis on speed, although admittedly
doesn't currently do any explicit SIMD optimizations.

Includes functions, constants and traits for 2D vector math
& transformations with matrices, cheap random number generation,
as well as a few basic numeric functions.

* [Documentation](https://docs.rs/ludomath)    
* [`crates.io` page](https://crates.io/crates/ludomath)
## Basic usage

In `Cargo.toml`:
```rust
[depencencies]
ludomath = "1.1"
```

And in your Rust source:
```rust
extern crate ludomath;

use ludomath::vec2d::*;

fn main() {
    let point = Point::new(2.0, 3.0);
    println!("{:?}", point);
}
```

## License
Copyright (c) 2017 Erik Bivrin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
The License is included in the repository, named LICENSE.txt.
You may also obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an **"AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND**, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
