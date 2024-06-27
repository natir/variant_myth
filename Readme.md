#  variant_myth 🧬 💻

![Test](https://github.com/natir/variant_myth/workflows/Test/badge.svg)
![Lints](https://github.com/natir/variant_myth/workflows/Lints/badge.svg)
![MSRV](https://github.com/natir/variant_myth/workflows/MSRV/badge.svg)
[![Documentation](https://github.com/natir/variant_myth/workflows/Documentation/badge.svg)](https://natir.github.io/variant_myth/variant_myth)
[![CodeCov](https://codecov.io/gh/natir/variant_myth/branch/main/graph/badge.svg)](https://codecov.io/gh/natir/variant_myth)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/variant_myth/blob/master/LICENSE)

A variant annotator.

## Rationale

Variant annotation take time, too much time. variant_myth try to solve this problem using modern and parallel tools.

To detect if a variant match an interval we use IntervalSet, information associated with intervals are stored in classic Swiss Table with AHash hash function.

## Installation

### From source

```bash
git clone https://github.com/natir/variant_myth.git
cd variant_myth
cargo install --path .
```

## Usage

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.74.
