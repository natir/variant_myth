#  variant_myth 🧬 💻

[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/variant_myth/blob/master/LICENSE)
![CI](https://github.com/natir/variant_myth/workflows/CI/badge.svg)
[![Documentation](https://github.com/natir/variant_myth/workflows/Documentation/badge.svg)](https://natir.github.io/variant_myth/variant_myth)
[![CodeCov](https://codecov.io/gh/natir/variant_myth/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/variant_myth)

A variant annotator.

## Rational

Variant annotation take time, too many time, variant_myth try to solve this problem by use modern and parallel tools.

To detect if a variant match interval we use IntervalSet, information associate with interval are store in classic Swiss Table with AHash hash function.

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
