# Variant Myth 🧬 💻

![tests](https://github.com/natir/variant_myth/workflows/tests/badge.svg)
![lints](https://github.com/natir/variant_myth/workflows/lints/badge.svg)
![msrv](https://github.com/natir/variant_myth/workflows/msrv/badge.svg)
[![website](https://github.com/natir/variant_myth/workflows/website/badge.svg)](https://natir.github.io/variant_myth)
[![dev-doc](https://github.com/natir/variant_myth/workflows/website/badge.svg)](https://natir.github.io/variant_myth/doc/variant_myth/)
[![license](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/variant_myth/blob/main/LICENSE)
[![copier](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-grayscale-inverted-border-orange.json)](https://github.com/copier-org/copier)

A fast genomic variant annotator.

## Rationale

Variant annotation take time, too much time, variant_myth try to solve this problem using modern and parallel tools.

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
