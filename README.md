# Insize

Insert size stats and plot.

## Getting started

Help message

```shell
insize 0.1.0
slyo <sean.lyo@outlook.com>
Fast insert size distribution plot from bam.

USAGE:
    insize [OPTIONS] <bam> -o <FILE>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o <FILE>          Output pic file path, support `.svg` and `.png` suffix.
    -m <NUMBER>        Maximum insert size to record, default 500. !Bigger number costs more memory!.

ARGS:
    <bam>    Input bam file.

```

Run test case with PNG output.

```
insize -o insert-size.png test.bam
```

SVG is also supported.

```shell
insize -o insert-size.svg tests/test.bam
```

## Benchmark

~ 20s/Gb
