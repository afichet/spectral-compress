# Spectral Compress

This repository contains the code to compress and decompress spectral images. The various tools are in the folder `apps` while the common code is in `lib/common`.

An additional folder `supplemental` is independent and contains the scripts to generate the supplemental material.

## Build

To build the compression / decompression utilities, you need a C++ compiler, CMake and `libpng` installed on your system.

Ensure you have all submodules checked out:

```bash
git submodule update --init --recursive
```

Then you can build the executables:

```bash
mkdir build
cd build
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release
ninja
```

All the executables are placed in `build/bin`.

## Compression

The `compress` utility allows compressing a spectral OpenEXR to spectral JPEG-XL. It compresses an input file `input.exr` to an output set of files with:

```bash
compress input.exr compressed.jxl
```

You can specify compression levels, executing the utility with `-h` flag gives the full list of the possible parameters.

For example, `compress` shall be called with the following arguments to compress `input.exr` with
- 16 bits integer for the AC components (arguments `-q 16 -r 0`),
- a deterministic compression curve (argument `-c deterministic`),
- a frame distance for DC component of 0.5 and a starting frame distance of 2 for the first AC component (arguments `-a 0.5 -b 2`),
- a subsampling ratio of AC components of 1:2 (argument `-s 2`).

```bash
compress -q 16 -r 0 -c deterministic -a 0.5 -b 2 -s 2 input.exr compressed.jxl
```

## Decompression

To decompress the compressed file to spectral OpenEXR, you can run the `decompress` utility:

```bash
decompress compressed.jxl decompressed.exr
```
