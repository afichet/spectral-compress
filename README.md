# Spectral Compress

This repository contains the code to compress and decompress spectral images. The various utilies are in the folder `apps` while the common code is in `lib/common`.

An additional folder `supplemental` is independent and contains the scripts to generate the supplemental material.

## Build

To build the compression / decompression utilities, you need a C++ compiler, CMake and OpenEXR installed on your system.

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

The `compress` utility allows compressing a spectral OpenEXR to spectral JPEG-XL. It compress an input file `input.exr` to an output set of files with:

```bash
compress input.exr compressed.jxl
```

## Decompression

To decompress the compressed file to spectral OpenEXR, you can run the decompress utility:

```bash
decompress compressed.jxl decompressed.exr
```
