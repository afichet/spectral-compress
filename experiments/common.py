#!/usr/bin/env python3

import os, subprocess

path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/'

def run_compressor(
    input_file, output_file,
    log_file, binlog_file,
    technique,
    n_bits_start, flat_quantization,
    compression_start, flat_compression):

    for f in [output_file, log_file, binlog_file]:
        fp, fn = os.path.split(binlog_file)
        os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'compress'),
        input_file, output_file,
        '-l', log_file,
        '-k', binlog_file,
        '-m', technique,
        '-q', str(n_bits_start),
        '-a', str(0.1),
        '-b', str(compression_start),
        ]

    if flat_quantization:
        args.append('--q_flat')

    if flat_compression:
        args.append('--c_flat')

    subprocess.run(args)
    # print(' '.join(args))


def run_decompressor(input_file, output_file):
    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'decompress'),
        input_file, output_file
        ]

    subprocess.run(args)


def run_converter_exr_png(input_file, output_file, exposure):
    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'exr-png'),
        input_file, output_file,
        '-e', str(exposure)
        ]

    subprocess.run(args)


def run_diff(file_a, file_b, max_err, output_file, diff_error_file):
    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'exr-diff'),
        file_a, file_b,
        output_file,
        '-u', str(max_err),
        '-e', diff_error_file
        ]

    subprocess.run(args)


def get_path_cave(start, dataset_name, technique, n_bits_start, flat_curves):
    return os.path.join(
        start,
        dataset_name,
        technique,
        str(n_bits_start),
        'flat' if flat_curves else 'dynamic'
    )


def get_path_bonn_in(folder, dataset_name, type):
    return os.path.join(
        folder,
        dataset_name,
        type + '.exr'
    )

def get_path_bonn_out(
    start, dataset_name, type,
    technique,
    n_bits_start,
    c_dc, c_ac,
    flat_q, flat_c):
    return os.path.join(
        start,
        dataset_name,
        type,
        technique,
        str(n_bits_start),
        str(c_dc), str(c_ac),
        'q_flat' if flat_q else 'q_dynamic',
        'c_flat' if flat_c else 'c_dynamic'
    )
