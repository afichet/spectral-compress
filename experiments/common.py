#!/usr/bin/env python3

import os, subprocess

path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/'

def run_compressor(
    input_file: str, output_file: str,
    log_file: str, binlog_file: str, dump_file: str,
    technique: str,
    n_bits_start: int, flat_quantization: bool,
    compression_ac: float, compression_start_dc: float, flat_compression: bool):

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    fp, fn = os.path.split(log_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'compress'),
        input_file, output_file,
        '-l', log_file,
        '-k', binlog_file,
        '-d', dump_file,
        '-m', technique,
        '-q', str(n_bits_start),
        '-a', str(compression_ac),
        '-b', str(compression_start_dc),
        '-e', str(7)
        ]

    if flat_quantization:
        args.append('--q_flat')

    if flat_compression:
        args.append('--c_flat')

    subprocess.run(args)


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

    print(' '.join(args))


def get_path_cave_in(folder, dataset_name):
    return os.path.join(folder, dataset_name + '.exr')


def get_path_cave_out(
    start, dataset_name,
    technique,
    n_bits_start,
    c_dc, c_ac,
    flat_q, flat_c):
    return os.path.join(
        start,
        dataset_name,
        technique,
        str(n_bits_start),
        str(c_dc), str(c_ac),
        'q_flat' if flat_q else 'q_dynamic'
        'c_flat' if flat_c else 'c_dynamic'
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


###############################################################################
# Parse output files for report
###############################################################################

import struct

def get_err_from_bin_log(path):
    with open(path, 'rb') as f:
        data = f.read()
        # n = struct.unpack_from('i', data)
        err = struct.unpack_from('8d', data, offset=struct.calcsize('I'))

        return err # RMSE


def get_err_from_diff_bin(path):
    with open(path, 'rb') as f:
        data = f.read()
        err = struct.unpack_from('d', data)

        return err[0]

# Get quantization curve saved in the text log
def get_q_curve_from_txt_log(path):
    q_curve = []
    c_curve = []
    rmse = 0

    with open(path, 'r') as f:
        lines = f.readlines()
        for n in lines[2].split():
            q_curve.append(int(n))

    return q_curve


# Get compression curve saved in the text log
def get_c_curve_from_txt_log(path):
    c_curve = []

    with open(path, 'r') as f:
        lines = f.readlines()
        for n in lines[8].split():
            c_curve.append(float(n))

    return c_curve


def get_jxl_dir_size(path):
    size = 0

    for f in os.listdir(path):
        if f.endswith('.jxl'):
            size += os.path.getsize(os.path.join(path, f))

    return size
