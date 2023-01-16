#!/usr/bin/env python3

import os, subprocess, struct

path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/'


###############################################################################
# Executable calls
###############################################################################

def run_compressor(
    input_file: str, output_file: str,
    log_file: str, binlog_file: str, dump_file: str,
    technique: str,
    n_bits_start: int, n_exponent_bits: int, flat_quantization: bool,
    compression_ac: float, compression_start_dc: float, flat_compression: bool,
    downsampling_ratio_ac: int = 1,
    effort: int = 7):

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
        '-r', str(n_exponent_bits),
        '-a', str(compression_ac),
        '-b', str(compression_start_dc),
        '-e', str(effort),
        '-s', str(downsampling_ratio_ac)
        ]

    if flat_quantization:
        args.append('--q_flat')

    if flat_compression:
        args.append('--c_flat')

    subprocess.run(args)


def run_decompressor(input_file: str, output_file: str):
    if os.path.exists(output_file):
        return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'decompress'),
        input_file, output_file
        ]

    subprocess.run(args)

    ' '.join(args)


def run_converter_exr_png(input_file: str, output_file: str, exposure: float):
    if os.path.exists(output_file):
        return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'exr-png'),
        input_file, output_file,
        '-e', str(exposure)
        ]

    subprocess.run(args)


def run_diff(file_a: str, file_b: str, max_err: float, output_file: str, diff_error_file: str):
    if os.path.exists(output_file):
        return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'exr-diff'),
        file_a, file_b,
        output_file,
        '-u', str(max_err),
        '-e', diff_error_file
        ]

    subprocess.run(args)


###############################################################################
# Get various paths for I/O
###############################################################################

def get_path_cave_in(folder: str, dataset_name: str):
    return os.path.join(folder, dataset_name + '.exr')


def get_path_cave_out(
    start: str, dataset_name: str,
    technique: str,
    n_bits_start: int,
    c_dc: float, c_ac: float,
    flat_q: bool, flat_c: bool):

    return os.path.join(
        start,
        dataset_name,
        technique,
        str(n_bits_start),
        str(c_dc), str(c_ac),
        'q_flat' if flat_q else 'q_dynamic',
        'c_flat' if flat_c else 'c_dynamic'
    )


def get_path_bonn_in(folder: str, dataset_name: str, type: str):
    return os.path.join(
        folder,
        dataset_name,
        type + '.exr'
    )


def get_path_bonn_out(
    start: str, dataset_name: str, type: str,
    technique: str,
    n_bits_start: int,
    c_dc: float, c_ac: float,
    flat_q: bool, flat_c: bool):

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

def get_err_from_bin_log(path: str):
    with open(path, 'rb') as f:
        data = f.read()
        # n = struct.unpack_from('i', data)
        err = struct.unpack_from('8d', data, offset=struct.calcsize('I'))

        return err # RMSE


def get_err_from_diff_bin(path: str):
    with open(path, 'rb') as f:
        data = f.read()
        err = struct.unpack_from('d', data)

        return err[0]

# Get quantization curve saved in the text log
def get_q_curve_from_txt_log(path: str):
    q_curve = []
    c_curve = []
    rmse = 0

    with open(path, 'r') as f:
        lines = f.readlines()
        for n in lines[2].split():
            q_curve.append(int(n))

    return q_curve


# Get compression curve saved in the text log
def get_c_curve_from_txt_log(path: str):
    c_curve = []

    with open(path, 'r') as f:
        lines = f.readlines()
        for n in lines[9].split():
            c_curve.append(float(n))

    return c_curve


def get_jxl_dir_size(path: str):
    size = 0

    for f in os.listdir(path):
        if f.endswith('.jxl'):
            size += os.path.getsize(os.path.join(path, f))

    return size
