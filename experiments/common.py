#!/usr/bin/env python3

import os, subprocess, struct
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/'

save_tex = True

if (save_tex):
    matplotlib.use("pgf")
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })


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
    start: str, downsampling_ratio_ac: int,
    dataset_name: str,
    technique: str,
    n_bits_start: int,
    c_dc: float, c_ac: float,
    flat_q: bool, flat_c: bool):

    return os.path.join(
        '{}_{}'.format(start, downsampling_ratio_ac),
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
    start: str, downsampling_ratio_ac: int,
    dataset_name: str, variant: str,
    technique: str,
    n_bits_start: int,
    c_dc: float, c_ac: float,
    flat_q: bool, flat_c: bool):

    return os.path.join(
        '{}_{}'.format(start, downsampling_ratio_ac),
        dataset_name,
        variant,
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


def get_duration_from_txt_log(path: str):
    duration = 0

    with open(path, 'r') as f:
        lines = f.readlines()
        duration = lines[18].split()[2]

    return float(duration)

def get_jxl_dir_size(path: str):
    size = 0

    for f in os.listdir(path):
        if f.endswith('.jxl'):
            size += os.path.getsize(os.path.join(path, f))

    return size


###############################################################################
# Plotting
###############################################################################

def plot_mode_curves_param(
    output_filename: str,
    stats: dict,
    technique: str,
    n_bits: int,
    downsampling_ratios_ac: list,
    framedistances: list,
    flat_compression: list,
    key: str,
    y_label: str):
    w = 0.8
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    x = np.arange(len(flat_compression))

    n_downsampling_ratio_ac = len(downsampling_ratios_ac)
    n_framedistances        = len(framedistances)
    n_el_per_group = n_downsampling_ratio_ac * n_framedistances
    x_offset = n_el_per_group / 2 - n_el_per_group

    for ratio, i in zip(downsampling_ratios_ac, range(n_downsampling_ratio_ac)):
        for (dc, ac), j  in zip(framedistances, range(n_framedistances)):
            y = [
                stats[ratio][technique][n_bits][dc][ac][True][True][key],
                stats[ratio][technique][n_bits][dc][ac][True][False][key],
            ]

            x_offset = (i * n_framedistances + j) - n_el_per_group / 2 + .5

            ax.bar(
                x + x_offset * w/n_el_per_group,
                y,
                width=w/n_el_per_group,
                label="dc = {}, ac = {}, chroma downsampling: 1:{}".format(dc, ac, ratio))

    ax.set_xticks(x, [
        'Flat',
        'Dynamic',
        ]
    )
    ax.set_xlabel('Framedistance (compression quality)')
    ax.set_ylabel(y_label)

    ax.legend()

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_mode_curve_error(output_filename:str, stats:dict, technique:str, n_bits:int, downsampling_ratio_ac:list, frame_distances:list, flat_compression: list):
    plot_mode_curves_param(output_filename, stats, technique, n_bits, downsampling_ratio_ac, frame_distances, flat_compression, 'error', 'Error')


def plot_mode_curve_size(output_filename:str, stats:dict, technique:str, n_bits:int, downsampling_ratio_ac:list, frame_distances:list, flat_compression: list):
   plot_mode_curves_param(output_filename, stats, technique, n_bits, downsampling_ratio_ac, frame_distances, flat_compression, 'size', 'File size')


def plot_mode_curve_ratio(output_filename:str, stats:dict, technique:str, n_bits:int, downsampling_ratio_ac:list, frame_distances:list, flat_compression: list):
   plot_mode_curves_param(output_filename, stats, technique, n_bits, downsampling_ratio_ac, frame_distances, flat_compression, 'ratio', 'Compression ratio')


def plot_mode_curve_duration(output_filename:str, stats:dict, technique:str, n_bits:int, downsampling_ratio_ac:list, frame_distances:list, flat_compression: list):
   plot_mode_curves_param(output_filename, stats, technique, n_bits, downsampling_ratio_ac, frame_distances, flat_compression, 'duration', 'Computation time (ms)')


def plot_q_curves(output_filename, stats, techniques, n_bits):
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    for tech in techniques:
        if tech == 'upperbound' or tech == 'twobounds':
            x = np.arange(len(stats[tech][n_bits][True][False]['q_curve'][1:-1])) + 1
            y =               stats[tech][n_bits][True][False]['q_curve'][1:-1]
        else:
            x = np.arange(len(stats[tech][n_bits][True][False]['q_curve'][1:])) + 1
            y =               stats[tech][n_bits][True][False]['q_curve'][1:]

        ax.plot(x, y, label=tech.replace('_', ' '))

    ax.legend()
    ax.set_title('Dynamic compression curves')
    ax.set_xlabel('Moment order')
    ax.set_ylabel('Number of bits')

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_c_curves(
    output_filename: str,
    stats: dict,
    technique: str,
    n_bits: int,
    downsampling_ratio_ac: list,
    framedistances: list):
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    for ratio in downsampling_ratio_ac:
        for (c_dc, c_ac) in framedistances:
            y = stats[ratio][technique][n_bits][c_dc][c_ac][True][False]['c_curve'][1:]
            x = np.arange(len(y)) + 1
            ax.plot(x, y, label='dc = {}, ac = {}, Chroma downsampling = 1:{}'.format(c_dc, c_ac, ratio))

    ax.legend()
    ax.set_xlabel('Moment order')
    ax.set_ylabel('Framedistance (compression parameter)')

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()
