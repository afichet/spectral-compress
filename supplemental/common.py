#!/usr/bin/env python3

import os, subprocess, struct
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

path_bin  = os.path.join('..', 'build', 'bin')

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
    compression_dc: float, compression_start_ac: float, compression_curve_type: str,
    downsampling_ratio_ac: int = 1,
    effort: int = 7):

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    fp, fn = os.path.split(log_file)
    os.makedirs(fp, exist_ok=True)

    if technique == 'simple':
        args = [os.path.join(path_bin, 'simple-compress'),
            input_file, output_file,
            '-d', str(compression_dc),
            '-l', log_file
        ]
    else:
        args = [os.path.join(path_bin, 'compress'),
            input_file, output_file,
            '-l', log_file,
            '-k', binlog_file,
            '-d', dump_file,
            '-m', technique,
            '-q', str(n_bits_start),
            '-r', str(n_exponent_bits),
            '-a', str(compression_dc),
            '-b', str(compression_start_ac),
            '-e', str(effort),
            '-s', str(downsampling_ratio_ac)
            ]

        if flat_quantization:
            args.append('--q_flat')

        compression_param = 'invalid'
        if compression_curve_type == 'c_flat':
            compression_param = 'flat'
        elif compression_curve_type == 'c_dynamic':
            compression_param = 'dynamic'
        elif compression_curve_type == 'c_deterministic':
            compression_param = 'deterministic'

        args.append('-c')
        args.append(compression_param)

    subprocess.run(args)


def run_decompressor(input_file: str, output_file: str, technique: str, log_file: str):
    if os.path.exists(output_file):
        return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    if technique == 'simple':
        exec_path = os.path.join(path_bin, 'simple-decompress')
    else:
        exec_path = os.path.join(path_bin, 'decompress')

    args = [exec_path, input_file, output_file, '-l', log_file]

    subprocess.run(args)


def run_sgeg_min_max(input_file: str, output_file: str):
    if os.path.exists(output_file):
        return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    subprocess.run([
        os.path.join(path_bin, 'sgeg-min-max'),
        input_file, output_file
    ])


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
    # print('Warning: this shall not be used in production!')
    # if os.path.exists(output_file):
    #     return

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [os.path.join(path_bin, 'exr-diff'),
        file_a, file_b,
        output_file,
        '-u', str(max_err),
        '-e', diff_error_file,
        # '-v',
        '-p', '0.005'
        ]

    subprocess.run(args)


def run_exr_change_compression(input_file: str, output_file: str, target_pixel_type: str, target_compression: str):
    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [
        os.path.join(path_bin, 'exr-change-compression'),
        input_file, output_file,
        '-t', target_pixel_type,
        '-m', target_compression
    ]

    subprocess.run(args)


def run_strip_rgb(input_file: str, output_file: str):
    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    args = [
        os.path.join(path_bin, 'exr-strip-rgb'),
        input_file, output_file
    ]

    subprocess.run(args)


###############################################################################
# Get various paths for I/O
###############################################################################

def get_path_cave_in(folder: str, dataset_name: str):
    return os.path.join(folder, dataset_name + '.exr')


def get_path_cave_out(
    start: str, subsampling_ratio_ac: int,
    dataset_name: str,
    technique: str,
    n_bits_start: int,
    c_dc: float, c_ac: float,
    q_flat: bool, c_type: str):

    return os.path.join(
        '{}_{}'.format(start, subsampling_ratio_ac),
        dataset_name,
        technique,
        str(n_bits_start),
        str(c_dc), str(c_ac),
        'q_flat' if q_flat else 'q_dynamic',
        c_type
    )


def get_path_cave_out_partial(
    start: str, downsampling_ratio_ac: int,
    dataset_name: str,
    technique: str):

    p = os.path.join(
        '{}_{}'.format(start, downsampling_ratio_ac),
        dataset_name,
        technique
    )

    os.makedirs(p, exist_ok=True)

    return p


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
    q_flat: bool, c_type: bool):

    return os.path.join(
        '{}_{}'.format(start, downsampling_ratio_ac),
        dataset_name,
        variant,
        technique,
        str(n_bits_start),
        str(c_dc), str(c_ac),
        'q_flat' if q_flat else 'q_dynamic',
        c_type
    )


def get_path_bonn_out_partial(
    start: str, downsampling_ratio_ac: int,
    dataset_name: str, variant: str, technique: str):

    p = os.path.join(
        '{}_{}'.format(start, downsampling_ratio_ac),
        dataset_name,
        variant,
        technique
    )

    os.makedirs(p, exist_ok=True)

    return p


###############################################################################
# Parse output files for report
###############################################################################

def get_err_from_bin_log(path: str):
    with open(path, 'rb') as f:
        data = f.read()
        # n = struct.unpack_from('i', data)
        err = struct.unpack_from('8d', data, offset=struct.calcsize('I'))

        return err # RMSE


def get_rmse_from_diff_bin(path: str):
    with open(path, 'rb') as f:
        data = f.read()
        err = struct.unpack_from('2d', data)

        return err[0]


def get_five_percentile_from_diff_bin(path: str):
    with open(path, 'rb') as f:
        data = f.read()
        err = struct.unpack_from('2d', data)

        return err[1]


def get_min_max_from(path:str):
    min_curve, max_curve = [], []

    with open(path, 'r') as f:
        lines = f.readlines()
        for l in lines:
            curr_min, curr_max = l.split()
            min_curve.append(float(curr_min))
            max_curve.append(float(curr_max))

    return min_curve, max_curve


# Get quantization curve saved in the text log
def get_q_curve_from_txt_log(path: str):
    q_curve = []

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


# The timing log is not at the same position depending on the compression
# utility used, so a parameter is added
def get_duration_from_txt_log(path: str, technique: str):
    duration = 0

    if technique == 'simple':
        idx = 0
    else:
        idx = 18

    with open(path, 'r') as f:
        lines = f.readlines()
        duration = lines[idx].split()[2]

    return float(duration) / 1000


def get_duration_decompression(path: str):
    duration = 0

    with open(path, 'r') as f:
        lines = f.readlines()
        duration = lines[0].split()[2]

    return float(duration) / 1000


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
    frame_distances_base: list,
    technique: str,
    n_bits: int,
    subsampling_ratios_ac: list,
    frame_distances: list,
    flat_compression: list,
    key: str,
    y_label: str,
    ratio: float,
    plot_w: float=15,
    plot_h: float=12):
    w = 0.9
    fig, ax = plt.subplots(1, 1, figsize=(ratio * plot_w, ratio * plot_h))

    x = np.arange(len(flat_compression) + 1)

    n_downsampling_ratio_ac = len(subsampling_ratios_ac)
    n_framedistances        = len(frame_distances)
    n_el_per_group = n_downsampling_ratio_ac * n_framedistances
    x_offset = n_el_per_group / 2 - n_el_per_group

    for ratio, i in zip(subsampling_ratios_ac, range(n_downsampling_ratio_ac)):
        for (dc, ac), j  in zip(frame_distances, range(n_framedistances)):
            y = [
                0, # Ugly harcoded: baseline placeholder
                stats[ratio][technique][n_bits][dc][ac][True]['c_flat'][key],
                stats[ratio][technique][n_bits][dc][ac][True]['c_deterministic'][key],
                stats[ratio][technique][n_bits][dc][ac][True]['c_dynamic'][key],
            ]

            x_offset = (i * n_framedistances + j) - n_el_per_group / 2 + .5

            ax.bar(
                x + x_offset * w/n_el_per_group,
                y,
                width=w/n_el_per_group,
                label='dc = {}, ac = {}, chroma subsampling: 1:{}'.format(dc, ac, ratio))

    # Ugly hardcoded: baseline values
    n_framedistances_base = len(frame_distances_base)

    for (f_distance, _), ij in zip(frame_distances_base, range(n_framedistances_base)):
        y = [
            stats[1]['simple'][32][f_distance][0][True]['c_flat'][key],
            0, 0, 0
        ]

        x_offset = ij - n_framedistances_base / 2 + .5

        ax.bar(
            x + x_offset * w/n_framedistances_base,
            y,
            width=w/n_framedistances_base,
            label='simple - frame distance {}'.format(f_distance)
        )

    ax.set_xticks(x, [
        'Simple', # Ugly hardcoded: label baseline
        'Flat',
        'Deterministic',
        'Dynamic'
        ]
    )
    # ax.set_xlabel('Distance level (compression parameter)')
    ax.set_ylabel(y_label)

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_rmse(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
    plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'rmse', 'RMSE', ratio, plot_w, plot_h)


def plot_size(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'size', 'File size', ratio, plot_w, plot_h)


def plot_compression_ratio(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'ratio', 'Compression ratio', ratio, plot_w, plot_h)


def plot_duration_compression(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'time_c', 'Compression time (s)', ratio, plot_w, plot_h)


def plot_duration_decompression(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'time_d', 'Decompression time (s)', ratio, plot_w, plot_h)


def plot_duration_compression_per_pixel(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'time_c', 'Compression time per pixel (ms)', ratio, plot_w, plot_h)


def plot_duration_decompression_per_pixel(output_filename:str, stats:dict, frame_distances_base: list, technique:str, n_bits:int, subsampling_ratios_ac:list, frame_distances:list, flat_compression: list, ratio: float, plot_w: float, plot_h: float):
   plot_mode_curves_param(output_filename, stats, frame_distances_base, technique, n_bits, subsampling_ratios_ac, frame_distances, flat_compression, 'time_d', 'Decompression time per pixel (ms)', ratio, plot_w, plot_h)


def plot_min_max(
    output_filename:str,
    stats:dict,
    frame_distances_base: list,
    technique:str,
    n_bits:int,
    subsampling_ratios_ac:list,
    frame_distances:list,
    flat_compression: list,
    ratio: float,
    plot_w: float=15,
    plot_h: float=12):
    fig, ax = plt.subplots(1, 1, figsize=(ratio * plot_w, ratio * plot_h))

    y_min = stats[1][technique][n_bits][0][1][True]['c_flat']['min_curve']
    y_max = stats[1][technique][n_bits][0][1][True]['c_flat']['max_curve']
    x_min = np.arange(len(y_min)) + 1
    x_max = np.arange(len(y_max)) + 1

    ax.bar(x_min, y_min, color='gray')
    ax.bar(x_max, y_max, color='gray')

    ax.set_ylabel('Value range')
    ax.set_xlabel('Moment order')

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_xy_ratio_error_curves_key(
    output_filename: str,
    stats: dict,
    dataset: list, variants: list,
    frame_distances_base: list,
    technique: str,
    n_bits: int,
    subsampling_ratios_ac: list,
    frame_distances: list,
    point_size: float,
    key_compression: str,
    title: str,
    x_max: float, y_max: float,
    ratio: float,
    plot_w: float=15,
    plot_h: float=12):
    fig, ax = plt.subplots(1, 1, figsize=(ratio * plot_w, ratio * plot_h))

    # for ratio in subsampling_ratios_ac:
    #     for (c_dc, c_ac) in framedistances:
    #         x = [v for v in stats[ratio][technique][c_dc][c_ac][True]['c_dynamic']['ratio']]
    #         y = [v for v in stats[ratio][technique][c_dc][c_ac][True]['c_dynamic']['rmse']]
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in frame_distances:
            x, y = [], []

            if len(variants) == 0:
                for d in dataset:
                    x.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True][key_compression]['ratio'])
                    y.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True][key_compression]['rmse'])
            else:
                for d in dataset:
                    for v in variants:
                        x.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True][key_compression]['ratio'])
                        y.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True][key_compression]['rmse'])
            ax.scatter(x, y, s=point_size, color=default_cols[idx])

            idx += 1

    ax.set_title(title)

    ax.set_xlabel('Compression ratio')
    ax.set_ylabel('RMSE')

    ax.set_xbound((0, x_max))
    ax.set_ybound((0, y_max))

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_xy_ratio_error_curves(
    output_filename: str,
    stats: dict,
    dataset: list, variants: list,
    frame_distances_base: list,
    technique: str,
    n_bits: int,
    subsampling_ratios_ac: list,
    frame_distances: list,
    point_size: float,
    ratio: float,
    plot_w: float=15,
    plot_h: float=12):
    fig, ax = plt.subplots(1, 1, figsize=(ratio * plot_w, ratio * plot_h))

    # for ratio in subsampling_ratios_ac:
    #     for (c_dc, c_ac) in framedistances:
    #         x = [v for v in stats[ratio][technique][c_dc][c_ac][True]['c_dynamic']['ratio']]
    #         y = [v for v in stats[ratio][technique][c_dc][c_ac][True]['c_dynamic']['rmse']]
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in frame_distances:
            x_flat, y_flat = [], []
            x_det, y_det = [], []
            x_dyn, y_dyn = [], []

            if len(variants) == 0:
                for d in dataset:
                    x_flat.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_flat']['ratio'])
                    y_flat.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_flat']['rmse'])

                    x_det.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_deterministic']['ratio'])
                    y_det.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_deterministic']['rmse'])

                    x_dyn.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_dynamic']['ratio'])
                    y_dyn.append(stats[d][ratio][technique][n_bits][c_dc][c_ac][True]['c_dynamic']['rmse'])
            else:
                for d in dataset:
                    for v in variants:
                        x_flat.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_flat']['ratio'])
                        y_flat.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_flat']['rmse'])

                        x_det.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_deterministic']['ratio'])
                        y_det.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_deterministic']['rmse'])

                        x_dyn.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_dynamic']['ratio'])
                        y_dyn.append(stats[d][v][ratio][technique][n_bits][c_dc][c_ac][True]['c_dynamic']['rmse'])
            ax.scatter(x_flat, y_flat, marker='s', s=point_size, color=default_cols[idx])
            ax.scatter(x_det, y_det, marker='^', s=point_size, color=default_cols[idx])
            ax.scatter(x_dyn, y_dyn, marker='X', s=point_size, color=default_cols[idx])

            idx += 1

    for (f, _) in frame_distances_base:
        x, y = [], []
        if len(variants) == 0:
            for d in dataset:
                x.append(stats[d][1]['simple'][32][f][0][True]['c_flat']['ratio'])
                y.append(stats[d][1]['simple'][32][f][0][True]['c_flat']['rmse'])
        else:
            for d in dataset:
                for v in variants:
                    x.append(stats[d][v][1]['simple'][32][f][0][True]['c_flat']['ratio'])
                    y.append(stats[d][v][1]['simple'][32][f][0][True]['c_flat']['rmse'])
        ax.scatter(x, y, marker='o', s=point_size, color=default_cols[idx])

        idx += 1

    ax.scatter([], [], marker='s', color='black', label='flat')
    ax.scatter([], [], marker='^', color='black', label='deterministic')
    ax.scatter([], [], marker='X', color='black', label='dynamic')

    ax.legend()

    ax.set_xlabel('Compression ratio')
    ax.set_ylabel('RMSE')

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()

# def plot_q_curves(output_filename, stats, techniques, n_bits):
#     fig, ax = plt.subplots(1, 1, figsize=(5, 4))

#     for tech in techniques:
#         if tech == 'upperbound' or tech == 'twobounds':
#             x = np.arange(len(stats[tech][n_bits][True][False]['q_curve'][1:-1])) + 1
#             y =               stats[tech][n_bits][True][False]['q_curve'][1:-1]
#         else:
#             x = np.arange(len(stats[tech][n_bits][True][False]['q_curve'][1:])) + 1
#             y =               stats[tech][n_bits][True][False]['q_curve'][1:]

#         ax.plot(x, y, label=tech.replace('_', ' '))

#     # ax.legend()
#     ax.set_title('Dynamic compression curves')
#     ax.set_xlabel('Moment order')
#     ax.set_ylabel('Number of bits')

#     fig.tight_layout()

#     if (save_tex):
#         plt.savefig(output_filename)
#         plt.close()
#     else:
#         plt.show()


def plot_c_curves(
    output_filename: str,
    stats: dict,
    technique: str,
    n_bits: int,
    subsampling_ratios_ac: list,
    framedistances: list,
    ratio: float,
    plot_w: float=15,
    plot_h: float=12):
    fig, ax = plt.subplots(1, 1, figsize=(ratio * plot_w, ratio * plot_h))
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in framedistances:
            y_dyn = stats[ratio][technique][n_bits][c_dc][c_ac][True]['c_dynamic']['c_curve'][1:]
            x = np.arange(len(y_dyn)) + 1
            ax.plot(x, y_dyn, color=default_cols[idx]) #, label='dc = {}, ac = {}, chroma subsampling = 1:{}'.format(c_dc, c_ac, ratio))

            idx += 1

    # WARN: Hardcoded to avoid color clash with the other plots...
    ax.plot([], [], color='black', label='dynamic')
    for (c_dc, c_ac), i in zip(framedistances, range(len(framedistances))):
        y_det = stats[1][technique][n_bits][c_dc][c_ac][True]['c_deterministic']['c_curve'][1:]
        x = np.arange(len(y_det)) + 1
        ax.plot(x, y_det, color=default_cols[-i - 1], linestyle='dotted', label='deterministic ac = {}'.format(c_ac))

    ax.legend(loc='lower right')
    ax.set_xlabel('Moment order')
    ax.set_ylabel('Distance level (compression parameter)')

    fig.tight_layout()
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=1, top=1, wspace=0, hspace=0)

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_legend_1(
    output_filename: str,
    frame_distances_base: list,
    subsampling_ratios_ac: list,
    frame_distances: list):
    # This function is hardcoded in seval places to have a nicer layout...
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots(1, 1)

    idx = len(subsampling_ratios_ac) * len(frame_distances)

    for (f, _) in frame_distances_base:
        ax.plot(
            [], [],
            marker='s', ls='none', color=default_cols[idx],
            label='simple - frame distance {}'.format(f))
        idx += 1

    # Dirty hack to align all other types on the same column
    ax.plot([], [], label=' ', color='white')

    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in frame_distances:
            ax.plot(
                [], [],
                marker='s', ls='none', color=default_cols[idx],
                label='chroma subsampling = 1:{}, dc = {}, ac = {}'.format(ratio, c_dc, c_ac))
            idx += 1

    legend = ax.legend(ncol=1, bbox_to_anchor=(1.05, 1), loc='upper left', fancybox=False)

    expand=[-5,-5,5,5]

    fig = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    if (save_tex):
        fig.savefig(output_filename, bbox_inches=bbox)
        plt.close()
    else:
        plt.show()


def plot_legend_2(
    output_filename: str,
    frame_distances_base: list,
    subsampling_ratios_ac: list,
    frame_distances: list):
    # This function is hardcoded in seval places to have a nicer layout...
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots(1, 1)

    idx = len(subsampling_ratios_ac) * len(frame_distances)

    for (f, _) in frame_distances_base:
        ax.plot(
            [], [],
            marker='s', ls='none', color=default_cols[idx],
            label='simple - frame distance {}'.format(f))
        idx += 1

    # Dirty hack to align all other types on the same column
    ax.plot([], [], label=' ', color='white')

    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in frame_distances:
            ax.plot(
                [], [],
                marker='s', ls='none', color=default_cols[idx],
                label='chroma subsampling = 1:{}, dc = {}, ac = {}'.format(ratio, c_dc, c_ac))
            idx += 1

    # Dirty hack to align all other types on the same column
    ax.plot([], [], label=' ', color='white')

    legend = ax.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left', fancybox=False)

    expand=[-5,-5,5,5]

    fig = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    if (save_tex):
        fig.savefig(output_filename, bbox_inches=bbox)
        plt.close()
    else:
        plt.show()


def plot_legend_2_strip(
    output_filename: str,
    frame_distances_base: list,
    subsampling_ratios_ac: list,
    frame_distances: list):
    # This function is hardcoded in seval places to have a nicer layout...
    default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots(1, 1)

    idx = 0

    for ratio in subsampling_ratios_ac:
        for (c_dc, c_ac) in frame_distances:
            ax.plot(
                [], [],
                marker='s', ls='none', color=default_cols[idx],
                label='chroma subsampling = 1:{}, dc = {}, ac = {}'.format(ratio, c_dc, c_ac))
            idx += 1

    legend = ax.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left', fancybox=False)

    expand=[-5,-5,5,5]

    fig = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    if (save_tex):
        fig.savefig(output_filename, bbox_inches=bbox)
        plt.close()
    else:
        plt.show()

def crop_png(png_filename_in: str, png_filename_out: str, cropped_size: int):
    img = Image.open(png_filename_in)

    left  = (img.width - cropped_size) / 2
    right = (img.width + cropped_size) / 2

    top    = (img.height - cropped_size) / 2
    bottom = (img.height + cropped_size) / 2

    img_cropped = img.crop((left, top, right, bottom))

    img_cropped.save(png_filename_out)
