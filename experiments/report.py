#!/usr/bin/env python3

from common import get_path_cave, run_decompressor, run_converter_exr_png, run_diff

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import os
import struct

cave_list = [
    'balloons_ms',
    # 'clay_ms',
    # 'fake_and_real_beers_ms',
    # 'fake_and_real_peppers_ms',
    # 'feathers_ms',
    # 'jelly_beans_ms',
    # 'pompoms_ms',
    # 'stuffed_toys_ms',
    # 'beads_ms',
    # 'cloth_ms',
    # 'fake_and_real_food_ms',
    # 'fake_and_real_strawberries_ms',
    # 'flowers_ms',
    # 'oil_painting_ms',
    # 'real_and_fake_apples_ms',
    # 'superballs_ms',
    # 'cd_ms',
    # 'egyptian_statue_ms',
    # 'fake_and_real_lemon_slices_ms',
    # 'fake_and_real_sushi_ms',
    # 'glass_tiles_ms',
    # 'paints_ms',
    # 'real_and_fake_peppers_ms',
    # 'thread_spools_ms',
    # 'chart_and_stuffed_toy_ms',
    # 'face_ms',
    # 'fake_and_real_lemons_ms',
    # 'fake_and_real_tomatoes_ms',
    # 'hairs_ms',
    # 'photo_and_face_ms',
    # 'sponges_ms'
]

techniques = ['linear', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits = [6, 8, 10, 12]

def main():
    # start_compress = [0.1]
    flat_quantization = [True, False]
    flat_compression  = [True, False]

    flat_curves = [True, False]

    prefix_cave = 'cave'

    matplotlib.use("pgf")
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })

    stats = {}
    g_stats = {}
    tex_stream = ''

    for tech in techniques:
        g_stats[tech] = {}

        g_stats[tech]['rmse_q'] = []
        g_stats[tech]['rmse_s'] = []
        g_stats[tech]['ratio']  = []


    for dirname in cave_list:
        stats [dirname] = {}

        path_export = os.path.join('export', prefix_cave)

        org_exr_file = get_original_dirname_cave(dirname)
        org_png_file = os.path.join(path_export, dirname + '.png')

        run_converter_exr_png(org_exr_file, org_png_file, -6.5)

        org_size = os.path.getsize(org_exr_file) / (1024*1024)

        for tech in techniques:
            stats [dirname][tech] = {}

            for flat in flat_curves:
                flat_name = 'flat' if flat else 'dyn'

                stats[dirname][tech][flat_name] = {}
                stats[dirname][tech][flat_name]['rmse_q'] = []
                stats[dirname][tech][flat_name]['rmse_c'] = []
                stats[dirname][tech][flat_name]['size'] = []
                stats[dirname][tech][flat_name]['ratio'] = []
                stats[dirname][tech][flat_name]['err_size'] = []

                for bits in start_bits:
                    path_curr_out = get_path_cave(prefix_cave, dirname, tech, bits, flat)
                    print(path_curr_out)

                    compressed_file = os.path.join(path_curr_out, dirname + '.jxl')
                    binlog_file     = os.path.join(path_curr_out, dirname + '.bin')
                    txtlog_file     = os.path.join(path_curr_out, dirname + '.txt')

                    size = get_jxl_dir_size(path_curr_out) / (1024*1024)
                    err = get_err_from_bin_log(binlog_file)
                    ratio = org_size / size

                    stats[dirname][tech][flat_name]['rmse_q'].append(err[0])
                    stats[dirname][tech][flat_name]['rmse_c'].append(err[3])
                    stats[dirname][tech][flat_name]['size'].append(size)
                    stats[dirname][tech][flat_name]['ratio'].append(ratio)
                    stats[dirname][tech][flat_name]['err_size'].append(err[3] / size)

                    path_curr_export = get_path_cave(path_export, dirname, tech, bits, flat)

                    decompressed_exr_file = os.path.join(path_curr_export, dirname + '.exr')
                    decompressed_png_file = os.path.join(path_curr_export, dirname + '.png')
                    diff_png_file         = os.path.join(path_curr_export, dirname + '_diff.png')

                    # run_decompressor(compressed_file, decompressed_exr_file)
                    # run_converter_exr_png(decompressed_exr_file, decompressed_png_file, -6.5)
                    # max_err = 7E-3
                    # run_diff(org_exr_file, decompressed_exr_file, max_err, diff_png_file)

        plot_err_q_file      = os.path.join(path_export, dirname + '_err_q.pgf')
        plot_err_c_file      = os.path.join(path_export, dirname + '_err_c.pgf')
        plot_ratio_file      = os.path.join(path_export, dirname + '_ratio.pgf')
        plot_err_size_file   = os.path.join(path_export, dirname + '_err_size.pgf')

        gen_plot_err_q(stats[dirname], plot_err_q_file)
        gen_plot_err_c(stats[dirname], plot_err_c_file)
        gen_plot_ratio(stats[dirname], plot_ratio_file)
        gen_plot_err_size(stats[dirname], plot_err_size_file)

        tex_stream += get_tex_stream(dirname)

    with open(os.path.join('export', 'cave.tex'), 'w') as f:
        f.write(tex_stream)


def get_tex_stream(dataset):
    with open(os.path.join('export', 'mat_template.tex'), 'r') as f:
        stream = ''.join(f.readlines())

        stream.replace('balloons_ms', dataset)

        return stream


def get_err_from_bin_log(path):
    with open(path, 'rb') as f:
        data = f.read()
        # n = struct.unpack_from('i', data)
        err = struct.unpack_from('6d', data, offset=struct.calcsize('I'))

        return err # RMSE


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
        for n in lines[7].split():
            c_curve.append(float(n))

    return c_curve


def get_jxl_dir_size(path):
    size = 0

    for f in os.listdir(path):
        if f.endswith('.jxl'):
            size += os.path.getsize(os.path.join(path, f))

    return size


def get_original_dirname_cave(dataset):
    return '/home/afichet/spectral_images/EXRs/CAVE/' + dataset + '.exr'




def gen_plot_err_q(data, output_file):
    width = 0.8
    curves = 'dyn'

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    ax.bar(np.arange(len(start_bits)) - 2 * width / 5, data['linear']              [curves]['rmse_q'], width=width / 5, label='linear')
    ax.bar(np.arange(len(start_bits)) - 1 * width / 5, data['unbounded']           [curves]['rmse_q'], width=width / 5, label='unbounded')
    ax.bar(np.arange(len(start_bits)) + 0 * width / 5, data['unbounded_to_bounded'][curves]['rmse_q'], width=width / 5, label='unbounded_to_bounded')
    ax.bar(np.arange(len(start_bits)) + 1 * width / 5, data['upperbound']          [curves]['rmse_q'], width=width / 5, label='upperbound')
    ax.bar(np.arange(len(start_bits)) + 2 * width / 5, data['twobounds']           [curves]['rmse_q'], width=width / 5, label='twobounds')

    ax.legend()

    ax.set_xticks(np.arange(len(start_bits)), start_bits)
    ax.set_xlabel('Start number of bits')
    ax.set_ylabel('RMSE')

    fig.tight_layout()
    plt.savefig(output_file)


def gen_plot_err_c(data, output_file):
    width = 0.8
    curves = 'dyn'

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    ax.bar(np.arange(len(start_bits)) - 2 * width / 5, data['linear']              [curves]['rmse_c'], width=width / 5, label='linear')
    ax.bar(np.arange(len(start_bits)) - 1 * width / 5, data['unbounded']           [curves]['rmse_c'], width=width / 5, label='unbounded')
    ax.bar(np.arange(len(start_bits)) + 0 * width / 5, data['unbounded_to_bounded'][curves]['rmse_c'], width=width / 5, label='unbounded_to_bounded')
    ax.bar(np.arange(len(start_bits)) + 1 * width / 5, data['upperbound']          [curves]['rmse_c'], width=width / 5, label='upperbound')
    ax.bar(np.arange(len(start_bits)) + 2 * width / 5, data['twobounds']           [curves]['rmse_c'], width=width / 5, label='twobounds')

    ax.legend()

    ax.set_xticks(np.arange(len(start_bits)), start_bits)
    ax.set_xlabel('Start number of bits')
    ax.set_ylabel('RMSE')

    fig.tight_layout()
    plt.savefig(output_file)


def gen_plot_ratio(data, output_file):
    width = 0.8
    curves = 'dyn'

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    ax.bar(np.arange(len(start_bits)) - 2 * width / 5, data['linear']              [curves]['ratio'], width=width / 5, label='linear')
    ax.bar(np.arange(len(start_bits)) - 1 * width / 5, data['unbounded']           [curves]['ratio'], width=width / 5, label='unbounded')
    ax.bar(np.arange(len(start_bits)) + 0 * width / 5, data['unbounded_to_bounded'][curves]['ratio'], width=width / 5, label='unbounded_to_bounded')
    ax.bar(np.arange(len(start_bits)) + 1 * width / 5, data['upperbound']          [curves]['ratio'], width=width / 5, label='upperbound')
    ax.bar(np.arange(len(start_bits)) + 2 * width / 5, data['twobounds']           [curves]['ratio'], width=width / 5, label='twobounds')

    ax.legend()

    ax.set_xticks(np.arange(len(start_bits)), start_bits)
    ax.set_xlabel('Start number of bits')
    ax.set_ylabel('Compression ratio')

    fig.tight_layout()
    plt.savefig(output_file)


def gen_plot_err_size(data, output_file):
    width = 0.8
    curves = 'dyn'

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    ax.bar(np.arange(len(start_bits)) - 2 * width / 5, data['linear']              [curves]['err_size'], width=width / 5, label='linear')
    ax.bar(np.arange(len(start_bits)) - 1 * width / 5, data['unbounded']           [curves]['err_size'], width=width / 5, label='unbounded')
    ax.bar(np.arange(len(start_bits)) + 0 * width / 5, data['unbounded_to_bounded'][curves]['err_size'], width=width / 5, label='unbounded_to_bounded')
    ax.bar(np.arange(len(start_bits)) + 1 * width / 5, data['upperbound']          [curves]['err_size'], width=width / 5, label='upperbound')
    ax.bar(np.arange(len(start_bits)) + 2 * width / 5, data['twobounds']           [curves]['err_size'], width=width / 5, label='twobounds')

    ax.legend()

    ax.set_xticks(np.arange(len(start_bits)), start_bits)
    ax.set_xlabel('Start number of bits')
    ax.set_ylabel('RMSE / size')

    fig.tight_layout()
    plt.savefig(output_file)



###############################################################################

if __name__ == '__main__':
    main()

###############################################################################

# curves = 'dyn'

# width = 0.5

# fig, ax = plt.subplots(1, 1)

# ax.scatter(stats[cave_list[0]]['linear']              [curves]['size'], stats[cave_list[0]]['linear']              [curves]['error'], label='linear')
# ax.scatter(stats[cave_list[0]]['unbounded']           [curves]['size'], stats[cave_list[0]]['unbounded']           [curves]['error'], label='unbounded')
# ax.scatter(stats[cave_list[0]]['unbounded_to_bounded'][curves]['size'], stats[cave_list[0]]['unbounded_to_bounded'][curves]['error'], label='unbounded_to_bounded')
# ax.scatter(stats[cave_list[0]]['upperbound']          [curves]['size'], stats[cave_list[0]]['upperbound']          [curves]['error'], label='upperbound')
# ax.scatter(stats[cave_list[0]]['twobounds']           [curves]['size'], stats[cave_list[0]]['twobounds']           [curves]['error'], label='twobounds')

# ax.legend()

# ax.set_xlabel('Size (MB)')
# ax.set_ylabel('RMSE')

# plt.show()
