#!/usr/bin/env python3

import os
import numpy as np

import common
import matplotlib
import matplotlib.pyplot as plt
import openexr.spectralexr as sexr

techniques            = ["linavg"] #['linear', 'linavg', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits            = [16]
flat_quantization     = [True] #[True, False]
flat_compression      = [True, False]
c_dc                  = 0
c_ac                  = 1
downsampling_ratio_ac = 1

path_data   = '/home/afichet/spectral_images/EXRs/CAVE/'
path_bin    = '/home/afichet/Repositories/spectral-compress/build/bin/compress'
path_out    = 'cave_{}'.format(downsampling_ratio_ac)
path_report = 'cave_{}'.format(downsampling_ratio_ac)

exposure_cave = -6.5

db = [ d for d in os.listdir(path_out) ]

variants = ['diffuse', 'specular']

save_tex = True

if (save_tex):
    matplotlib.use("pgf")
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })


def get_tex_stream(path_report, dataset):
    with open(os.path.join('export', 'mat_template_cave.tex'), 'r') as f:
        stream = ''.join(f.readlines())

        stream = stream.replace('\\exportdir', path_report)
        stream = stream.replace('\\material', dataset)

        return stream


def get_avg_stats(
    stats,
    dataset,
    techniques,
    n_bits):

    y = {}
    div = 1 / len(dataset)

    for t in techniques:
        y[t] = {}
        y[t][n_bits] = {}
        for q_curve in flat_quantization:
            y[t][n_bits][q_curve] = {}
            for c_curve in flat_compression:
                y[t][n_bits][q_curve][c_curve] = {}
                y[t][n_bits][q_curve][c_curve]['error'] = 0
                y[t][n_bits][q_curve][c_curve]['ratio'] = 0
                y[t][n_bits][q_curve][c_curve]['q_curve'] = []
                y[t][n_bits][q_curve][c_curve]['c_curve'] = []

                for i in range(len(stats[dataset[0]][t][n_bits][q_curve][c_curve]['q_curve'])):
                    y[t][n_bits][q_curve][c_curve]['q_curve'].append(0)
                for i in range(len(stats[dataset[0]][t][n_bits][q_curve][c_curve]['c_curve'])):
                    y[t][n_bits][q_curve][c_curve]['c_curve'].append(0)

    for d in dataset:
        for t in techniques:
            for q_curve in flat_quantization:
                for c_curve in flat_compression:
                    y[t][n_bits][q_curve][c_curve]['error'] += div * stats[d][t][n_bits][q_curve][c_curve]['error']
                    y[t][n_bits][q_curve][c_curve]['ratio'] += div * stats[d][t][n_bits][q_curve][c_curve]['ratio']

                    for i in range(len(stats[d][t][n_bits][q_curve][c_curve]['q_curve'])):
                        y[t][n_bits][q_curve][c_curve]['q_curve'][i] += div * stats[d][t][n_bits][q_curve][c_curve]['q_curve'][i]

                    for i in range(len(stats[d][t][n_bits][q_curve][c_curve]['c_curve'])):
                        y[t][n_bits][q_curve][c_curve]['c_curve'][i] += div * stats[d][t][n_bits][q_curve][c_curve]['c_curve'][i]

    return y


def plot_mode_curves_param(
    output_filename,
    stats,
    techniques, n_bits,
    key, y_label):
    w = 0.8
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    x = np.arange(len(flat_compression) * len(flat_quantization))

    n_techinques = len(techniques)
    x_offset = n_techinques / 2 - n_techinques

    for tech, i in zip(techniques, range(n_techinques)):
        y = [
            stats[tech][n_bits][True][True][key],
            # stats[tech][n_bits][False][True][key],
            stats[tech][n_bits][True][False][key],
            # stats[tech][n_bits][False][False][key],
        ]

        x_offset = i - n_techinques / 2 + .5

        ax.bar(x + x_offset * w/n_techinques, y, width=w / n_techinques, label=tech.replace('_', ' '))

    ax.set_xticks(x, [
        'flat quant\nflat comp',
        # 'dyn quant\nflat comp',
        'flat quant\ndyn comp',
        # 'dyn quant\ndyn comp'
        ]
    )
    ax.set_xlabel('Curves')
    ax.set_ylabel(y_label)

    ax.legend()

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def plot_mode_curve_error(output_filename, stats, techniques, n_bits):
    plot_mode_curves_param(output_filename, stats, techniques, n_bits, 'error', 'Error')


def plot_mode_curve_size(output_filename, stats, techniques, n_bits):
   plot_mode_curves_param(output_filename, stats, techniques, n_bits, 'size', 'File size')


def plot_mode_curve_ratio(output_filename, stats, techniques, n_bits):
   plot_mode_curves_param(output_filename, stats, techniques, n_bits, 'ratio', 'Compression ratio')


# def plot_q_curves(output_filename, stats, techniques, n_bits):
#     fig, ax = plt.subplots(1, 1, figsize=(5, 4))

#     for tech in techniques:
#         if tech == 'upperbound' or tech == 'twobounds':
#             x = np.arange(len(stats[tech][n_bits][False][False]['q_curve'][1:-1])) + 1
#             y =               stats[tech][n_bits][False][False]['q_curve'][1:-1]
#         else:
#             x = np.arange(len(stats[tech][n_bits][False][False]['q_curve'][1:])) + 1
#             y =               stats[tech][n_bits][False][False]['q_curve'][1:]

#         ax.plot(x, y, label=tech.replace('_', ' '))

#     ax.legend()

#     ax.set_xlabel('Moment order')
#     ax.set_ylabel('Number of bits')

#     fig.tight_layout()

#     if (save_tex):
#         plt.savefig(output_filename)
#         plt.close()
#     else:
#         plt.show()


def plot_c_curves(output_filename, stats, techniques, n_bits):
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))

    for tech in techniques:
        if tech == 'upperbound' or tech == 'twobounds':
            x  = np.arange(len(stats[tech][n_bits][True][False]['c_curve'][1:-1])) + 1
            y0 =               stats[tech][n_bits][True][False]['c_curve'][1:-1]
            # y1 =               stats[tech][n_bits][False][False]['c_curve'][1:-1]
        else:
            x  = np.arange(len(stats[tech][n_bits][True][False]['c_curve'][1:])) + 1
            y0 =               stats[tech][n_bits][True][False]['c_curve'][1:]
            # y1 =               stats[tech][n_bits][False][False]['c_curve'][1:]

        ax.plot(x, y0, label=tech.replace('_', ' '))
        # ax[1].plot(x, y1, label=tech.replace('_', ' '))

    ax.legend()

    # ax[0].set_title('Without quantization')
    # ax[1].set_title('With quantization')

    ax.set_xlabel('Moment order')
    ax.set_ylabel('Framedistance (compression parameter)')

    fig.tight_layout()

    if (save_tex):
        plt.savefig(output_filename)
        plt.close()
    else:
        plt.show()


def main():
    stats = {}

    path_export = os.path.join('export', path_report)

    n_data = 0
    tex_stream = ''

    for d in db:
        stats[d] = {}
        tex_stream += '\n\\subsubsection{' + d.replace('_', ' ') + '}\n'

        org_exr_file = common.get_path_cave_in(path_data, d)
        org_png_file = os.path.join(path_export, d, d + '.png')
        org_file_size = os.path.getsize(org_exr_file)
        curr_max_err = 1E-3

        common.run_converter_exr_png(org_exr_file, org_png_file, exposure_cave)

        print(d)

        for tech in techniques:
            stats[d][tech] = {}

            print(' ', tech)

            for bits in start_bits:
                stats[d][tech][bits] = {}

                for q_flat in flat_quantization:
                    stats[d][tech][bits][q_flat] = {}

                    for c_flat in flat_compression:
                        stats[d][tech][bits][q_flat][c_flat] = {}
                        n_data += 1

                        path_curr_in = common.get_path_cave_out(
                            path_out,
                            d,
                            tech,
                            bits,
                            c_dc, c_ac,
                            q_flat, c_flat)

                        path_curr_out = common.get_path_cave_out(
                            path_export,
                            d,
                            tech,
                            bits,
                            c_dc, c_ac,
                            q_flat, c_flat)

                        # Inputs
                        compressed_file = os.path.join(path_curr_in, d + '.jxl')
                        log_file        = os.path.join(path_curr_in, d + '.txt')
                        binlog_file     = os.path.join(path_curr_in, d + '.bin')
                        dump_file       = os.path.join(path_curr_in, d + '.dat')

                        # Outputs
                        decompressed_exr_file = os.path.join(path_curr_out, d + '.exr')
                        decompressed_png_file = os.path.join(path_curr_out, d + '.png')
                        diff_png_file         = os.path.join(path_curr_out, d + '_diff.png')
                        diff_error_file       = os.path.join(path_curr_out, d + '_err.bin')

                        common.run_decompressor(compressed_file, decompressed_exr_file)
                        common.run_converter_exr_png(decompressed_exr_file, decompressed_png_file, exposure_cave)
                        common.run_diff(org_exr_file, decompressed_exr_file, curr_max_err, diff_png_file, diff_error_file)

                        size    = common.get_jxl_dir_size(path_curr_in)
                        ratio   = org_file_size / size
                        error   = common.get_err_from_diff_bin(diff_error_file)
                        q_curve = common.get_q_curve_from_txt_log(log_file)
                        c_curve = common.get_c_curve_from_txt_log(log_file)

                        stats[d][tech][bits][q_flat][c_flat]['size']    = size
                        stats[d][tech][bits][q_flat][c_flat]['ratio']   = ratio
                        stats[d][tech][bits][q_flat][c_flat]['error']   = error
                        stats[d][tech][bits][q_flat][c_flat]['q_curve'] = q_curve
                        stats[d][tech][bits][q_flat][c_flat]['c_curve'] = c_curve


        plot_curve_error_file  = os.path.join(path_export, d, d + '_error.pgf')
        plot_curve_size_file   = os.path.join(path_export, d, d + '_size.pgf')
        plot_curve_ratio_file  = os.path.join(path_export, d, d + '_ratio.pgf')
        plot_q_curve_file      = os.path.join(path_export, d, d + '_q_curve.pgf')
        plot_c_curve_file      = os.path.join(path_export, d, d + '_c_curve.pgf')

        meta_org_file_size_file = os.path.join(path_export, d, d + '_org_size.txt')
        meta_n_bands_file       = os.path.join(path_export, d, d + '_n_bands.txt')
        meta_spectrum_type_file = os.path.join(path_export, d, d + '_spectrum_type.txt')
        meta_max_error_file     = os.path.join(path_export, d, d + '_max_err.txt')

        plot_mode_curve_error(plot_curve_error_file, stats[d], techniques, 16)
        plot_mode_curve_size (plot_curve_size_file , stats[d], techniques, 16)
        plot_mode_curve_ratio(plot_curve_ratio_file, stats[d], techniques, 16)
        # plot_q_curves(plot_q_curve_file, stats[d], techniques, 8)
        plot_c_curves(plot_c_curve_file, stats[d], techniques, 16)

        with open(meta_org_file_size_file, 'w') as f:
            f.write('{:.2f} MiB'.format(org_file_size / (1000*1000)))

        # Retrive misc. info about the original OpenEXR file
        org_exr_image = sexr.SpectralEXRFile(org_exr_file)

        spectrum_type = ''
        n_bands = 0

        # Two ifs not else to account for cases where the file contains
        # both emissive and reflective channels
        if org_exr_image.is_emissive:
            spectrum_type = 'Emissive '
            n_bands = len(org_exr_image.emissive_wavelengths_nm[0])
        if org_exr_image.is_reflective:
            spectrum_type += 'Reflective'
            n_bands += len(org_exr_image.reflective_wavelengths_nm)

        with open(meta_spectrum_type_file, 'w') as f:
            f.write(spectrum_type)

        with open(meta_n_bands_file, 'w') as f:
            f.write('{}'.format(n_bands))

        with open(meta_max_error_file, 'w') as f:
            string_err = '{:.1E}'.format(curr_max_err)
            string_err = string_err.replace('E', '\cdot 10^{') + '}'
            string_err = string_err.replace('1.0\cdot', '')
            f.write(string_err)

        # Populate LaTeX stream
        tex_stream += get_tex_stream(path_report, d)
        tex_stream += '\n\\clearpage\n'

    plot_avg_curve_error_file = os.path.join(path_export, 'avg_error.pgf')
    plot_avg_curve_ratio_file = os.path.join(path_export, 'avg_ratio.pgf')
    plot_avg_q_curve_file     = os.path.join(path_export, 'avg_q_curve.pgf')
    plot_avg_c_curve_file     = os.path.join(path_export, 'avg_c_curve.pgf')

    avg_stats = get_avg_stats(stats, db, techniques, 16)
    plot_mode_curve_error(plot_avg_curve_error_file, avg_stats, techniques, 16)
    plot_mode_curve_ratio(plot_avg_curve_ratio_file, avg_stats, techniques, 16)
    # plot_q_curves(plot_avg_q_curve_file, avg_stats, techniques, 16)
    plot_c_curves(plot_avg_c_curve_file, avg_stats, techniques, 16)

    with open(os.path.join('export', path_report + '.tex'), 'w') as f:
        f.write(tex_stream)


if __name__ == '__main__':
    main()
