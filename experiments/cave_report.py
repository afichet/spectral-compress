#!/usr/bin/env python3

import os, math
import common
import openexr.spectralexr as sexr

techniques             = ['simple', 'linavg']
start_bits             = [16]
flat_quantization      = [True]
flat_compression       = ['c_flat', 'c_dynamic', 'c_deterministic']
frame_distances        = [(0, 1), (0.5, 2)]
frame_distances_simple = [(0.1, 0), (0.5, 0), (1, 0), (1.5, 0), (2, 0), (2.5, 0)]
subsampling_ratios_ac  = [1, 2]

crop_size      = 50
path_data      = '/home/afichet/spectral_images/EXRs/CAVE/'
path_out       = 'cave'
subpath_report = 'cave'
path_report    = os.path.join('export', subpath_report)
db             = [ d[:-4] for d in os.listdir(path_data) ]
db.sort()

exposure_cave = -6.5

def get_tex_stream(path_report, dataset):
    with open(os.path.join('export', 'mat_template.tex'), 'r') as f:
        stream = ''.join(f.readlines())

        stream = stream.replace('\\exportdir', path_report)
        stream = stream.replace('\\database', 'cave')
        stream = stream.replace('\\matvariant', dataset)
        stream = stream.replace('\\material', dataset)

        return stream


def get_avg_stats(
    stats: dict,
    dataset: list,
    techniques: list,
    start_bits: list,
    subsampling_ratios_ac: list,
    frame_distances: list,
    flat_quantization: list,
    flat_compression: list):
    # Initialize all fields
    avg_stats = {}
    for subsampling in subsampling_ratios_ac:
        avg_stats[subsampling] = {}

        for tech in techniques:
            avg_stats[subsampling][tech] = {}

            if tech != 'simple':
                curr_start_bits        = start_bits
                curr_frame_distances   = frame_distances
                curr_flat_quantization = flat_quantization
                curr_flat_compression  = flat_compression
            else:
                curr_start_bits        = [32]
                curr_frame_distances   = frame_distances_simple
                curr_flat_quantization = [True]
                curr_flat_compression  = ['c_flat']
                if subsampling != 1:
                    continue

            for bits in curr_start_bits:
                avg_stats[subsampling][tech][bits] = {}

                for c_dc, c_ac in curr_frame_distances:
                    if not c_dc in avg_stats[subsampling][tech][bits]:
                        avg_stats[subsampling][tech][bits][c_dc] = {}

                    avg_stats[subsampling][tech][bits][c_dc][c_ac] = {}

                    for q_flat in curr_flat_quantization:
                        avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat] = {}

                        for c_flat in curr_flat_compression:
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat] = {}

                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse']   = 0
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio' ] = 0
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_c'] = 0
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_d'] = 0

                            for k in ['q_curve', 'c_curve', 'min_curve', 'max_curve']:
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k] = []
                                for i in range(len(stats[dataset[0]][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k])):
                                    avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k].append(0)

    div = 1 / len(dataset)

    for d in dataset:
        for subsampling in subsampling_ratios_ac:
            for tech in techniques:
                if tech != 'simple':
                    curr_start_bits        = start_bits
                    curr_frame_distances   = frame_distances
                    curr_flat_quantization = flat_quantization
                    curr_flat_compression  = flat_compression
                else:
                    curr_start_bits        = [32]
                    curr_frame_distances   = frame_distances_simple
                    curr_flat_quantization = [True]
                    curr_flat_compression  = ['c_flat']
                    if subsampling != 1:
                        continue

                for bits in curr_start_bits:
                    for c_dc, c_ac in curr_frame_distances:
                        for q_flat in curr_flat_quantization:
                            for c_flat in curr_flat_compression:
                                width    = stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['width']
                                height   = stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['height']
                                n_pixels = width * height

                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse']     += div * stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse']
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio']    += div * stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio']
                                # I want miliseconds for this one
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_c'] += 1000 * div * stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_c'] / n_pixels
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_d'] += 1000 * div * stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['time_d'] / n_pixels

                                for k in ['q_curve', 'c_curve', 'min_curve', 'max_curve']:
                                    for i in range(len(stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k])):
                                        avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k][i] += div * stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat][k][i]

    return avg_stats


def run_for(
    stats: dict,
    spectral_image:str, width: int, height: int,
    path_out:str, path_report:str,
    subsampling: int,
    dataset:str,
    technique: str,
    bits: int,
    c_dc: float, c_ac: float,
    q_flat: bool, c_type: str,
    curr_max_err: float):
    path_curr_in = common.get_path_cave_out(
        path_out, subsampling,
        dataset,
        technique,
        bits,
        c_dc, c_ac,
        q_flat, c_type)

    path_curr_out = common.get_path_cave_out(
        path_report, subsampling,
        dataset,
        technique,
        bits,
        c_dc, c_ac,
        q_flat, c_type)

    # Inputs
    compressed_file = os.path.join(path_curr_in, dataset + '.jxl')
    log_file        = os.path.join(path_curr_in, dataset + '.txt')
    binlog_file     = os.path.join(path_curr_in, dataset + '.bin')
    dump_file       = os.path.join(path_curr_in, dataset + '.dat')

    # Outputs
    decompressed_exr_file  = os.path.join(path_curr_out, dataset + '.exr')
    decompressed_png_file  = os.path.join(path_curr_out, dataset + '.png')
    timing_decompress_file = os.path.join(path_curr_out, dataset + '_decompress_time.txt')
    min_max_file           = os.path.join(path_curr_out, dataset + '_min_max.txt')
    diff_png_file          = os.path.join(path_curr_out, dataset + '_diff.png')
    diff_error_file        = os.path.join(path_curr_out, dataset + '_diff.bin')
    meta_file_size_file    = os.path.join(path_curr_out, dataset + '_size.txt')

    cropped_decompressed_png_file = os.path.join(path_curr_out, dataset + '_cropped.png')
    cropped_diff_png_file         = os.path.join(path_curr_out, dataset + '_diff_cropped.png')

    common.run_decompressor(compressed_file, decompressed_exr_file, technique, timing_decompress_file)
    common.run_sgeg_min_max(compressed_file, min_max_file)
    common.run_converter_exr_png(decompressed_exr_file, decompressed_png_file, exposure_cave)
    common.run_diff(spectral_image, decompressed_exr_file, curr_max_err, diff_png_file, diff_error_file)

    common.crop_png(decompressed_png_file, cropped_decompressed_png_file, crop_size)
    common.crop_png(diff_png_file, cropped_diff_png_file, crop_size)

    org_file_size = os.path.getsize(spectral_image)

    size     = common.get_jxl_dir_size(path_curr_in)
    ratio    = org_file_size / size
    rmse     = common.get_rmse_from_diff_bin(diff_error_file)
    time_c   = common.get_duration_from_txt_log(log_file, technique)
    time_d   = common.get_duration_decompression(timing_decompress_file)

    if technique != 'simple':
        q_curve  = common.get_q_curve_from_txt_log(log_file)
        c_curve  = common.get_c_curve_from_txt_log(log_file)
        min_curve, max_curve = common.get_min_max_from(min_max_file)
    else:
        q_curve = []
        c_curve = []
        min_curve, max_curve = [], []

    # Dynamically compute max bound
    curr_max_err = max(curr_max_err, common.get_five_percentile_from_diff_bin(diff_error_file))

    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['width']     = width
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['height']    = height
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['size']      = size
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['ratio']     = ratio
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['rmse']      = rmse
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['q_curve']   = q_curve
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['c_curve']   = c_curve
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['time_c']    = time_c
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['time_d']    = time_d
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['min_curve'] = min_curve
    stats[dataset][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['max_curve'] = max_curve

    with open(meta_file_size_file, 'w') as f:
        f.write('{:.2f} MiB'.format(size / (1000 * 1000)))

    return curr_max_err


def main():
    stats = {}
    for d in db:
        stats[d] = {}

        for subsampling in subsampling_ratios_ac:
            stats[d][subsampling] = {}

            for tech in techniques:
                stats[d][subsampling][tech] = {}

                if tech != 'simple':
                    curr_start_bits        = start_bits
                    curr_frame_distances   = frame_distances
                    curr_flat_quantization = flat_quantization
                    curr_flat_compression  = flat_compression
                else:
                    curr_start_bits        = [32]
                    curr_frame_distances   = frame_distances_simple
                    curr_flat_quantization = [True]
                    curr_flat_compression  = ['c_flat']

                for bits in curr_start_bits:
                    stats[d][subsampling][tech][bits] = {}

                    for c_dc, c_ac in curr_frame_distances:
                        if not c_dc in stats[d][subsampling][tech][bits]:
                            stats[d][subsampling][tech][bits][c_dc] = {}

                        stats[d][subsampling][tech][bits][c_dc][c_ac] = {}

                        for q_flat in curr_flat_quantization:
                            stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat] = {}

                            for c_flat in curr_flat_compression:
                                stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat] = {}

    tex_stream = ''
    x_max, y_max = 0, 0

    for d in db:
        print('dataset:', d)

        tex_stream += '\n\\subsection{' + d.replace('_', ' ') + '}\n'

        org_exr_file = common.get_path_cave_in(path_data, d)
        org_png_file = os.path.join(path_report, d, d + '.png')
        org_file_size = os.path.getsize(org_exr_file)

        cropped_org_png_file = os.path.join(path_report, d, d +'_cropped.png')

        # Retrive misc. info about the original OpenEXR file
        org_exr_image = sexr.SpectralEXRFile(org_exr_file)

        spectrum_type = ''
        wl_min, wl_max = 0, 0
        n_bands = 0
        width, height = org_exr_image.width, org_exr_image.height

        # Two ifs not else to account for cases where the file contains
        # both emissive and reflective channels
        if org_exr_image.is_emissive:
            spectrum_type = 'Emissive '
            n_bands = len(org_exr_image.emissive_wavelengths_nm[0])
            wl_min = org_exr_image.emissive_wavelengths_nm[0][0]
            wl_max = org_exr_image.emissive_wavelengths_nm[0][-1]
        if org_exr_image.is_reflective:
            spectrum_type += 'Reflective'
            n_bands += len(org_exr_image.reflective_wavelengths_nm)
            wl_min = org_exr_image.reflective_wavelengths_nm[0]
            wl_max = org_exr_image.reflective_wavelengths_nm[-1]

        common.run_converter_exr_png(org_exr_file, org_png_file, exposure_cave)
        common.crop_png(org_png_file, cropped_org_png_file, crop_size)

        for subsampling in subsampling_ratios_ac:
            print('  subsampling:', subsampling)
            for tech in techniques:
                curr_max_err = 0

                if tech != 'simple':
                    curr_start_bits        = start_bits
                    curr_frame_distances   = frame_distances
                    curr_flat_quantization = flat_quantization
                    curr_flat_compression  = flat_compression
                else:
                    curr_start_bits        = [32]
                    curr_frame_distances   = frame_distances_simple
                    curr_flat_quantization = [True]
                    curr_flat_compression  = ['c_flat']
                    if subsampling != 1:
                        continue

                print('    technique:', tech)

                for bits in curr_start_bits:
                    for c_dc, c_ac in curr_frame_distances:
                        for q_flat in curr_flat_quantization:
                            for c_flat in curr_flat_compression:
                                m = run_for(
                                    stats,
                                    org_exr_file, width, height,
                                    path_out, path_report,
                                    subsampling,
                                    d,
                                    tech,
                                    bits,
                                    c_dc, c_ac,
                                    q_flat, c_flat,
                                    curr_max_err)

                                x_max = max(x_max, stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio'])
                                y_max = max(y_max, stats[d][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse'])

                                if not math.isnan(m) and not math.isinf(m) and m < 1:
                                    curr_max_err = m
                                else:
                                    print("      Err max err - Ignoring this max value for heatmap scale")


                # Run diff a second time with the dynamically computed max bound value
                # path_curr_out_partial = common.get_path_cave_out_partial(path_report, subsampling, d, tech)
                # meta_file_max_diff_file = os.path.join(path_curr_out_partial, d + '_max_err.txt')

                # with open(meta_file_max_diff_file, 'w') as f:
                #     a, b = '{:.1E}'.format(curr_max_err).split('E')
                #     string_err = '{}\\cdot 10^{{{}}}'.format(float(a), int(b))
                #     # string_err = string_err.replace('1.0\cdot', '')
                #     f.write(string_err)

                # for bits in curr_start_bits:
                #     for c_dc, c_ac in curr_frame_distances:
                #         for q_flat in curr_flat_quantization:
                #             for c_flat in curr_flat_compression:
                #                 path_curr_out = common.get_path_cave_out(
                #                     path_report, subsampling,
                #                     d,
                #                     tech,
                #                     bits,
                #                     c_dc, c_ac,
                #                     q_flat, c_flat)

                #                 decompressed_exr_file   = os.path.join(path_curr_out, d + '.exr')
                #                 diff_png_file           = os.path.join(path_curr_out, d + '_diff.png')
                #                 diff_error_file         = os.path.join(path_curr_out, d + '_diff.bin')
                #                 cropped_diff_png_file   = os.path.join(path_curr_out, d + '_diff_cropped.png')

                #                 common.run_diff(org_exr_file, decompressed_exr_file, curr_max_err, diff_png_file, diff_error_file)
                #                 common.crop_png(diff_png_file, cropped_diff_png_file, crop_size)


        plot_rmse_file           = os.path.join(path_report, d, d + '_rmse.pgf')
        plot_size_file           = os.path.join(path_report, d, d + '_size.pgf')
        plot_ratio_file          = os.path.join(path_report, d, d + '_ratio.pgf')
        plot_time_c_file         = os.path.join(path_report, d, d + '_time_c.pgf')
        plot_time_d_file         = os.path.join(path_report, d, d + '_time_d.pgf')
        plot_min_max_file        = os.path.join(path_report, d, d + '_min_max.pgf')
        plot_xy_ratio_error_file = os.path.join(path_report, d, d + '_xy_ratio_err.pgf')
        plot_c_curve_file        = os.path.join(path_report, d, d + '_c_curve.pgf')
        plot_legend              = os.path.join(path_report, d, d + '_legend.pgf')

        meta_org_file_size_file  = os.path.join(path_report, d, d + '_org_size.txt')
        meta_org_file_dim_file   = os.path.join(path_report, d, d + '_dim.txt')
        meta_wl_range_file       = os.path.join(path_report, d, d + '_wl_range.txt')
        meta_n_bands_file        = os.path.join(path_report, d, d + '_n_bands.txt')
        meta_spectrum_type_file  = os.path.join(path_report, d, d + '_spectrum_type.txt')

        common.plot_rmse(                  plot_rmse_file   , stats[d], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .25, 15, 12)
        common.plot_compression_ratio(     plot_ratio_file  , stats[d], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .25, 15, 12)
        common.plot_duration_compression(  plot_time_c_file , stats[d], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .25, 15, 12)
        common.plot_duration_decompression(plot_time_d_file , stats[d], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .25, 15, 12)
        common.plot_min_max(               plot_min_max_file, stats[d], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5, 15, 8)

        common.plot_xy_ratio_error_curves(plot_xy_ratio_error_file, stats, [d], [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 80, .65, 15, 12)
        common.plot_c_curves(             plot_c_curve_file       , stats[d]                              , 'linavg', 16, subsampling_ratios_ac, frame_distances, .5, 15, 8)

        common.plot_legend_1(plot_legend, frame_distances_simple, subsampling_ratios_ac, frame_distances)

        with open(meta_org_file_size_file, 'w') as f:
            f.write('{:.2f} MiB'.format(org_file_size / (1000 * 1000)))

        with open(meta_spectrum_type_file, 'w') as f:
            f.write(spectrum_type)

        with open(meta_wl_range_file, 'w') as f:
            f.write('[{}, {}] nm'.format(wl_min, wl_max))

        with open(meta_n_bands_file, 'w') as f:
            f.write('{}'.format(n_bands))

        with open(meta_org_file_dim_file, 'w') as f:
            f.write('{}x{} px'.format(width, height))

        # Populate LaTeX stream
        tex_stream += get_tex_stream(subpath_report, d)
        tex_stream += '\n\\clearpage\n'

    plot_avg_rmse_file                = os.path.join(path_report, 'avg_rmse.pgf')
    plot_avg_ratio_file               = os.path.join(path_report, 'avg_ratio.pgf')
    plot_avg_time_c_file              = os.path.join(path_report, 'avg_time_c.pgf')
    plot_avg_time_d_file              = os.path.join(path_report, 'avg_time_d.pgf')
    plot_avg_min_max_file             = os.path.join(path_report, 'avg_min_max.pgf')
    # plot_avg_xy_ratio_error_file      = os.path.join(path_report, 'avg_xy_ratio_err.pgf')
    plot_avg_xy_ratio_error_flat_file = os.path.join(path_report, 'avg_xy_ratio_err_flat.pgf')
    plot_avg_xy_ratio_error_det_file  = os.path.join(path_report, 'avg_xy_ratio_err_det.pgf')
    plot_avg_xy_ratio_error_dyn_file  = os.path.join(path_report, 'avg_xy_ratio_err_dyn.pgf')
    plot_avg_c_curve_file             = os.path.join(path_report, 'avg_c_curve.pgf')
    plot_avg_legend                   = os.path.join(path_report, 'avg_legend.pgf')

    # For the paper, the figures does not have the same size, we regenerate those with the proper scaling
    plot_paper_rmse_file                = os.path.join(path_report, 'rmse.pgf')
    plot_paper_ratio_file               = os.path.join(path_report, 'ratio.pgf')
    plot_paper_time_c_file              = os.path.join(path_report, 'time_c.pgf')
    plot_paper_time_d_file              = os.path.join(path_report, 'time_d.pgf')
    plot_paper_xy_ratio_error_flat_file = os.path.join(path_report, 'xy_ratio_err_flat.pgf')
    plot_paper_xy_ratio_error_det_file  = os.path.join(path_report, 'xy_ratio_err_det.pgf')
    plot_paper_xy_ratio_error_dyn_file  = os.path.join(path_report, 'xy_ratio_err_dyn.pgf')
    plot_paper_c_curve_file             = os.path.join(path_report, 'c_curve.pgf')
    plot_paper_legend                   = os.path.join(path_report, 'legend.pgf')

    avg_stats = get_avg_stats(stats, db, techniques, start_bits, subsampling_ratios_ac, frame_distances, flat_quantization, flat_compression)

    x_max += 5/100 * x_max
    y_max += 5/100 * y_max

    common.plot_rmse(                            plot_avg_rmse_file   , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .33, 15, 12)
    common.plot_compression_ratio(               plot_avg_ratio_file  , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .33, 15, 12)
    common.plot_duration_compression_per_pixel(  plot_avg_time_c_file , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .33, 15, 12)
    common.plot_duration_decompression_per_pixel(plot_avg_time_d_file , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .33, 15, 12)
    common.plot_min_max(                         plot_avg_min_max_file, avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5, 15, 10)

    # common.plot_xy_ratio_error_curves(    plot_avg_xy_ratio_error_file       , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 15)
    common.plot_xy_ratio_error_curves_key(plot_avg_xy_ratio_error_flat_file  , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 15, 'c_flat', 'Flat curves', x_max, y_max, .33, 15, 12)
    common.plot_xy_ratio_error_curves_key(plot_avg_xy_ratio_error_det_file   , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 15, 'c_deterministic', 'Deterministic curves', x_max, y_max, .33, 15, 12)
    common.plot_xy_ratio_error_curves_key(plot_avg_xy_ratio_error_dyn_file   , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 15, 'c_dynamic', 'Dynamic curves', x_max, y_max, .33, 15, 12)
    common.plot_c_curves(                 plot_avg_c_curve_file              , avg_stats                            , 'linavg', 16, subsampling_ratios_ac, frame_distances, .5, 15, 10)

    common.plot_legend_2(plot_avg_legend, frame_distances_simple, subsampling_ratios_ac, frame_distances)


    # Paper version
    common.plot_rmse(                            plot_paper_rmse_file   , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5 * .5, 15, 10)
    common.plot_compression_ratio(               plot_paper_ratio_file  , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5 * .5, 15, 10)
    common.plot_duration_compression_per_pixel(  plot_paper_time_c_file , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5 * .5, 15, 12)
    common.plot_duration_decompression_per_pixel(plot_paper_time_d_file , avg_stats, frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression, .5 * .5, 15, 12)

    common.plot_xy_ratio_error_curves_key(plot_paper_xy_ratio_error_flat_file  , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 2, 'c_flat', 'Flat curves', x_max, y_max, .33 * .5, 15, 13)
    common.plot_xy_ratio_error_curves_key(plot_paper_xy_ratio_error_det_file   , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 2, 'c_deterministic', 'Deterministic curves', x_max, y_max, .33 * .5, 15, 13)
    common.plot_xy_ratio_error_curves_key(plot_paper_xy_ratio_error_dyn_file   , stats, db, [], frame_distances_simple, 'linavg', 16, subsampling_ratios_ac, frame_distances, 2, 'c_dynamic', 'Dynamic curves', x_max, y_max, .33 * .5, 15, 13)
    common.plot_c_curves(                 plot_paper_c_curve_file              , avg_stats                            , 'linavg', 16, subsampling_ratios_ac, frame_distances, .5 * .5, 15, 12)

    common.plot_legend_2(plot_paper_legend, frame_distances_simple, subsampling_ratios_ac, frame_distances)

    with open(os.path.join(path_report + '.tex'), 'w') as f:
        f.write(tex_stream)


if __name__ == '__main__':
    main()
