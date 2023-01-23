#!/usr/bin/env python3

import os, math
import common
import openexr.spectralexr as sexr

techniques             = ['simple', 'linavg'] #['linear', 'linavg', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits             = [16]
flat_quantization      = [True] #[True, False]
flat_compression       = ['c_flat', 'c_dynamic', 'c_deterministic']
frame_distances        = [(0, 1), (0.5, 2)]
frame_distances_simple = [(0.1, 0), (1, 0), (1.5, 0), (2, 0), (3, 0), (4, 0)]
subsampling_ratios_ac  = [1, 2]

path_data      = '/home/afichet/spectral_images/EXRs/Bonn/'
path_bin       = '/home/afichet/Repositories/spectral-compress/build/bin/compress'
path_out       = 'bonn'
subpath_report = 'bonn'
path_report    = os.path.join('export', subpath_report)
variants       = ['diffuse', 'specular']
db             = [ d for d in os.listdir(path_data) ]
db.sort()

exposure_bonn = 0


def get_tex_stream(subpath_report: str, dataset: str, variant: str):
    with open(os.path.join('export', 'mat_template.tex'), 'r') as f:
        stream = ''.join(f.readlines())

        stream = stream.replace('\\exportdir', subpath_report)
        stream = stream.replace('\\database', 'bonn')
        stream = stream.replace('\\matvariant', os.path.join(dataset, variant))
        stream = stream.replace('\\material', dataset)

        return stream


def get_avg_stats(
    stats: dict,
    dataset: list,
    variants: list,
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

                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse'] = 0
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio'] = 0
                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['duration'] = 0

                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'] = []
                            for i in range(len(stats[dataset[0]][variants[0]][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'])):
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'].append(0)

                            avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'] = []
                            for i in range(len(stats[dataset[0]][variants[0]][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'])):
                                avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'].append(0)

    # TODO: we have "- 1" because one material fails the compression
    div = 1 / (len(dataset) * len(variants) - 1)
    # div = 1 / (len(dataset) * len(variants))

    for d in dataset:
        for v in variants:
            # TODO: remove, one material fails the compression
            if v == 'diffuse' and d == 'Brokat_Sorbonne_pink_aniso_high_gloss_':
                continue
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
                                    width    = stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['width']
                                    height   = stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['height']
                                    n_pixels = width * height

                                    avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse']     += div * stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['rmse']
                                    avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio']    += div * stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['ratio']
                                    avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['duration'] += div * stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['duration'] / n_pixels

                                    for i in range(len(stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'])):
                                        avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'][i] += div * stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['q_curve'][i]

                                    for i in range(len(stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'])):
                                        avg_stats[subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'][i] += div * stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat]['c_curve'][i]

    return avg_stats


def run_for(
    stats: dict,
    spectral_image: str, width: int, height: int,
    path_out:str, path_report:str,
    subsampling: int,
    dataset:str, variant:str,
    technique: str,
    bits: int,
    c_dc: float, c_ac: float,
    q_flat: bool, c_type: str,
    curr_max_err: float):
    path_curr_in = common.get_path_bonn_out(
        path_out, subsampling,
        dataset, variant,
        technique,
        bits,
        c_dc, c_ac,
        q_flat, c_type)

    path_curr_out = common.get_path_bonn_out(
        path_report, subsampling,
        dataset, variant,
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
    decompressed_exr_file = os.path.join(path_curr_out, dataset + '.exr')
    decompressed_png_file = os.path.join(path_curr_out, dataset + '.png')
    diff_png_file         = os.path.join(path_curr_out, dataset + '_diff.png')
    diff_error_file       = os.path.join(path_curr_out, dataset + '_diff.bin')
    meta_file_size_file   = os.path.join(path_curr_out, dataset + '_size.txt')

    cropped_decompressed_png_file = os.path.join(path_curr_out, dataset = '_cropped.png')
    cropped_diff_png_file         = os.path.join(path_curr_out, dataset + '_diff_cropped.png')

    common.run_decompressor(compressed_file, decompressed_exr_file, technique)
    common.run_converter_exr_png(decompressed_exr_file, decompressed_png_file, exposure_bonn)
    common.run_diff(spectral_image, decompressed_exr_file, curr_max_err, diff_png_file, diff_error_file)

    common.crop_png(decompressed_png_file, cropped_decompressed_png_file, 50)
    common.crop_png(diff_png_file, cropped_diff_png_file, 50)

    org_file_size = os.path.getsize(spectral_image)

    size     = common.get_jxl_dir_size(path_curr_in)
    ratio    = org_file_size / size
    rmse     = common.get_rmse_from_diff_bin(diff_error_file)
    duration = common.get_duration_from_txt_log(log_file, technique)

    if technique != 'simple':
        q_curve  = common.get_q_curve_from_txt_log(log_file)
        c_curve  = common.get_c_curve_from_txt_log(log_file)
    else:
        q_curve = [0]
        c_curve = [0]

    # Dynamically compute max bound
    curr_max_err = max(curr_max_err, common.get_five_percentile_from_diff_bin(diff_error_file))

    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['width']    = width
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['height']   = height
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['size']     = size
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['ratio']    = ratio
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['rmse']     = rmse
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['q_curve']  = q_curve
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['c_curve']  = c_curve
    stats[dataset][variant][subsampling][technique][bits][c_dc][c_ac][q_flat][c_type]['duration'] = duration

    with open(meta_file_size_file, 'w') as f:
        f.write('{:.2f} MiB'.format(size / (1000 * 1000)))

    return curr_max_err

def main():
    stats = {}
    for d in db:
        stats[d] = {}

        for v in variants:
            stats[d][v] = {}

            for subsampling in subsampling_ratios_ac:
                stats[d][v][subsampling] = {}

                for tech in techniques:
                    stats[d][v][subsampling][tech] = {}

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
                        stats[d][v][subsampling][tech][bits] = {}

                        for c_dc, c_ac in curr_frame_distances:
                            if not c_dc in stats[d][v][subsampling][tech][bits]:
                                stats[d][v][subsampling][tech][bits][c_dc] = {}

                            stats[d][v][subsampling][tech][bits][c_dc][c_ac] = {}

                            for q_flat in curr_flat_quantization:
                                stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat] = {}

                                for c_flat in curr_flat_compression:
                                    stats[d][v][subsampling][tech][bits][c_dc][c_ac][q_flat][c_flat] = {}

    tex_stream = ''

    for d in db:
        print('dataset:', d)

        tex_stream += '\n\\subsection{' + d.replace('_', ' ') + '}\n'

        for v in variants:
            # TODO: remove, one material fails the conversion
            if v == 'diffuse' and d == 'Brokat_Sorbonne_pink_aniso_high_gloss_':
                continue

            print('  variant:', v)

            tex_stream += '\n\\subsubsection{' + v + '}\n'

            org_exr_file = common.get_path_bonn_in(path_data, d, v)
            org_png_file = os.path.join(path_report, d, v, d + '.png')
            org_file_size = os.path.getsize(org_exr_file)
            curr_max_err = 1E-2

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

            common.run_converter_exr_png(org_exr_file, org_png_file, exposure_bonn)

            for subsampling in subsampling_ratios_ac:
                print('    subsampling:', subsampling)

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

                    print('      technique:', tech)

                    for bits in curr_start_bits:
                        for c_dc, c_ac in curr_frame_distances:
                            for q_flat in curr_flat_quantization:
                                for c_flat in curr_flat_compression:
                                    m = run_for(
                                        stats,
                                        org_exr_file, width, height,
                                        path_out, path_report,
                                        subsampling,
                                        d, v,
                                        tech,
                                        bits,
                                        c_dc, c_ac,
                                        q_flat, c_flat,
                                        curr_max_err)
                                    if not math.isnan(m) and not math.isinf(m) and m < 1:
                                        curr_max_err = m
                                    else:
                                        print("      Err max err - Ignoring this max value for heatmap scale")

                    # Run diff a second time with the dynamically computed max bound value
                    path_curr_out_partial = common.get_path_bonn_out_partial(path_report, subsampling, d, v, tech)
                    meta_file_max_diff_file = os.path.join(path_curr_out_partial, d + '_max_err.txt')

                    with open(meta_file_max_diff_file, 'w') as f:
                        a, b = '{:.1E}'.format(curr_max_err).split('E')
                        string_err = '{}\\cdot 10^{{{}}}'.format(float(a), int(b))
                        # string_err = string_err.replace('1.0\cdot', '')
                        f.write(string_err)

                    for bits in curr_start_bits:
                        for c_dc, c_ac in curr_frame_distances:
                            for q_flat in curr_flat_quantization:
                                for c_flat in curr_flat_compression:
                                    path_curr_out = common.get_path_bonn_out(
                                        path_report, subsampling,
                                        d,
                                        v,
                                        tech,
                                        bits,
                                        c_dc, c_ac,
                                        q_flat, c_flat)

                                    decompressed_exr_file   = os.path.join(path_curr_out, d + '.exr')
                                    diff_png_file           = os.path.join(path_curr_out, d + '_diff.png')
                                    diff_error_file         = os.path.join(path_curr_out, d + '_err.bin')

                                    common.run_diff(org_exr_file, decompressed_exr_file, curr_max_err, diff_png_file, diff_error_file)

            plot_curve_rmse_file     = os.path.join(path_report, d, v, d + '_rmse.pgf')
            plot_curve_size_file     = os.path.join(path_report, d, v, d + '_size.pgf')
            plot_curve_ratio_file    = os.path.join(path_report, d, v, d + '_ratio.pgf')
            plot_curve_duration_file = os.path.join(path_report, d, v, d + '_duration.pgf')
            plot_c_curve_file        = os.path.join(path_report, d, v, d + '_c_curve.pgf')
            plot_legend              = os.path.join(path_report, d, v, d + '_legend.pgf')

            meta_org_file_size_file = os.path.join(path_report, d, v, d + '_org_size.txt')
            meta_org_file_dim_file  = os.path.join(path_report, d, v, d + '_dim.txt')
            meta_wl_range_file      = os.path.join(path_report, d, v, d + '_wl_range.txt')
            meta_n_bands_file       = os.path.join(path_report, d, v, d + '_n_bands.txt')
            meta_spectrum_type_file = os.path.join(path_report, d, v, d + '_spectrum_type.txt')

            common.plot_mode_curve_rmse(plot_curve_rmse_file        , stats[d][v], 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
            common.plot_mode_curve_size (plot_curve_size_file       , stats[d][v], 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
            common.plot_mode_curve_ratio(plot_curve_ratio_file      , stats[d][v], 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
            common.plot_mode_curve_duration(plot_curve_duration_file, stats[d][v], 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
            common.plot_c_curves(plot_c_curve_file                  , stats[d][v], 'linavg', 16, subsampling_ratios_ac, frame_distances)
            common.plot_legend(plot_legend, subsampling_ratios_ac, frame_distances)

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
            tex_stream += get_tex_stream(subpath_report, d, v)
            tex_stream += '\n\\clearpage\n'

    plot_avg_curve_rmse_file    = os.path.join(path_report, 'avg_rmse.pgf')
    plot_avg_curve_ratio_file    = os.path.join(path_report, 'avg_ratio.pgf')
    plot_avg_curve_duration_file = os.path.join(path_report, 'avg_duration.pgf')
    plot_avg_c_curve_file        = os.path.join(path_report, 'avg_c_curve.pgf')
    plot_avg_legend              = os.path.join(path_report, 'avg_legend.pgf')

    avg_stats = get_avg_stats(stats, db, variants, techniques, start_bits, subsampling_ratios_ac, frame_distances, flat_quantization, flat_compression)

    common.plot_mode_curve_rmse(plot_avg_curve_rmse_file                  , avg_stats, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
    common.plot_mode_curve_ratio(plot_avg_curve_ratio_file                , avg_stats, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
    common.plot_mode_curve_duration_per_pixel(plot_avg_curve_duration_file, avg_stats, 'linavg', 16, subsampling_ratios_ac, frame_distances, flat_compression)
    common.plot_c_curves(plot_avg_c_curve_file                            , avg_stats, 'linavg', 16, subsampling_ratios_ac, frame_distances)
    common.plot_legend(plot_avg_legend, subsampling_ratios_ac, frame_distances)

    with open(os.path.join(path_report + '.tex'), 'w') as f:
        f.write(tex_stream)


if __name__ == '__main__':
    main()
