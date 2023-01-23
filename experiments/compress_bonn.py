#!/usr/bin/env python3

import os
import common


techniques             = ['simple', 'linavg']
start_bits             = [16]
flat_quantization      = [True] #[True, False]
flat_compression       = ['c_flat', 'c_dynamic', 'c_deterministic']
frame_distances        = [(0, 1), (0.5, 2)]
frame_distances_simple = [(0.1, 0), (1, 0), (1.5, 0), (2, 0), (3, 0), (4, 0)]
subsampling_ratio_ac   = [1, 2]

path_data = '/home/afichet/spectral_images/EXRs/Bonn/'
path_out  = 'bonn'


def run_for(
    spectral_image,
    prefix_path_out: str, subsampling: int,
    dataset_name: str, variant: str,
    technique: str,
    n_bits_start: int, n_exponent_bits: int,
    c_dc: float, c_ac: float,
    q_flat: bool, c_type: str):
    path_final = common.get_path_bonn_out(
        prefix_path_out, subsampling,
        dataset_name,
        variant,
        technique,
        n_bits_start,
        c_dc, c_ac,
        q_flat, c_type)

    output_file = os.path.join(path_final, dataset_name + '.jxl')
    log_file    = os.path.join(path_final, dataset_name + '.txt')
    binlog_file = os.path.join(path_final, dataset_name + '.bin')
    dump_file   = os.path.join(path_final, dataset_name + '.dat')

    common.run_compressor(
        spectral_image,
        output_file,
        log_file,
        binlog_file,
        dump_file,
        technique,
        n_bits_start, n_exponent_bits, q_flat,
        c_dc, c_ac, c_type,
        subsampling)


def main():
    for dataset_name in os.listdir(path_data):
        print(dataset_name)

        for variant in ['diffuse', 'specular']:
            print(' ', variant)
            spectral_image = common.get_path_bonn_in(path_data, dataset_name, variant)

            for tech in techniques:
                print('  ', tech)

                if tech != 'simple':
                    curr_start_bits         = start_bits
                    curr_frame_distances    = frame_distances
                    curr_flat_quantization  = flat_quantization
                    curr_flat_compression   = flat_compression
                    curr_subsampling_ratios = subsampling_ratio_ac
                else:
                    curr_start_bits         = [32]
                    curr_frame_distances    = frame_distances_simple
                    curr_flat_quantization  = [True]
                    curr_flat_compression   = ['c_flat']
                    curr_subsampling_ratios = [1]

                for subsampling in curr_subsampling_ratios:
                    for bits in curr_start_bits:
                        for c_dc, c_ac in curr_frame_distances:
                            for q_flat in curr_flat_quantization:
                                for c_flat in curr_flat_compression:
                                    run_for(
                                        spectral_image,
                                        path_out, subsampling,
                                        dataset_name, variant,
                                        tech,
                                        bits, 0,
                                        c_dc, c_ac,
                                        q_flat, c_flat
                                    )

if __name__ == '__main__':
    main()
