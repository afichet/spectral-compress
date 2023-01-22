#!/usr/bin/env python3

import os
import common


techniques             = ['simple', 'linavg']
start_bits             = [16]
flat_quantization      = [True] #[True, False]
flat_compression       = [True, False]
frame_distances        = [(0, 1), (0.5, 2)]
frame_distances_simple = [0.1, 1, 2]
subsampling_ratio_ac   = [1, 2]

path_data = '/home/afichet/spectral_images/EXRs/CAVE/'
prefix_path_out  = 'cave'


def run_for(
    spectral_image,
    prefix_path_out: str, subsampling: int,
    dataset_name: str,
    technique: str,
    n_bits_start: int, n_exponent_bits: int,
    c_dc: float, c_ac: float,
    q_flat: bool, c_flat: bool):
    path_curr_out = common.get_path_cave_out(
        prefix_path_out, subsampling,
        dataset_name,
        technique,
        n_bits_start,
        c_dc, c_ac,
        q_flat, c_flat)

    output_file = os.path.join(path_curr_out, dataset_name + '.jxl')
    log_file    = os.path.join(path_curr_out, dataset_name + '.txt')
    binlog_file = os.path.join(path_curr_out, dataset_name + '.bin')
    dump_file   = os.path.join(path_curr_out, dataset_name + '.dat')

    common.run_compressor(
        spectral_image,
        output_file,
        log_file,
        binlog_file,
        dump_file,
        technique,
        n_bits_start, n_exponent_bits, q_flat,
        c_dc, c_ac, c_flat,
        subsampling)


def main():
    for filename in os.listdir(path_data):
        dataset_name = filename[:-4]
        spectral_image = common.get_path_cave_in(path_data, dataset_name)

        print(dataset_name)

        for tech in techniques:
            print('  ', tech)

            if tech != 'simple':
                for subsampling in subsampling_ratio_ac:
                    for bits in start_bits:
                        for c_dc, c_ac in frame_distances:
                            for q_flat in flat_quantization:
                                for c_flat in flat_compression:
                                    run_for(
                                        spectral_image,
                                        prefix_path_out, subsampling,
                                        dataset_name,
                                        tech,
                                        bits, 0,
                                        c_dc, c_ac,
                                        q_flat, c_flat)
            else:
                # placeholders
                subsampling = 1
                bits = 32
                c_ac = 0
                q_flat = True
                c_flat = True

                for framedistance in frame_distances_simple:
                    run_for(
                        spectral_image,
                        prefix_path_out, subsampling,
                        dataset_name,
                        tech,
                        bits, 0,
                        framedistance, c_ac,
                        q_flat, c_flat)


if __name__ == '__main__':
    main()
