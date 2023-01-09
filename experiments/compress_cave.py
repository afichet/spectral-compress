#!/usr/bin/env python3

import os

import common

techniques            = ['linear', 'linavg', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits            = [8]
flat_quantization     = [True, False]
flat_compression      = [True, False]
c_dc                  = 0
c_ac                  = 1
downsampling_ratio_ac = 1

path_data = '/home/afichet/spectral_images/EXRs/CAVE/'
path_out  = 'cave'

for filename in os.listdir(path_data):
    dataset_name = filename[:-4]
    spectral_image = common.get_path_cave_in(path_data, dataset_name)

    print(dataset_name)

    for tech in techniques:
        print('  ', tech)

        for bits in start_bits:
            for q_flat in flat_quantization:
                for c_flat in flat_compression:
                    path_curr_out = common.get_path_cave_out(
                        path_out,
                        dataset_name,
                        tech,
                        bits,
                        c_dc, c_ac,
                        q_flat, c_flat)

                    output_file = os.path.join(path_curr_out, dataset_name + '.jxl')
                    log_file    = os.path.join(path_curr_out, dataset_name + '.txt')
                    binlog_file = os.path.join(path_curr_out, dataset_name + '.bin')
                    dump_file   = os.path.join(path_curr_out, dataset_name + '.dat')

                    common.run_compressor(
                        spectral_image, output_file,
                        log_file, binlog_file, dump_file,
                        tech,
                        bits, q_flat,
                        c_dc, c_ac, c_flat,
                        downsampling_ratio_ac)
