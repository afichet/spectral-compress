#!/usr/bin/env python3

import os, subprocess

import common

techniques            = ['linear', 'linavg', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits            = [8]
flat_quantization     = [True, False]
flat_compression      = [True, False]
c_dc                  = 0
c_ac                  = 1
downsampling_ratio_ac = 1

path_data = '/home/afichet/spectral_images/EXRs/Bonn/'
path_out  = 'bonn_2'

for d in os.listdir(path_data):
    print(d)

    for type in ['diffuse', 'specular']:
        spectral_image = common.get_path_bonn_in(path_data, d, type)

        print(' ', type)

        for tech in techniques:
            print('  ', tech)

            for bits in start_bits:
                for q_flat in flat_quantization:
                    for c_flat in flat_compression:
                        path_final = common.get_path_bonn_out(
                            path_out,
                            d,
                            type,
                            tech,
                            bits,
                            c_dc, c_ac,
                            q_flat, c_flat)

                        output_file = os.path.join(path_final, d[:-4] + '.jxl')
                        log_file    = os.path.join(path_final, d[:-4] + '.txt')
                        binlog_file = os.path.join(path_final, d[:-4] + '.bin')
                        dump_file   = os.path.join(path_final, d[:-4] + '.dat')

                        common.run_compressor(
                            spectral_image, output_file,
                            log_file, binlog_file, dump_file,
                            tech,
                            bits, q_flat,
                            c_dc, c_ac, c_flat,
                            downsampling_ratio_ac)
