#!/usr/bin/env python3

import os

from common import get_path_cave, run_compressor

path_data = '/home/afichet/spectral_images/EXRs/CAVE/'
path_out  = 'cave'

techniques = ['linear', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits = [6, 8, 10, 12]
# start_compress = [0.1]
flat_quantization = [True, False]
flat_compression  = [True, False]

flat_curves = [True, False]


for filename in os.listdir(path_data):
    spectral_image = os.path.join(path_data, filename)
    dataset_name = filename[:-4]

    for tech in techniques:
        for bits in start_bits:
            for flat in flat_curves:
                path_curr_out = get_path_cave(path_out, dataset_name, tech, bits, flat)

                output_file = os.path.join(path_curr_out, dataset_name + '.jxl')
                log_file    = os.path.join(path_curr_out, dataset_name + '.txt')
                binlog_file = os.path.join(path_curr_out, dataset_name + '.bin')

                print(spectral_image, output_file)

                run_compressor(
                    spectral_image, output_file,
                    log_file, binlog_file,
                    tech, bits, flat, 0.1, flat)
