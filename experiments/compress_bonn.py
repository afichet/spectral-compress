#!/usr/bin/env python3

import os, subprocess

from common import get_path_bonn_in, get_path_bonn_out

path_data = '/home/afichet/spectral_images/EXRs/Bonn/'
path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/compress'
path_out  = 'bonn'

techniques = ['linear', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits = [8]
flat_compression = [True, False]
flat_quantization = [True, False]

c_dc = 0
c_ac = 1

def run_compressor(
    input_file, output_file,
    log_file, binlog_file, dump_file,
    technique,
    n_bits_start, flat_quantization,
    compression_ac, compression_start_dc, flat_compression):

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    fp, fn = os.path.split(log_file)
    os.makedirs(fp, exist_ok=True)

    args = [path_bin,
        spectral_image, output_file,
        '-l', log_file,
        '-k', binlog_file,
        '-d', dump_file,
        '-m', tech,
        '-q', str(n_bits_start),
        '-a', str(compression_ac),
        '-b', str(compression_start_dc),
        '-e', str(7)
        ]

    if flat_quantization:
        args.append('--q_flat')

    if flat_compression:
        args.append('--c_flat')

    subprocess.run(args)


for d in os.listdir(path_data):
    for type in ['diffuse', 'specular']:
        spectral_image = get_path_bonn_in(path_data, d, type)

        for tech in techniques:
            for bits in start_bits:
                for q_flat in flat_quantization:
                    for c_flat in flat_compression:
                        path_final = get_path_bonn_out(
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

                        print(spectral_image, output_file)

                        run_compressor(
                            spectral_image, output_file,
                            log_file, binlog_file, dump_file,
                            tech,
                            bits, q_flat,
                            c_dc, c_ac, c_flat)
