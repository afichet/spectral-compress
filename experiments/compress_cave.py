#!/usr/bin/env python3

import os, subprocess

path_data = '/home/afichet/spectral_images/EXRs/CAVE/'
path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/compress'
path_out  = 'cave'

if not os.path.exists(path_out):
    os.mkdir(path_out)

techniques = ['linear', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits = [8, 6, 8, 10, 12]
# start_compress = [0.1]
flat_quantization = [True, False]
flat_compression  = [True, False]

flat_curves = [True, False]

def run_compressor(
    input_file, output_file, log_file,
    technique,
    n_bits_start, flat_quantization,
    compression_start, flat_compression):

    fp, fn = os.path.split(output_file)
    os.makedirs(fp, exist_ok=True)

    fp, fn = os.path.split(log_file)
    os.makedirs(fp, exist_ok=True)

    args = [path_bin,
        spectral_image, output_file,
        '-l', log_file,
        '-m', tech,
        '-q', str(n_bits_start),
        '-a', str(0.1),
        '-b', str(compression_start),
        ]

    if flat_quantization:
        args.append('--q_flat')

    if flat_compression:
        args.append('--c_flat')

    subprocess.run(args)


for d in os.listdir(path_data):
    spectral_image = os.path.join(path_data, d)

    path_set = os.path.join(path_out, d[:-4])

    for tech in techniques:
        path_tech = os.path.join(path_set, tech)

        for bits in start_bits:
            path_bits = os.path.join(path_tech, str(bits))

            for flat in flat_curves:
                if flat:
                    path_curves = os.path.join(path_bits, 'flat')
                else:
                    path_curves = os.path.join(path_bits, 'dynamic')


                output_file = os.path.join(path_curves, d[:-4] + '.jxl')
                log_file    = os.path.join(path_curves, d[:-4] + '.txt')

                print(spectral_image, output_file)

                run_compressor(
                    spectral_image, output_file, log_file,
                    tech, bits, flat, 0.1, flat)
