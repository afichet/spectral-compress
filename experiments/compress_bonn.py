#!/usr/bin/env python3

import os, subprocess

path_data = '/home/afichet/spectral_images/EXRs/Bonn/'
path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/compress'
path_out  = 'bonn'

if not os.path.exists(path_out):
    os.mkdir(path_out)

techniques = ['linear', 'unbounded', 'unbounded_to_bounded', 'upperbound', 'twobounds']
start_bits = [6, 8, 10, 12]


for d in os.listdir(path_data):
    path_set = os.path.join(path_out, d[:-4])

    if not os.path.exists(path_set):
        os.mkdir(path_set)

    for s in ['diffuse', 'specular']:
        spectral_image = os.path.join(path_data, d, s + '.exr')

        path_img = os.path.join(path_set, s)

        if not os.path.exists(path_img):
            os.mkdir(path_img)

        for tech in techniques:
            path_tech = os.path.join(path_img, tech)

            if not os.path.exists(path_tech):
                os.mkdir(path_tech)

            for bits in start_bits:
                path_bits = os.path.join(path_tech, str(bits))

                if not os.path.exists(path_bits):
                    os.mkdir(path_bits)

                output_file = os.path.join(path_bits, s) + '.jxl'
                print(spectral_image, output_file)
                subprocess.run([path_bin, spectral_image, output_file, '-m', tech])
                