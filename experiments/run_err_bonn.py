#!/usr/bin/env python3

import os, subprocess

path_data = '/home/afichet/spectral_images/EXRs/Bonn/'
path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/eval-quantization'
path_out  = 'output'

if not os.path.exists(path_out):
    os.mkdir(path_out)

for d in os.listdir(path_data):
    # run for Bonn material maps
    for s in ['diffuse', 'specular']:
        spectral_image = os.path.join(path_data, d, s + '.exr')
        output_file = os.path.join(path_out, d + '_' + s)

        flag_r = 'y'
        if s == 'specular':
            flag_r = 'n'

        subprocess.run([path_bin, spectral_image, output_file, flag_r])
