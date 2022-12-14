#!/usr/bin/env python3

import os, subprocess

path_data = '/home/afichet/spectral_images/EXRs/CAVE/'
path_bin  = '/home/afichet/Repositories/spectral-compress/build/bin/eval-quantization'
path_out  = 'CAVE'

if not os.path.exists(path_out):
    os.mkdir(path_out)

for d in os.listdir(path_data):
    spectral_image = os.path.join(path_data, d)
    output_file = os.path.join(path_out, d[:-4])
    flag_r = 'n'

    subprocess.run([path_bin, spectral_image, output_file, flag_r])
