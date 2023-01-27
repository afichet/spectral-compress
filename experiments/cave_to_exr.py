#!/usr/bin/env python3

import os
import subprocess

cave_png_dirs = '/home/afichet/spectral_images/CAVE'
output_exr_dir = '/home/afichet/spectral_images/EXRs/CAVE'
bin_dir = os.path.join('..', 'build', 'bin')


def main():
    os.makedirs(output_exr_dir, exist_ok=True)

    for f in os.listdir(cave_png_dirs):
        print(f)
        path_first_png = os.path.join(cave_png_dirs, f, f, f + '_01.png')
        path_output    = os.path.join(output_exr_dir, f + '.exr')

        subprocess.run([
            os.path.join(bin_dir, 'cave-exr'),
            path_first_png,
            path_output,
            '-t', 'float'
            ]
        )


if __name__ == '__main__':
    main()
