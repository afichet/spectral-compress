#!/usr/bin/env python3

import os
import common


path_data      = '/home/afichet/spectral_images/EXRs/CAVE/'
db             = [ d for d in os.listdir(path_data) ]


def main():
    for d in db:
        filename = os.path.join(path_data, d)
        common.run_strip_rgb(filename, filename)


if __name__ == '__main__':
    main()
