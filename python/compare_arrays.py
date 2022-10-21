#!/usr/bin/env python3

import struct
import numpy as np


def main():
    with open('test.dat', 'rb') as f_c:
        data = f_c.read()
        sz, = struct.unpack_from('I', data)
        print(sz)

        c_data = struct.unpack_from(
            str(sz) + 'd',
            data,
            offset=struct.calcsize('I'))

    print(c_data)


if __name__ == '__main__':
    main()