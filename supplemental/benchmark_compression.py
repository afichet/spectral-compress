#!/usr/bin/env python3

import os
import common

path_data = '/home/afichet/spectral_images/EXRs/'
path_out = 'exr'

def format_sci(val: float):
    a, b = '{:.1E}'.format(val).split('E')
    b = int(b)

    return '\\footnotesize{$' + str(a) + '\mathrm{E}^{' + str(b) + '}$}'
    # line += '& {:.2E}'.format(err)

def main():

    # -------------------------------------------------------------------------
    # CAVE database
    # -------------------------------------------------------------------------

    f_ratio = open(os.path.join('export', 'cave_ratio.tex'), 'w')
    f_error = open(os.path.join('export', 'cave_error.tex'), 'w')

    path_cave = os.path.join(path_data, 'CAVE')
    db_cave = [ d for d in os.listdir(path_cave) ]
    db_cave.sort()

    avg_ratio_pxr24       = 0
    avg_ratio_b44         = 0
    avg_ratio_b44a        = 0
    avg_ratio_simple01    = 0
    avg_ratio_simple05    = 0
    avg_ratio_simple10    = 0
    avg_ratio_simple15    = 0
    avg_ratio_simple20    = 0
    avg_ratio_simple25    = 0
    avg_ratio_dc00ac1_11  = 0
    avg_ratio_dc05ac2_11  = 0
    avg_ratio_dc00ac1_12  = 0
    avg_ratio_dc05ac2_12  = 0

    avg_rmse_pxr24      = 0
    avg_rmse_b44        = 0
    avg_rmse_simple01   = 0
    avg_rmse_simple05   = 0
    avg_rmse_simple10   = 0
    avg_rmse_simple15   = 0
    avg_rmse_simple20   = 0
    avg_rmse_simple25   = 0
    avg_rmse_dc00ac1_11 = 0
    avg_rmse_dc05ac2_11 = 0
    avg_rmse_dc00ac1_12 = 0
    avg_rmse_dc05ac2_12 = 0

    for d, i in zip(db_cave, range(len(db_cave))):
        filename_in  = os.path.join(path_cave, d)
        element_name = os.path.basename(filename_in)[:-4]

        # File Sizes

        path_pxr24      = os.path.join(path_out, 'cave', element_name + '_float_pxr24.exr')
        path_b44        = os.path.join(path_out, 'cave', element_name + '_half_b44.exr')
        path_simple01   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 0.1, 0, True, 'c_flat')
        path_simple05   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 0.5, 0, True, 'c_flat')
        path_simple10   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 1  , 0, True, 'c_flat')
        path_simple15   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 1.5, 0, True, 'c_flat')
        path_simple20   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 2  , 0, True, 'c_flat')
        path_simple25   = common.get_path_cave_out('cave', 1, element_name, 'simple', 32, 2.5, 0, True, 'c_flat')
        path_dc00ac1_11 = common.get_path_cave_out('cave', 1, element_name, 'linavg', 16, 0  , 1, True, 'c_deterministic')
        path_dc05ac2_11 = common.get_path_cave_out('cave', 1, element_name, 'linavg', 16, 0.5, 2, True, 'c_deterministic')
        path_dc00ac1_12 = common.get_path_cave_out('cave', 2, element_name, 'linavg', 16, 0  , 1, True, 'c_deterministic')
        path_dc05ac2_12 = common.get_path_cave_out('cave', 2, element_name, 'linavg', 16, 0.5, 2, True, 'c_deterministic')

        # common.run_exr_change_compression(filename_in, path_pxr24, 'float', 'pxr24')
        # common.run_exr_change_compression(filename_in, path_b44  , 'half' , 'b44')

        file_size_org        = os.path.getsize(filename_in)
        file_size_pxr24      = os.path.getsize(path_pxr24)
        file_size_b44        = os.path.getsize(path_b44)
        file_size_simple01   = common.get_jxl_dir_size(path_simple01)
        file_size_simple05   = common.get_jxl_dir_size(path_simple05)
        file_size_simple10   = common.get_jxl_dir_size(path_simple10)
        file_size_simple15   = common.get_jxl_dir_size(path_simple15)
        file_size_simple20   = common.get_jxl_dir_size(path_simple20)
        file_size_simple25   = common.get_jxl_dir_size(path_simple25)
        file_size_dc00ac1_11 = common.get_jxl_dir_size(path_dc00ac1_11)
        file_size_dc05ac2_11 = common.get_jxl_dir_size(path_dc05ac2_11)
        file_size_dc00ac1_12 = common.get_jxl_dir_size(path_dc00ac1_12)
        file_size_dc05ac2_12 = common.get_jxl_dir_size(path_dc05ac2_12)

        ratio_pxr24      = file_size_org / file_size_pxr24
        ratio_b44        = file_size_org / file_size_b44
        ratio_simple01   = file_size_org / file_size_simple01
        ratio_simple05   = file_size_org / file_size_simple05
        ratio_simple10   = file_size_org / file_size_simple10
        ratio_simple15   = file_size_org / file_size_simple15
        ratio_simple20   = file_size_org / file_size_simple20
        ratio_simple25   = file_size_org / file_size_simple25
        ratio_dc00ac1_11 = file_size_org / file_size_dc00ac1_11
        ratio_dc05ac2_11 = file_size_org / file_size_dc05ac2_11
        ratio_dc00ac1_12 = file_size_org / file_size_dc00ac1_12
        ratio_dc05ac2_12 = file_size_org / file_size_dc05ac2_12

        avg_ratio_pxr24      += ratio_pxr24
        avg_ratio_b44        += ratio_b44
        avg_ratio_simple01   += ratio_simple01
        avg_ratio_simple05   += ratio_simple05
        avg_ratio_simple10   += ratio_simple10
        avg_ratio_simple15   += ratio_simple15
        avg_ratio_simple20   += ratio_simple20
        avg_ratio_simple25   += ratio_simple25
        avg_ratio_dc00ac1_11 += ratio_dc00ac1_11
        avg_ratio_dc05ac2_11 += ratio_dc05ac2_11
        avg_ratio_dc00ac1_12 += ratio_dc00ac1_12
        avg_ratio_dc05ac2_12 += ratio_dc05ac2_12

        values_to_write = [
            ratio_pxr24,
            ratio_b44,
            ratio_simple01,
            ratio_simple05,
            ratio_simple10,
            ratio_simple15,
            ratio_simple20,
            ratio_simple25,
            ratio_dc00ac1_11,
            ratio_dc05ac2_11,
            ratio_dc00ac1_12,
            ratio_dc05ac2_12
        ]

        line = element_name.replace('_', ' ')
        for ratio in values_to_write:
            line += '& {:.2f}'.format(ratio)
        f_ratio.write(line + '\\\\\n')

        # RMSE

        filename_rmse_pxr24     = os.path.join(path_out, 'cave', element_name + '_pxr24_rmse.bin')
        filename_rmse_b44       = os.path.join(path_out, 'cave', element_name + '_b44_rmse.bin')
        filename_png_rmse_pxr24 = os.path.join(path_out, 'cave', element_name + '_pxr24_rmse.png')
        filename_png_rmse_b44   = os.path.join(path_out, 'cave', element_name + '_b44_rmse.png')

        common.run_diff(filename_in, path_pxr24, 1, filename_png_rmse_pxr24, filename_rmse_pxr24)
        common.run_diff(filename_in, path_b44  , 1, filename_png_rmse_b44  , filename_rmse_b44)

        filename_rmse_simple01   = os.path.join('export', path_simple01, element_name + '_diff.bin')
        filename_rmse_simple05   = os.path.join('export', path_simple05, element_name + '_diff.bin')
        filename_rmse_simple10   = os.path.join('export', path_simple10, element_name + '_diff.bin')
        filename_rmse_simple15   = os.path.join('export', path_simple15, element_name + '_diff.bin')
        filename_rmse_simple20   = os.path.join('export', path_simple20, element_name + '_diff.bin')
        filename_rmse_simple25   = os.path.join('export', path_simple25, element_name + '_diff.bin')
        filename_rmse_dc00ac1_11 = os.path.join('export', path_dc00ac1_11, element_name + '_diff.bin')
        filename_rmse_dc05ac2_11 = os.path.join('export', path_dc05ac2_11, element_name + '_diff.bin')
        filename_rmse_dc00ac1_12 = os.path.join('export', path_dc00ac1_12, element_name + '_diff.bin')
        filename_rmse_dc05ac2_12 = os.path.join('export', path_dc05ac2_12, element_name + '_diff.bin')

        rmse_pxr24      = common.get_rmse_from_diff_bin(filename_rmse_pxr24)
        rmse_b44        = common.get_rmse_from_diff_bin(filename_rmse_b44)
        rmse_simple01   = common.get_rmse_from_diff_bin(filename_rmse_simple01)
        rmse_simple05   = common.get_rmse_from_diff_bin(filename_rmse_simple05)
        rmse_simple10   = common.get_rmse_from_diff_bin(filename_rmse_simple10)
        rmse_simple15   = common.get_rmse_from_diff_bin(filename_rmse_simple15)
        rmse_simple20   = common.get_rmse_from_diff_bin(filename_rmse_simple20)
        rmse_simple25   = common.get_rmse_from_diff_bin(filename_rmse_simple25)
        rmse_dc00ac1_11 = common.get_rmse_from_diff_bin(filename_rmse_dc00ac1_11)
        rmse_dc05ac2_11 = common.get_rmse_from_diff_bin(filename_rmse_dc05ac2_11)
        rmse_dc00ac1_12 = common.get_rmse_from_diff_bin(filename_rmse_dc00ac1_12)
        rmse_dc05ac2_12 = common.get_rmse_from_diff_bin(filename_rmse_dc05ac2_12)

        avg_rmse_pxr24      += rmse_pxr24
        avg_rmse_b44        += rmse_b44
        avg_rmse_simple01   += rmse_simple01
        avg_rmse_simple05   += rmse_simple05
        avg_rmse_simple10   += rmse_simple10
        avg_rmse_simple15   += rmse_simple15
        avg_rmse_simple20   += rmse_simple20
        avg_rmse_simple25   += rmse_simple25
        avg_rmse_dc00ac1_11 += rmse_dc00ac1_11
        avg_rmse_dc05ac2_11 += rmse_dc05ac2_11
        avg_rmse_dc00ac1_12 += rmse_dc00ac1_12
        avg_rmse_dc05ac2_12 += rmse_dc05ac2_12

        values_to_write = [
            rmse_pxr24,
            rmse_b44,
            rmse_simple01,
            rmse_simple05,
            rmse_simple10,
            rmse_simple15,
            rmse_simple20,
            rmse_simple25,
            rmse_dc00ac1_11,
            rmse_dc05ac2_11,
            rmse_dc00ac1_12,
            rmse_dc05ac2_12
        ]

        line = element_name.replace('_', ' ')
        for err in values_to_write:
            line += '& ' + format_sci(err)
        f_error.write(line + '\\\\\n')

    n_elements = len(db_cave)

    avg_ratio_pxr24     /= n_elements
    avg_ratio_b44       /= n_elements
    avg_ratio_b44a      /= n_elements

    avg_ratio_simple01    /= n_elements
    avg_ratio_simple05    /= n_elements
    avg_ratio_simple10    /= n_elements
    avg_ratio_simple15    /= n_elements
    avg_ratio_simple20    /= n_elements
    avg_ratio_simple25    /= n_elements
    avg_ratio_dc00ac1_11  /= n_elements
    avg_ratio_dc05ac2_11  /= n_elements
    avg_ratio_dc00ac1_12  /= n_elements
    avg_ratio_dc05ac2_12  /= n_elements

    avg_rmse_pxr24      /= n_elements
    avg_rmse_b44        /= n_elements
    avg_rmse_simple01   /= n_elements
    avg_rmse_simple05   /= n_elements
    avg_rmse_simple10   /= n_elements
    avg_rmse_simple15   /= n_elements
    avg_rmse_simple20   /= n_elements
    avg_rmse_simple25   /= n_elements
    avg_rmse_dc00ac1_11 /= n_elements
    avg_rmse_dc05ac2_11 /= n_elements
    avg_rmse_dc00ac1_12 /= n_elements
    avg_rmse_dc05ac2_12 /= n_elements

    values_to_write = [
        avg_ratio_pxr24,
        avg_ratio_b44,
        # avg_ratio_b44a,
        avg_ratio_simple01,
        avg_ratio_simple05,
        avg_ratio_simple10,
        avg_ratio_simple15,
        avg_ratio_simple20,
        avg_ratio_simple25,
        avg_ratio_dc00ac1_11,
        avg_ratio_dc05ac2_11,
        avg_ratio_dc00ac1_12,
        avg_ratio_dc05ac2_12
    ]

    f_ratio.write('\\midrule\n')
    line = 'Average'
    for ratio in values_to_write:
        line += '& {:.2f}'.format(ratio)
    f_ratio.write(line + '\\\\\n')
    f_ratio.write('\\bottomrule\n')
    f_ratio.close()

    values_to_write = [
        avg_rmse_pxr24,
        avg_rmse_b44,
        avg_rmse_simple01,
        avg_rmse_simple05,
        avg_rmse_simple10,
        avg_rmse_simple15,
        avg_rmse_simple20,
        avg_rmse_simple25,
        avg_rmse_dc00ac1_11,
        avg_rmse_dc05ac2_11,
        avg_rmse_dc00ac1_12,
        avg_rmse_dc05ac2_12
    ]

    f_error.write('\\midrule\n')
    line = 'Average'
    for err in values_to_write:
        line += '& ' + format_sci(err)
    f_error.write(line + '\\\\\n')
    f_error.write('\\bottomrule\n')
    f_error.close()

    # -------------------------------------------------------------------------
    # Bonn database
    # -------------------------------------------------------------------------

    f_ratio = open(os.path.join('export', 'bonn_ratio.tex'), 'w')
    f_error = open(os.path.join('export', 'bonn_error.tex'), 'w')

    path_bonn = os.path.join(path_data, 'Bonn')
    db_bonn = [ d for d in os.listdir(path_bonn) ]
    db_bonn.sort()

    avg_ratio_pxr24 = 0
    avg_ratio_b44   = 0

    avg_ratio_simple01    = 0
    avg_ratio_simple05    = 0
    avg_ratio_simple10    = 0
    avg_ratio_simple15    = 0
    avg_ratio_simple20    = 0
    avg_ratio_simple25    = 0
    avg_ratio_dc00ac1_11  = 0
    avg_ratio_dc05ac2_11  = 0
    avg_ratio_dc00ac1_12  = 0
    avg_ratio_dc05ac2_12  = 0

    avg_rmse_pxr24      = 0
    avg_rmse_b44        = 0
    avg_rmse_simple01   = 0
    avg_rmse_simple05   = 0
    avg_rmse_simple10   = 0
    avg_rmse_simple15   = 0
    avg_rmse_simple20   = 0
    avg_rmse_simple25   = 0
    avg_rmse_dc00ac1_11 = 0
    avg_rmse_dc05ac2_11 = 0
    avg_rmse_dc00ac1_12 = 0
    avg_rmse_dc05ac2_12 = 0

    for v in ['diffuse', 'specular']:
        for d, i in zip(db_bonn, range(len(db_bonn))):
            if v == 'diffuse' and d == 'Brokat_Sorbonne_pink_aniso_high_gloss_':
                continue
            filename_in = os.path.join(path_bonn, d, v + '.exr')
            element_name = d + v

            # File Sizes

            path_pxr24      = os.path.join(path_out, 'bonn', element_name + '_float_pxr24.exr')
            path_b44        = os.path.join(path_out, 'bonn', element_name + '_half_b44.exr')
            path_simple01   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 0.1, 0, True, 'c_flat')
            path_simple05   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 0.5, 0, True, 'c_flat')
            path_simple10   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 1  , 0, True, 'c_flat')
            path_simple15   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 1.5, 0, True, 'c_flat')
            path_simple20   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 2  , 0, True, 'c_flat')
            path_simple25   = common.get_path_bonn_out('bonn', 1, d, v, 'simple', 32, 2.5, 0, True, 'c_flat')
            path_dc00ac1_11 = common.get_path_bonn_out('bonn', 1, d, v, 'linavg', 16, 0  , 1, True, 'c_deterministic')
            path_dc05ac2_11 = common.get_path_bonn_out('bonn', 1, d, v, 'linavg', 16, 0.5, 2, True, 'c_deterministic')
            path_dc00ac1_12 = common.get_path_bonn_out('bonn', 2, d, v, 'linavg', 16, 0  , 1, True, 'c_deterministic')
            path_dc05ac2_12 = common.get_path_bonn_out('bonn', 2, d, v, 'linavg', 16, 0.5, 2, True, 'c_deterministic')

            common.run_exr_change_compression(filename_in, path_pxr24, 'float', 'pxr24')
            common.run_exr_change_compression(filename_in, path_b44  , 'half' , 'b44')

            file_size_org        = os.path.getsize(filename_in)
            file_size_pxr24      = os.path.getsize(path_pxr24)
            file_size_b44        = os.path.getsize(path_b44)
            file_size_simple_01  = common.get_jxl_dir_size(path_simple01)
            file_size_simple_05  = common.get_jxl_dir_size(path_simple05)
            file_size_simple_10  = common.get_jxl_dir_size(path_simple10)
            file_size_simple_15  = common.get_jxl_dir_size(path_simple15)
            file_size_simple_20  = common.get_jxl_dir_size(path_simple20)
            file_size_simple_25  = common.get_jxl_dir_size(path_simple25)
            file_size_dc00ac1_11 = common.get_jxl_dir_size(path_dc00ac1_11)
            file_size_dc05ac2_11 = common.get_jxl_dir_size(path_dc05ac2_11)
            file_size_dc00ac1_12 = common.get_jxl_dir_size(path_dc00ac1_12)
            file_size_dc05ac2_12 = common.get_jxl_dir_size(path_dc05ac2_12)

            ratio_pxr24      = file_size_org / file_size_pxr24
            ratio_b44        = file_size_org / file_size_b44
            ratio_simple01   = file_size_org / file_size_simple_01
            ratio_simple05   = file_size_org / file_size_simple_05
            ratio_simple10   = file_size_org / file_size_simple_10
            ratio_simple15   = file_size_org / file_size_simple_15
            ratio_simple20   = file_size_org / file_size_simple_20
            ratio_simple25   = file_size_org / file_size_simple_25
            ratio_dc00ac1_11 = file_size_org / file_size_dc00ac1_11
            ratio_dc05ac2_11 = file_size_org / file_size_dc05ac2_11
            ratio_dc00ac1_12 = file_size_org / file_size_dc00ac1_12
            ratio_dc05ac2_12 = file_size_org / file_size_dc05ac2_12

            avg_ratio_pxr24      += ratio_pxr24
            avg_ratio_b44        += ratio_b44
            avg_ratio_simple01   += ratio_simple01
            avg_ratio_simple05   += ratio_simple05
            avg_ratio_simple10   += ratio_simple10
            avg_ratio_simple15   += ratio_simple15
            avg_ratio_simple20   += ratio_simple20
            avg_ratio_simple25   += ratio_simple25
            avg_ratio_dc00ac1_11 += ratio_dc00ac1_11
            avg_ratio_dc05ac2_11 += ratio_dc05ac2_11
            avg_ratio_dc00ac1_12 += ratio_dc00ac1_12
            avg_ratio_dc05ac2_12 += ratio_dc05ac2_12

            values_to_write = [
                ratio_pxr24,
                ratio_b44,
                ratio_simple01,
                ratio_simple05,
                ratio_simple10,
                ratio_simple15,
                ratio_simple20,
                ratio_simple25,
                ratio_dc00ac1_11,
                ratio_dc05ac2_11,
                ratio_dc00ac1_12,
                ratio_dc05ac2_12
            ]

            line = element_name.replace('_', ' ')
            for ratio in values_to_write:
                line += '& {:.2f}'.format(ratio)
            f_ratio.write(line + '\\\\\n')

            # RMSE

            filename_rmse_pxr24     = os.path.join(path_out, 'bonn', element_name + '_pxr24_rmse.bin')
            filename_rmse_b44       = os.path.join(path_out, 'bonn', element_name + '_b44_rmse.bin')
            filename_png_rmse_pxr24 = os.path.join(path_out, 'bonn', element_name + '_pxr24_rmse.png')
            filename_png_rmse_b44   = os.path.join(path_out, 'bonn', element_name + '_b44_rmse.png')

            common.run_diff(filename_in, path_pxr24, 1, filename_png_rmse_pxr24, filename_rmse_pxr24)
            common.run_diff(filename_in, path_b44  , 1, filename_png_rmse_b44  , filename_rmse_b44)

            filename_rmse_simple01   = os.path.join('export', path_simple01, d + '_diff.bin')
            filename_rmse_simple05   = os.path.join('export', path_simple05, d + '_diff.bin')
            filename_rmse_simple10   = os.path.join('export', path_simple10, d + '_diff.bin')
            filename_rmse_simple15   = os.path.join('export', path_simple15, d + '_diff.bin')
            filename_rmse_simple20   = os.path.join('export', path_simple20, d + '_diff.bin')
            filename_rmse_simple25   = os.path.join('export', path_simple25, d + '_diff.bin')
            filename_rmse_dc00ac1_11 = os.path.join('export', path_dc00ac1_11, d + '_diff.bin')
            filename_rmse_dc05ac2_11 = os.path.join('export', path_dc05ac2_11, d + '_diff.bin')
            filename_rmse_dc00ac1_12 = os.path.join('export', path_dc00ac1_12, d + '_diff.bin')
            filename_rmse_dc05ac2_12 = os.path.join('export', path_dc05ac2_12, d + '_diff.bin')

            rmse_pxr24      = common.get_rmse_from_diff_bin(filename_rmse_pxr24)
            rmse_b44        = common.get_rmse_from_diff_bin(filename_rmse_b44)
            rmse_simple01   = common.get_rmse_from_diff_bin(filename_rmse_simple01)
            rmse_simple05   = common.get_rmse_from_diff_bin(filename_rmse_simple05)
            rmse_simple10   = common.get_rmse_from_diff_bin(filename_rmse_simple10)
            rmse_simple15   = common.get_rmse_from_diff_bin(filename_rmse_simple15)
            rmse_simple20   = common.get_rmse_from_diff_bin(filename_rmse_simple20)
            rmse_simple25   = common.get_rmse_from_diff_bin(filename_rmse_simple25)
            rmse_dc00ac1_11 = common.get_rmse_from_diff_bin(filename_rmse_dc00ac1_11)
            rmse_dc05ac2_11 = common.get_rmse_from_diff_bin(filename_rmse_dc05ac2_11)
            rmse_dc00ac1_12 = common.get_rmse_from_diff_bin(filename_rmse_dc00ac1_12)
            rmse_dc05ac2_12 = common.get_rmse_from_diff_bin(filename_rmse_dc05ac2_12)

            avg_rmse_pxr24      += rmse_pxr24
            avg_rmse_b44        += rmse_b44
            avg_rmse_simple01   += rmse_simple01
            avg_rmse_simple05   += rmse_simple05
            avg_rmse_simple10   += rmse_simple10
            avg_rmse_simple15   += rmse_simple15
            avg_rmse_simple20   += rmse_simple20
            avg_rmse_simple25   += rmse_simple25
            avg_rmse_dc00ac1_11 += rmse_dc00ac1_11
            avg_rmse_dc05ac2_11 += rmse_dc05ac2_11
            avg_rmse_dc00ac1_12 += rmse_dc00ac1_12
            avg_rmse_dc05ac2_12 += rmse_dc05ac2_12

            values_to_write = [
                rmse_pxr24,
                rmse_b44,
                rmse_simple01,
                rmse_simple05,
                rmse_simple10,
                rmse_simple15,
                rmse_simple20,
                rmse_simple25,
                rmse_dc00ac1_11,
                rmse_dc05ac2_11,
                rmse_dc00ac1_12,
                rmse_dc05ac2_12
            ]

            line = element_name.replace('_', ' ')
            for err in values_to_write:
                line += '& ' + format_sci(err)
            f_error.write(line + '\\\\\n')

    n_elements = len(db_bonn) * 2 - 1

    avg_ratio_pxr24     /= n_elements
    avg_ratio_b44       /= n_elements

    avg_ratio_simple01    /= n_elements
    avg_ratio_simple05    /= n_elements
    avg_ratio_simple10    /= n_elements
    avg_ratio_simple15    /= n_elements
    avg_ratio_simple20    /= n_elements
    avg_ratio_simple25    /= n_elements
    avg_ratio_dc00ac1_11  /= n_elements
    avg_ratio_dc05ac2_11  /= n_elements
    avg_ratio_dc00ac1_12  /= n_elements
    avg_ratio_dc05ac2_12  /= n_elements

    avg_rmse_pxr24      /= n_elements
    avg_rmse_b44        /= n_elements
    avg_rmse_simple01   /= n_elements
    avg_rmse_simple05   /= n_elements
    avg_rmse_simple10   /= n_elements
    avg_rmse_simple15   /= n_elements
    avg_rmse_simple20   /= n_elements
    avg_rmse_simple25   /= n_elements
    avg_rmse_dc00ac1_11 /= n_elements
    avg_rmse_dc05ac2_11 /= n_elements
    avg_rmse_dc00ac1_12 /= n_elements
    avg_rmse_dc05ac2_12 /= n_elements


    values_to_write = [
        avg_ratio_pxr24,
        avg_ratio_b44,
        avg_ratio_simple01,
        avg_ratio_simple05,
        avg_ratio_simple10,
        avg_ratio_simple15,
        avg_ratio_simple20,
        avg_ratio_simple25,
        avg_ratio_dc00ac1_11,
        avg_ratio_dc05ac2_11,
        avg_ratio_dc00ac1_12,
        avg_ratio_dc05ac2_12
    ]

    f_ratio.write('\\midrule\n')
    line = 'Average'
    for ratio in values_to_write:
        line += '& {:.2f}'.format(ratio)
    f_ratio.write(line + '\\\\\n')
    f_ratio.write('\\bottomrule\n')
    f_ratio.close()

    values_to_write = [
        avg_rmse_pxr24,
        avg_rmse_b44,
        avg_rmse_simple01,
        avg_rmse_simple05,
        avg_rmse_simple10,
        avg_rmse_simple15,
        avg_rmse_simple20,
        avg_rmse_simple25,
        avg_rmse_dc00ac1_11,
        avg_rmse_dc05ac2_11,
        avg_rmse_dc00ac1_12,
        avg_rmse_dc05ac2_12
    ]

    f_error.write('\\midrule\n')
    line = 'Average'
    for err in values_to_write:
        line += '& ' + format_sci(err)
    f_error.write(line + '\\\\\n')
    f_error.write('\\bottomrule\n')
    f_error.close()


if __name__ == '__main__':
    main()
