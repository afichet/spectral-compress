#!/usr/bin/env python3

import os, sys
import OpenImageIO as oiio
import numpy as np
from openexr.spectralexr import SpectralEXR
from radiometry.cmf import CMF


def export_1param(filename: str , img: np.array):
    width, height = img.shape[:2]

    spec = oiio.ImageSpec(width, height, 1, 'float')

    out = oiio.ImageOutput.create(filename)
    out.open(filename, spec)
    out.write_image(img)
    out.close()


def export_2params(filename: str , img: np.array):
    width, height = img.shape[:2]

    foo_fb = np.zeros((width, height, 1))
    fb = np.dstack([img, foo_fb])

    spec = oiio.ImageSpec(width, height, 3, 'float')

    out = oiio.ImageOutput.create(filename)
    out.open(filename, spec)
    out.write_image(fb)
    out.close()


def export_3params(filename: str, img: np.array):
    width, height = img.shape[:2]

    spec = oiio.ImageSpec(width, height, 3, 'float')

    out = oiio.ImageOutput.create(filename)
    out.open(filename, spec)
    out.write_image(img)
    out.close()


def export_material(material_dir: str,
                    material_name: str,
                    cmf: CMF,
                    export_dir: str):
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)

    mat_sRGB = [[ 3.2404542, -0.9692660, 0.0556434],
                [-1.5371385, 1.8760108, -0.2040259],
                [-0.4985314, 0.0415560, 1.0572252]]

    data = np.load(os.path.join(material_dir, material_name))

    wavelengths = (data['channels'].reshape(1, data['channels'].shape[0])[0])
    width, height = data['diffuse'].shape[:2]

    diffuse  = SpectralEXR(width, height, wavelengths, True, False, False)

    diffuse.reflective_image = data['diffuse']
    diffuse.save(os.path.join(export_dir, 'diffuse.exr'))

    specular = SpectralEXR(width, height, wavelengths, True, False, False)

    specular.reflective_image = data['specular']
    specular.save(os.path.join(export_dir, 'specular.exr'))

    fresnel_filename = os.path.join(export_dir, 'fresnel.exr')
    export_1param(fresnel_filename, data['fresnel'])

    normal_filename = os.path.join(export_dir, 'normal.exr')
    export_3params(normal_filename, 0.5 * (1 + data['normal']))

    aniso_filename = os.path.join(export_dir, 'aniso.exr')
    export_1param(aniso_filename, (np.pi/2 + data['aniso']) / (2 * np.pi))

    lobes_filename = os.path.join(export_dir, 'lobes.exr')
    export_2params(lobes_filename, data['lobes'])

    # Geometry
    mat_obj = os.path.join(export_dir, 'geometry.obj')

    with open(mat_obj, 'w') as f:
        for vertex in data['vertices']:
            f.write('v {} {} {}\n'.format(vertex[0], vertex[1], vertex[2]))

        for uv in data['texture_coords']:
            f.write('vt {} {}\n'.format(uv[0], uv[1]))

        for face in data['faces']:
            f.write('f {}/{} {}/{} {}/{}\n'.format(
                face[0] + 1, face[0] + 1,
                face[1] + 1, face[1] + 1,
                face[2] + 1, face[2] + 1))


def main():
    if len(sys.argv) < 3:
        print('Usage:')
        print('------')
        print(sys.argv[0], '<material> <export_dir>')
        exit(0)

    path       = sys.argv[1]
    export_dir = sys.argv[2]

    if not os.path.exists(export_dir):
        os.makedirs(export_dir)

    material_name = os.path.basename(path)
    material_dir  = path[:len(path) - len(material_name)]
    base = os.path.abspath(os.path.dirname(__file__))
    cmf = CMF(os.path.join(base, '..', 'data', 'cmf', 'ciexyz06_2deg.csv'))

    export_material(material_dir, material_name, cmf, export_dir)


if __name__ == '__main__':
    main()
