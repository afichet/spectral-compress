#!/usr/bin/env python3

import OpenImageIO as oiio
import numpy as np
import re

from radiometry.cmf import CMF
import os

class SpectralEXR:
    width = 0
    height = 0

    emissive_wavelengths_nm = [[], [], [], []]
    emissive_images = [None, None, None, None]

    reflective_wavelengths_nm = []
    reflective_image = None

    is_reflective = False
    is_emissive = False
    is_polarised = False

    def __init__(self, width, height, wl_nm, is_reflective, is_emissive, is_polarised):
        self.width, self.height = width, height
        
        self.is_reflective = is_reflective
        self.is_emissive   = is_emissive
        self.is_polarised  = is_polarised

        if is_reflective:
            self.reflective_wavelengths_nm = wl_nm.copy()
            self.reflective_image = np.zeros((height, width, len(wl_nm)))
        
        if is_emissive:
            n_stokes = 4 if is_polarised else 1

            for s in range(n_stokes):
                self.emissive_wavelengths_nm[s] = wl_nm.copy()
                self.emissive_images[s] = np.zeros((height, width, len(wl_nm)))


    def save(self, filename: str):
        d65_wl = [300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0]
        d65_sp = [0.0341, 1.6643, 3.2945, 11.7652, 20.236, 28.6447, 37.0535, 38.5011, 39.9488, 42.4302, 44.9117, 45.775, 46.6383, 49.3637, 52.0891, 51.0323, 49.9755, 52.3118, 54.6482, 68.7015, 82.7549, 87.1204, 91.486, 92.4589, 93.4318, 90.057, 86.6823, 95.7736, 104.865, 110.936, 117.008, 117.41, 117.812, 116.336, 114.861, 115.392, 115.923, 112.367, 108.811, 109.082, 109.354, 108.578, 107.802, 106.296, 104.79, 106.239, 107.689, 106.047, 104.405, 104.225, 104.046, 102.023, 100.0, 98.1671, 96.3342, 96.0611, 95.788, 92.2368, 88.6856, 89.3459, 90.0062, 89.8026, 89.5991, 88.6489, 87.6987, 85.4936, 83.2886, 83.4939, 83.6992, 81.863, 80.0268, 80.1207, 80.2146, 81.2462, 82.2778, 80.281, 78.2842, 74.0027, 69.7213, 70.6652, 71.6091, 72.979, 74.349, 67.9765, 61.604, 65.7448, 69.8856, 72.4863, 75.087, 69.3398, 63.5927, 55.0054, 46.4182, 56.6118, 66.8054, 65.0941, 63.3828]
        
        channels = ['R', 'G', 'B']

        # Create the RGB framebuffer
        base = os.path.abspath(os.path.dirname(__file__))
        cmf = CMF(os.path.join(base, '..', '..', 'data', 'cmf', 'ciexyz06_2deg.csv'))

        mat_sRGB = [[ 3.2404542, -0.9692660, 0.0556434],
                    [-1.5371385, 1.8760108, -0.2040259],
                    [-0.4985314, 0.0415560, 1.0572252]]
        
        framebuffers = np.zeros((self.width, self.height, 3))

        if self.is_emissive:
            framebuffers = cmf.get_xyz_emissive_img(self.emissive_wavelengths_nm[0], self.emissive_images[0])
        else:
            framebuffers = cmf.get_xyz_reflective_img(d65_wl, d65_sp, self.reflective_wavelengths_nm, self.reflective_image)
        
        framebuffers = framebuffers @ mat_sRGB
        
        # Append emissive layers
        if self.is_emissive:
            n_stokes = 4 if self.is_polarised else 1

            for s in range(n_stokes):
                framebuffers = np.dstack((framebuffers, self.emissive_images[s]))

                for wl in self.emissive_wavelengths_nm[s]:
                    channels.append(SpectralEXR.get_channel_name_emissive(s, wl))

        # Append reflective layers
        if self.is_reflective:
            framebuffers = np.dstack((framebuffers, self.reflective_image))

            for wl in self.reflective_wavelengths_nm:
                channels.append(SpectralEXR.get_channel_name_reflective(wl))

        # Write the image
        out = oiio.ImageOutput.create(filename)
        spec = oiio.ImageSpec(self.width, self.height, len(channels), 'float')
        spec.channelnames = channels
        
        out.open(filename, spec)
        out.write_image(framebuffers)
        out.close()


    def print_info(self):
        print('Reflective:', self.is_reflective)
        print('  Emissive:', self.is_emissive)
        print(' Polarised:', self.is_polarised)
        print('Dimensions:', self.width, 'px,', self.height, 'px')


    def get_reflective(self):
        return self.reflective_wavelengths_nm, self.reflective_image


    def get_emissive(self, stokes_component = 0):
        return self.reflective_wavelengths_nm[stokes_component], self.emissive_images[stokes_component]
    
    @staticmethod
    def str_to_nm(value, multiplier: str, units: str):
        v = float(value.replace(',', '.'))

        # avoid rounding errors if the value is already in nanometres
        if multiplier == 'n' and units == 'm':
            return v

        unit_prefix = { 'Y':1e24, 'Z':1e21, 'E':1e18, 'P':1e15, 'T':1e12, 'G':1e9, 'M':1e6, 'k':1e3, 'h':1e2, 'd':1e1, '': 1, 'd':1e-1, 'c':1e-2, 'm':1e-3, 'u':1e-6, 'n':1e-9, 'p':1e-12}

        v *= unit_prefix[multiplier]

        if units == 'Hz':
            v = 299792458. / v * 1e9
        elif units == 'm':
            v *= 1e9
        else:
            raise ValueError

        return v

    @staticmethod
    def get_emissive_channels_idx(im: oiio.ImageInput):
        wl_nm = [[], [], [], []]
        idx = [[], [], [], []]

        # regex = r'^S([0-3])\.*(\d*,?\d*([eE][-+]?\d+)?)(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p|f|a|z|y)?(m|Hz)$'
        regex = r'^S([0-3])\.*(((\d+(,\d*)?)|(,\d+))([eE][+-]?\d+)?)(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p|f|a|z|y)?(m|Hz)$'
        p = re.compile(regex)
        channels = im.spec().channelnames

        for c, i in zip(channels, range(len(channels))):
            m = p.match(c)

            if m:
                stokes_component = int(m.group(1))
                value = m.group(2)
                multiplier, units = m.group(8, 9)

                value = SpectralEXR.str_to_nm(value, multiplier, units)

                wl_nm[stokes_component].append(value)
                idx[stokes_component].append(i)

        return wl_nm, idx

    @staticmethod
    def get_reflective_channels_idx(im: oiio.ImageInput):
        wl_nm = []
        idx = []

        regex = r'^T\.*(((\d+(,\d*)?)|(,\d+))([eE][+-]?\d+)?)(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p|f|a|z|y)?(m|Hz)$'

        p = re.compile(regex)
        channels = im.spec().channelnames

        for c, i in zip(channels, range(len(channels))):
            m = p.match(c)

            if m:
                value = m.group(1)
                multiplier, units = m.group(7, 8)

                value = SpectralEXR.str_to_nm(value, multiplier, units)

                wl_nm.append(value)
                idx.append(i)

        return wl_nm, idx

    @staticmethod
    def get_channel_name_reflective(wavelength_nm):
        return 'T.{}nm'.format(str(wavelength_nm).replace('.', ','))

    @staticmethod
    def get_channel_name_emissive(stokes_component: int, wavelength_nm):
        return 'S{}.{}nm'.format(stokes_component, str(wavelength_nm).replace('.', ','))


class SpectralEXRFile(SpectralEXR):
    def __init__(self, filename: str):
        im = oiio.ImageInput.open(filename)
        px = im.read_image()

        # self.width, self.height = im.spec().width, im.spec().height
        self.height, self.width = px.shape[:2]

        r_wl, r_idx = SpectralEXR.get_reflective_channels_idx(im)
        e_wl, e_idx = SpectralEXR.get_emissive_channels_idx(im)

        self.reflective_wavelengths_nm = r_wl
        self.emissive_wavelengths_nm   = e_wl

        if len(r_idx) > 0:
            self.reflective_image = px[:, :, r_idx]
        else:
            self.reflective_image = None

        for i in range(4):
            if len(e_idx[i]) > 0:
                self.emissive_images[i] = px[:, :, e_idx[i]]
            else:
                self.emissive_images[i] = None

        self.is_reflective = self.reflective_image != None
        self.is_emissive = len(self.emissive_wavelengths_nm[0]) > 0
        self.is_polarised = (
            len(self.emissive_wavelengths_nm[0]) > 0 and 
            len(self.emissive_wavelengths_nm[1]) > 0 and 
            len(self.emissive_wavelengths_nm[2]) > 0 and 
            len(self.emissive_wavelengths_nm[3]) > 0)
