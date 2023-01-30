#!/usr/bin/env python3

import numpy as np
from radiometry.spectrum import Spectrum, TabSpectrum


class FluoSpectrum:
    wavelength_i_start = 0
    wavelength_i_end = 0
    wavelength_i_sampling = 1
    wavelength_i_n_samples = 0
    
    wavelength_o_start = 0
    wavelength_o_end = 0
    wavelength_o_sampling = 1
    wavelength_o_n_samples = 0

    data = np.zeros((wavelength_i_n_samples, wavelength_o_n_samples, 1))

    def __init__(self, filename):
        tmp_data = []

        with open(filename, 'r') as f:
            if filename.endswith('.bfc') or filename.endswith('.BFC'):

                line_number = 0
                for line in f:
                    line_number += 1

                    end_of_data = line.startswith('EOD') # Last line of BFC

                    if line_number < 11:
                        pass
                    elif line_number == 11:
                        # We parse the sampling and boundaries
                        (self.wavelength_o_start, 
                         self.wavelength_o_end, 
                         self.wavelength_o_sampling,
                         self.wavelength_i_n_samples,
                         self.wavelength_i_start,
                         self.wavelength_i_sampling) = [int(el) for el in line.split()]

                        self.wavelength_o_n_samples = int((self.wavelength_o_end - self.wavelength_o_start) / self.wavelength_o_sampling) + 1
                        self.wavelength_i_end = self.wavelength_i_start + (self.wavelength_i_n_samples - 1) * self.wavelength_i_sampling

                        # print('i: start = {}, end = {}, sampling = {}, n_samples = {}'.format(
                        #      self.wavelength_i_start, self.wavelength_i_end, self.wavelength_i_sampling, self.wavelength_i_n_samples))
                        # print('o: start = {}, end = {}, sampling = {}, n_samples = {}'.format(
                        #      self.wavelength_o_start, self.wavelength_o_end, self.wavelength_o_sampling, self.wavelength_o_n_samples))

                    elif line_number > 12 and not end_of_data:
                        # Populate the data
                        read_data = [float(el) for el in line.split()[1:]]

                        # This avoids a bug in file having an empty line in the middle of data
                        if (len(read_data) > 0):
                            tmp_data.append(read_data)

                self.data = np.reshape(tmp_data, (self.wavelength_o_n_samples, self.wavelength_i_n_samples))

            else:
                header_pos = 0

                for line in f:

                    # Read header info
                    if line.startswith('#'):
                        a, b = line[1:].split()
                        if header_pos == 0:
                            self.wavelength_i_n_samples = int(a)
                            self.wavelength_o_n_samples = int(b)
                        elif header_pos == 1:
                            self.wavelength_i_start = float(a)
                            self.wavelength_i_sampling = float(b)
                        elif header_pos == 2:
                            self.wavelength_o_start = float(a)
                            self.wavelength_o_sampling = float(b)
                        header_pos += 1
                    else:
                        tmp_data.append([float(el) for el in line.split()])

                self.data = np.reshape(tmp_data, (self.wavelength_o_n_samples, self.wavelength_i_n_samples))
                
                self.wavelength_i_end = self.wavelength_i_start + (self.wavelength_i_n_samples - 1) * self.wavelength_i_sampling
                self.wavelength_o_end = self.wavelength_o_start + (self.wavelength_o_n_samples - 1) * self.wavelength_o_sampling
                

        start_wl = max(self.wavelength_i_start, self.wavelength_o_start)
        start_wl_i_idx = self.idx_for_wl_in(start_wl)
        start_wl_o_idx = self.idx_for_wl_out(start_wl)

        # Surface reflectance (diagonal values) can not be negative
        for idx_i, idx_o in zip(range(start_wl_i_idx, self.wavelength_i_n_samples), 
                                range(start_wl_o_idx, self.wavelength_o_n_samples)):
            if self.data[idx_o, idx_i] < 0:
                self.data[idx_o, idx_i] = 0
                
    def export_bfc(self, filename):
        with open(filename, 'w') as f:
            f.write('VEC_01	5167\n')
            f.write('BFC-450 Matrix File\n')
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                int(self.wavelength_o_start),
                int(self.wavelength_o_end),
                int(self.wavelength_o_sampling),
                int(self.wavelength_i_n_samples),
                int(self.wavelength_i_start),
                int(self.wavelength_i_sampling)))
                    
            f.write('r:c:\t')
            
            output_wl = [int(w) for w in self.get_reflectance_wavelengths()]
            input_wl  = [int(w) for w in self.get_excitation_wavelengths()]
            
            for wl in range(input_wl[0], input_wl[-1], int(self.wavelength_i_sampling)):
                f.write('{} '.format(wl))
            
            f.write('{}\n'.format(input_wl[-1]))
            
            for wl_o, idx_o in zip(output_wl, range(len(output_wl))):
                f.write('{} '.format(wl_o))
                
                for wl_i, idx_i in zip(input_wl, range(len(input_wl))):
                    f.write('{}'.format(self.data[idx_o, idx_i]))
                    
                    if idx_i < len(input_wl) - 1:
                        f.write('\t')
                        
                f.write('\n')
                
            f.write('EOD')

    def get_pure_fluo(self):
        if self.wavelength_i_sampling != self.wavelength_o_sampling:
            print('We do not support manipulation of reradiation matrix when')
            print('different sampling precisions are used for input and output!')

        start_wl = max(self.wavelength_i_start, self.wavelength_o_start)

        start_wl_i_idx = self.idx_for_wl_in(start_wl)
        start_wl_o_idx = self.idx_for_wl_out(start_wl)

        ret_array = self.data.copy()

        # Zeroing out the diagonal
        for idx_i, idx_o in zip(range(start_wl_i_idx, self.wavelength_i_n_samples), 
                                range(start_wl_o_idx, self.wavelength_o_n_samples)):
            ret_array[idx_o, idx_i] = 0

        return ret_array


    def get_pure_fluo_filtered(self):
        fluo = self.get_pure_fluo()

        for o in range(fluo.shape[0]):
            for i in range(fluo.shape[1]):
                # Zeroing out invalid negative values
                if fluo[o, i] < 0:
                    fluo[o, i] = 0

                # Zeroing out invalid (nonzero) values under the diagonal
                if self.wl_for_idx_in(i) >= self.wl_for_idx_out(o):
                    fluo[o, i] = 0

        return fluo
    
    
    # Returns diagonal (keep it there for compatibility, get_non_fluo_spectrum shall be used instead)
    def get_non_fluo(self):
        if self.wavelength_i_sampling != self.wavelength_o_sampling:
            print('We do not support manipulation of reradiation matrix when')
            print('different sampling precisions are used for input and output!')

        start_wl = max(self.wavelength_i_start, self.wavelength_o_start)

        start_wl_i_idx = self.idx_for_wl_in(start_wl)
        start_wl_o_idx = self.idx_for_wl_out(start_wl)

        ret_array = []

        for idx_i, idx_o in zip(range(start_wl_i_idx, self.wavelength_i_n_samples), 
                                range(start_wl_o_idx, self.wavelength_o_n_samples)):
            ret_array.append(self.data[idx_o, idx_i])

        return ret_array

    
    # Returns the diagonal with the associated wavelengths
    def get_non_fluo_spectrum(self) -> Spectrum:
        if self.wavelength_i_sampling != self.wavelength_o_sampling:
            print('We do not support manipulation of reradiation matrix when')
            print('different sampling precisions are used for input and output!')

        start_wl = max(self.wavelength_i_start, self.wavelength_o_start)

        start_wl_i_idx = self.idx_for_wl_in(start_wl)
        start_wl_o_idx = self.idx_for_wl_out(start_wl)

        wavelengths = []
        values = []

        for idx_i, idx_o in zip(range(start_wl_i_idx, self.wavelength_i_n_samples), 
                                range(start_wl_o_idx, self.wavelength_o_n_samples)):
            wavelengths.append(self.wl_for_idx_in(idx_i))
            values.append(self.data[idx_o, idx_i])

        return TabSpectrum(wavelengths, values)
    

    def get_pure_reradiation_spectrum(self, wl_i):
        wl_idx_i = self.idx_for_wl_in(wl_i)
        refl_fluo_spectrum = []

        for wl_idx_o in range(self.data.shape[0]):
            wl_o = self.wl_for_idx_out(wl_idx_o)

            # Discard anything under the diagonal including the diagonal
            if wl_o > wl_i:
                refl_fluo_spectrum.append(self.data[wl_idx_o, wl_idx_i])
            else:
                refl_fluo_spectrum.append(0)

        return refl_fluo_spectrum


    def get_pure_absorption_spectrum(self, wl_o):
        wl_idx_o = self.idx_for_wl_in(wl_o)
        absorption_fluo_spectrum = []

        for wl_idx_i in range(self.data.shape[1]):
            wl_i = self.wl_for_idx_out(wl_idx_i)

            # Discard anything over the diagonal including the diagonal
            if wl_i > wl_o:
                absorption_fluo_spectrum.append(self.data[wl_idx_o, wl_idx_i])
            else:
                absorption_fluo_spectrum.append(0)

        return absorption_fluo_spectrum


    def get_average_absorption_spectrum(self):
        return np.sum(self.get_pure_fluo_filtered(), axis=0) / self.wavelength_o_n_samples


    def get_average_reemission_spectrum(self):
        return np.sum(self.get_pure_fluo_filtered(), axis=1) / self.wavelength_i_n_samples


    def gnuplot_reradiation(self, output):
        rerad = self.get_pure_fluo_filtered()
        
        with open(output, "w") as f:
            for i in range(rerad.shape[0]):
                for o in range(rerad.shape[1]):
                    f.write('{} '.format(rerad[i, o]))

                f.write('\n')


    def gnuplot_non_fluo(self, output):
        # Ensure the sampling are the same
        if self.wavelength_i_sampling != self.wavelength_o_sampling:
            raise RuntimeError("Cannot export file with different sampling precision for wl_i and wl_o!")

        # Get the starting point
        wl_start = max(self.wavelength_i_start, self.wavelength_o_start)
        wl_end   = min(self.wavelength_i_end, self.wavelength_o_end)

        idx_wl_i = self.idx_for_wl_in(wl_start)
        idx_wl_o = self.idx_for_wl_out(wl_start)

        curr_wl = wl_start

        with open(output, 'w') as f:
            while curr_wl <= wl_end:
                f.write('{} {}\n'.format(curr_wl, self.data[idx_wl_o, idx_wl_i]))
                curr_wl += self.wavelength_i_sampling
                idx_wl_o += 1
                idx_wl_i += 1


    def gnuplot_average_absorption(self, output):
        absorption_spectrum = self.get_average_absorption_spectrum()

        with open(output, 'w') as f:
            for wl, v in zip(self.get_excitation_wavelengths(), absorption_spectrum):
                f.write('{} {}\n'.format(wl, v))

    
    def gnuplot_average_reemission(self, output):
        reemission_spectrum = self.get_average_reemission_spectrum()

        with open(output, 'w') as f:
            for wl, v in zip(self.get_reflectance_wavelengths(), reemission_spectrum):
                f.write('{} {}\n'.format(wl, v))


    def tikz_reradiation(self, output):
        rerad = self.get_pure_fluo()

        with open(output, "w") as f:
            for i in range(rerad.shape[1]):
                wl_i = self.wl_for_idx_in(i)

                for o in range(rerad.shape[0]):
                    wl_o = self.wl_for_idx_out(o)
                
                    if rerad[o, i] > 0:
                        f.write('{} {} {}\n'.format(wl_i, wl_o, rerad[o, i]))
                    else:
                        f.write('{} {} 0\n'.format(wl_i, wl_o))

                f.write('\n')


    def get_excitation_wavelengths(self):
        return np.linspace(self.wavelength_i_start, self.wavelength_i_end, self.wavelength_i_n_samples)


    def get_reflectance_wavelengths(self):
        return np.linspace(self.wavelength_o_start, self.wavelength_o_end, self.wavelength_o_n_samples)


    def idx_for_wl_in(self, wl_in):
        return int(np.floor(wl_in - self.wavelength_i_start) // self.wavelength_i_sampling)


    def idx_for_wl_out(self, wl_out):
        return int(np.floor(wl_out - self.wavelength_o_start) // self.wavelength_o_sampling)


    def wl_for_idx_in(self, wl_idx_in):
        return self.wavelength_i_start + wl_idx_in * self.wavelength_i_sampling
    

    def wl_for_idx_out(self, wl_idx_out):
        return self.wavelength_o_start + wl_idx_out * self.wavelength_o_sampling



class RSFluoSpectrum(FluoSpectrum):
    def __init__(self, 
        wl_i_start, wl_i_sampling,
        wl_o_start, wl_o_sampling, 
        data):

        self.wavelength_i_n_samples = data.shape[1]
        self.wavelength_i_sampling = wl_i_sampling
        self.wavelength_i_start = wl_i_start
        self.wavelength_i_end = self.wavelength_i_start + (self.wavelength_i_n_samples - 1) * self.wavelength_i_sampling

        self.wavelength_o_n_samples = data.shape[0]
        self.wavelength_o_sampling = wl_o_sampling
        self.wavelength_o_start = wl_o_start
        self.wavelength_o_end = self.wavelength_o_start + (self.wavelength_o_n_samples - 1) * self.wavelength_o_sampling

        self.data = data