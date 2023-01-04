import numpy as np
import scipy.interpolate

class CMF:
    def __init__(self, filename):
        x_bar = []
        y_bar = []
        z_bar = []
        
        # Read the CSV file
        with open(filename) as f:    
            for l in f:
                wl, x, y, z = [float(el) for el in l.split(',')]
                x_bar.append([wl, x])
                y_bar.append([wl, y])
                z_bar.append([wl, z])

        x_bar = np.array(x_bar)
        y_bar = np.array(y_bar)
        z_bar = np.array(z_bar)
        
        # Interpolate every 1nm to ease computations
        self.wl_start = np.min(x_bar[:, 0])
        self.wl_end   = np.max(x_bar[:, 0])
        # self.wl       = np.linspace(self.wl_start, self.wl_end, num=int(self.wl_end - self.wl_start + 1))

        # let's be a little bit less hardcore here...
        self.wl       = np.arange(self.wl_start, self.wl_end, 10)

        x_bar_y = np.interp(self.wl, x_bar[:, 0], x_bar[:, 1], left=0, right=0)
        y_bar_y = np.interp(self.wl, y_bar[:, 0], y_bar[:, 1], left=0, right=0)
        z_bar_y = np.interp(self.wl, z_bar[:, 0], z_bar[:, 1], left=0, right=0)

        self.x_bar = np.stack((self.wl, x_bar_y), axis=1)
        self.y_bar = np.stack((self.wl, y_bar_y), axis=1)
        self.z_bar = np.stack((self.wl, z_bar_y), axis=1)
        
    
    # Get XYZ values for a given emissive spectrum
    def get_xyz_emissive(self, spectrum_wl, spectrum_values):                   
        # Interpolate spectrum to match the CMFs datapoints
        spectrum_values = np.interp(self.wl, spectrum_wl, spectrum_values, left=0, right=0)
        
        spectrum_x = np.trapz([ x_b * s for x_b, s in zip(self.x_bar[:, 1], spectrum_values) ], self.wl)
        spectrum_y = np.trapz([ y_b * s for y_b, s in zip(self.y_bar[:, 1], spectrum_values) ], self.wl)
        spectrum_z = np.trapz([ z_b * s for z_b, s in zip(self.z_bar[:, 1], spectrum_values) ], self.wl)
    
        return spectrum_x, spectrum_y, spectrum_z


    def get_sRGB_lin_emissive(self, spectrum_wl, spectrum_values):
        xyz = self.get_xyz_emissive(spectrum_wl, spectrum_values)

        mat_sRGB = [
            [ 3.2404542, -0.9692660,  0.0556434],
            [-1.5371385,  1.8760108, -0.2040259],
            [-0.4985314,  0.0415560,  1.0572252]]

        return xyz @ mat_sRGB


    def get_xyz_emissive_img(self, spectrum_wl, spectral_image):
        interp_image_f = scipy.interpolate.interp1d(spectrum_wl, spectral_image, bounds_error=False, fill_value=(0, 0), )
        interp_image = np.maximum(interp_image_f(self.wl), 0)

        s_x = np.trapz(interp_image * self.x_bar[:, 1], x=self.wl)
        s_y = np.trapz(interp_image * self.y_bar[:, 1], x=self.wl)
        s_z = np.trapz(interp_image * self.z_bar[:, 1], x=self.wl)

        return np.dstack((s_x, s_y, s_z))


    def get_sRGB_lin_emissive_img(self, spectrum_wl, spectral_image):
        image_xyz = self.get_xyz_emissive_img(spectrum_wl, spectral_image)

        mat_sRGB = [
            [ 3.2404542, -0.9692660,  0.0556434],
            [-1.5371385,  1.8760108, -0.2040259],
            [-0.4985314,  0.0415560,  1.0572252]]

        return image_xyz @ mat_sRGB


    def get_xyz_reflective_img(self, wavelength_illu, spectrum_illu, image_wl, spectral_image):
        """
        Converts an reflective spectral image stored in a numpy array 
        (y, x, bands) to a colour image.
        :param wavelnegth_illu list of the wavlengths of the illuminant 
               spectrum
        :param spectrum_illu list of the corresponding radiance of the 
               illuminant spectrum
        :param image_wl list containing the wavelength values in nanometers
               for the provided spectral image
        :param spectral_image numpy array storing the spectral image
        """
        illu_values = np.interp(self.wl, wavelength_illu, spectrum_illu, left=0, right=0)

        interp_image_f = scipy.interpolate.interp1d(image_wl, spectral_image, bounds_error=False, fill_value=(0, 0))
        interp_image = np.maximum(interp_image_f(self.wl), 0)

        Y_illu = np.trapz(illu_values * self.y_bar[:, 1], x=self.wl)

        s_x = np.trapz(interp_image * illu_values * self.x_bar[:, 1], x=self.wl) / Y_illu
        s_y = np.trapz(interp_image * illu_values * self.y_bar[:, 1], x=self.wl) / Y_illu
        s_z = np.trapz(interp_image * illu_values * self.z_bar[:, 1], x=self.wl) / Y_illu

        return np.dstack((s_x, s_y, s_z))


    def get_sRGB_lin_reflective_img(self, wavelength_illu, spectrum_illu, image_wl, spectral_image):
        image_xyz = self.get_xyz_reflective_img(wavelength_illu, spectrum_illu, image_wl, spectral_image)

        mat_sRGB = [
            [ 3.2404542, -0.9692660,  0.0556434],
            [-1.5371385,  1.8760108, -0.2040259],
            [-0.4985314,  0.0415560,  1.0572252]]

        return image_xyz @ mat_sRGB