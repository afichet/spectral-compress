#!/usr/bin/env python3

import numpy as np


# Converts a XYZ value to linear sRGB
def XYZ_to_sRGB(xyz):
    mat_sRGB = [[ 3.2404542, -1.5371385, -0.4985314],
                [-0.9692660, 1.8760108, 0.0415560],
                [0.0556434, -0.2040259, 1.0572252]]

    return np.matmul(mat_sRGB, xyz)
    

# Converts a linear sRGB value to XYZ
def sRGB_to_XYZ(rgb):
    mat_XYZ = [[0.4124564, 0.3575761, 0.1804375],
               [0.2126729, 0.7151522, 0.0721750],
               [0.0193339, 0.1191920, 0.9503041]]
    
    return np.matmul(mat_XYZ, rgb)
    

# Converts an XYZ value to Lab
def XYZ_to_Lab(xyz):
    epsilon = 0.008856
    kappa   = 903.3
    coefs   = [0.950489, 1., 1.08840]
    
    xyz_n = [ v / c for v, c in zip(xyz, coefs) ]
    
    f = [ v**(1/3) if v > epsilon else (kappa * v + 16) / 116 for v in xyz_n ]

    Lab = [
        116 * f[1] - 16,
        500 * (f[0] - f[1]),
        200 * (f[1] - f[2])
    ]
    
    return Lab


# Converts a linear sRGB vakye to Lab
def sRGB_to_Lab(rgb):
    return XYZ_to_Lab(sRGB_to_XYZ(rgb))


# Gives the gamma corrected sRGB value from a linear value
def to_sRGB_gamma(C):
    C = np.clip(C, 0., 1.)
    
    if abs(C) < 0.0031308:
        return 12.92 * C
    return 1.055 * C**0.41666 - 0.055


# Gives the linear sRGB value from a gamma corrected value
def from_sRGB_gamma(C):
    if C < 0.04045:
        return C / 12.92
    
    return ((C + 0.055) / 1.055)**2.4



# Computes $Delta_E_{2000}^$ of two Lab values
# https://hajim.rochester.edu/ece/sites/gsharma/ciede2000/ciede2000noteCRNA.pdf
def deltaE2000_Lab(Lab_1, Lab_2):
    L1, a1, b1 = Lab_1
    L2, a2, b2 = Lab_2

    # Compute C_prime_{1, 2} and h_prime_{1, 2}
    C_star_1 = np.sqrt(a1*a1 + b1*b1)
    C_star_2 = np.sqrt(a2*a2 + b2*b2)
    bar_C_star = (C_star_1 + C_star_2) / 2
    
    G = .5 * (1 - np.sqrt(bar_C_star**7 / (bar_C_star**7 + 25**7)))
    a_prime_1 = (1 + G) * a1
    a_prime_2 = (1 + G) * a2

    C_prime_1 = np.sqrt(a_prime_1**2 + b1**2)
    C_prime_2 = np.sqrt(a_prime_2**2 + b2**2)

    h_prime_1_rad, h_prime_2_rad = 0, 0

    if b1 != 0 or a_prime_1 != 0:
        h_prime_1_rad = np.arctan2(b1, a_prime_1)

        if h_prime_1_rad < 0:
            h_prime_1_rad += 2 * np.pi

    if b2 != 0 or a_prime_2 != 0:
        h_prime_2_rad = np.arctan2(b2, a_prime_2)

        if h_prime_2_rad < 0:
            h_prime_2_rad += 2 * np.pi

    # Compute Delta_L_prime, Delta_C_prime, Delta_H_prime
    p_C_prime_12 = C_prime_1 * C_prime_2

    Delta_h_prime_rad = 0

    if p_C_prime_12 != 0:
        if np.abs(h_prime_2_rad - h_prime_1_rad) <= np.pi:
            Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad
        elif h_prime_2_rad - h_prime_1_rad > np.pi:
            Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad - 2*np.pi 
        else:
            Delta_h_prime_rad = h_prime_2_rad - h_prime_1_rad + 2*np.pi

    Delta_L_prime = L2 - L1
    Delta_C_prime = C_prime_2 - C_prime_1
    Delta_H_prime = 2 * np.sqrt(C_prime_1 * C_prime_2) * np.sin(Delta_h_prime_rad / 2)

    # Compute CIE_Delta_E_{2000}
    bar_L_prime = (L1 + L2) / 2
    bar_C_prime = (C_prime_1 + C_prime_2) / 2

    bar_h_prime_rad = 0

    if p_C_prime_12 == 0:
        bar_h_prime_rad = h_prime_1_rad + h_prime_2_rad
    elif np.abs(h_prime_1_rad - h_prime_2_rad) <= np.pi:
        bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad) / 2
    elif h_prime_1_rad + h_prime_2_rad < 2 * np.pi:
        bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad + 2 * np.pi) / 2
    else:
        bar_h_prime_rad = (h_prime_1_rad + h_prime_2_rad - 2 * np.pi) / 2

    T = (1. 
        - 0.17 * np.cos(bar_h_prime_rad - np.pi / 6)
        + 0.24 * np.cos(2 * bar_h_prime_rad)
        + 0.32 * np.cos(3 * bar_h_prime_rad + np.pi / 30)
        - 0.20 * np.cos(4 * bar_h_prime_rad - 7 * np.pi / 20))

    Delta_theta_rad = np.pi / 6 * np.exp(- ((bar_h_prime_rad * 180 / np.pi - 275) / 25)**2)

    R_C = 2 * np.sqrt(bar_C_prime**7 / (bar_C_prime**7 + 25**7))
    S_L = 1 + 0.015 * (bar_L_prime - 50)**2 / np.sqrt(20 + (bar_L_prime - 50)**2)
    S_C = 1 + 0.045 * bar_C_prime
    S_H = 1 + 0.015 * bar_C_prime * T
    R_T = -np.sin(2 * Delta_theta_rad) * R_C

    K_L = 1
    K_C = 1
    K_H = 1

    delta_L_r = Delta_L_prime / (K_L * S_L)
    delta_C_r = Delta_C_prime / (K_C * S_C)
    delta_H_r = Delta_H_prime / (K_H * S_H)

    return np.sqrt(
          delta_L_r * delta_L_r 
        + delta_C_r * delta_C_r 
        + delta_H_r * delta_H_r
        + R_T * delta_C_r * delta_H_r)


# Computes $Delta_E_{2000}^$ of two XYZ values
def deltaE2000_XYZ(xyz1, xyz2):
    Lab1 = XYZ_to_Lab(xyz1)
    Lab2 = XYZ_to_Lab(xyz2)
    
    return deltaE2000_Lab(Lab1, Lab2)


# Computes $Delta_E_{2000}^$ of two linear sRGB values
def deltaE2000_sRGB(rgb1, rgb2):
    Lab1 = sRGB_to_Lab(rgb1)
    Lab2 = sRGB_to_Lab(rgb2)
    
    return deltaE2000_Lab(Lab1, Lab2)