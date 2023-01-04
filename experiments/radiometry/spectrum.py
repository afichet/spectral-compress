import numpy as np

class Spectrum:
    def __init__(self, path):
        wl = []
        val = []

        with open(path) as f:
            for l in f:
                wl_el, v_el = [ float(el) for el in l.split(',') ]
                wl.append(wl_el)
                val.append(v_el)

        self.data = np.stack((wl, val), axis=1)

class TabSpectrum(Spectrum):
    def __init__(self, wl, values):
        self.data = np.stack((wl, values), axis=1)
