import os
names = []


for f in os.listdir(os.path.join('..', 'data', 'spectra', 'illuminants')):
    str_test_init = ''
    with open(os.path.join('..', 'data', 'spectra', 'illuminants', f), 'r') as filein:
        name = f[11:-4].replace(' ', '_').replace('.', '_').replace('-', '_')
        
        wavelengths = []
        values = []

        for line in filein:
            wl, v = line.split(',')
            wavelengths.append(float(wl))
            values.append(str(float(v)))


        print('const float {}_start = {};'.format(name, wavelengths[0]))
        print('const float {}_end = {};'.format(name, wavelengths[-1]))
        print('const int {}_n_samples = {};'.format(name, len(wavelengths)))
        print('const std::vector<double> {} = {{{}}};'.format(name, ', '.join(values)))
        print()

        names.append(name)

    
for name in names:
    print('Spectrum s_{name} = {{ false, {name}_start, {name}_end, {name}_n_samples, {name}, "{name}" }};'.format(name=name))


print(',\n'.join(names))