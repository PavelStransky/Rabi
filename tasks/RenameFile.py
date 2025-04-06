import os

fnames = []
energies = []

dir = 'D:/results/Rabi/schnellbruder/3 brothers'

for fname in os.listdir(dir):
    start = fname.find("0.5)_") + 5
    end = fname.find(".png")

    fnames.append(fname)
    energies.append(float(fname[start:end]))

fnames = [x for _, x in sorted(zip(energies, fnames))]
print(fnames)

for i, fname in enumerate(fnames):
    src = f'{dir}/{fname}'
    dst = f'{dir}/{i}.png'

    os.rename(src, dst)

