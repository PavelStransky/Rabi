import os

fnames = []
energies = []

dir = 'D:/results/Rabi/schnellbruder/2 brothers detail lf=-2.0'

for fname in os.listdir(dir):
    start = fname.find("0.5)_")
    # start = fname.find("here_") + 5
    end = fname.find(".png")

    if start == -1 or end == -1:
        continue

    start += 5

    fnames.append(fname)
    energies.append(float(fname[start:end]))

fnames = [x for _, x in sorted(zip(energies, fnames))]
print(fnames)

for i, fname in enumerate(fnames):
    src = f'{dir}/{fname}'
    dst = f'{dir}/brothers_{i+1}.png'
    os.rename(src, dst)

