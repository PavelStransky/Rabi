import os
from PIL import Image
from alive_progress import alive_bar

dir = 'D:/results/Rabi/schnellbruder/3 brothers'

fnamesSphere = []
indexSphere = []

for fname in os.listdir(dir):
    start = fname.find("here_")
    end = fname.find(".png")

    if start == -1 or end == -1:
        continue

    start += 5

    fnamesSphere.append(fname)
    indexSphere.append(float(fname[start:end]))

fnamesSphere = [x for _, x in sorted(zip(indexSphere, fnamesSphere))]

fnamesWigner = []
indexWigner = []

for fname in os.listdir(dir):
    start = fname.find("hers_")
    end = fname.find(".png")

    if start == -1 or end == -1:
        continue

    start += 5

    fnamesWigner.append(fname)
    indexWigner.append(float(fname[start:end]))

fnamesWigner = [x for _, x in sorted(zip(indexWigner, fnamesWigner))]

with alive_bar(len(fnamesWigner), title='Combining images') as bar:
    for i, fname in enumerate(fnamesWigner):
        srcWigner = f'{dir}/{fname}'
        srcSphere = f'{dir}/{fnamesSphere[i]}'
        dst = f'{dir}/{i}.png'

        # Load the images
        imgWigner = Image.open(srcWigner)
        imgSphere = Image.open(srcSphere)

        imgSphere = imgSphere.resize((1000, 1000))
        imgSphere = imgSphere.crop((250, 250, 750, 750))
        # imgWigner = imgWigner.resize((1920, 1080))

        # Create a new image with the combined width
        combined = Image.new("RGBA", (1920, 1080))
        combined.paste(imgWigner, (0, 0))
        combined.paste(imgSphere, (1920 - 550, 1080 - 520))

        # Save the result
        combined.save(dst)
        bar()