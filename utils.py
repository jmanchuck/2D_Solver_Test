import imageio
import os

images = []
files = sorted(os.listdir('figs'))
print(files)
exit()
for filename in files:
    if filename.endswith('.png'):
        filepath = os.path.join('figs', filename)
        images.append(imageio.imread(filepath))
        print(filepath)

# imageio.mimsave('movie.gif', images, duration=1)

print("Done")