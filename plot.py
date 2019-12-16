import matplotlib as mpl
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from PIL import Image

F1 = 30 # Axes font size
F2 = 15 # Legend font size
line = 1.5 # Line width
alp = 1 # alpha
scale= 'linear'
colorshift=3

# Setting the font structure

rc = mpl.rcParams # Font structure is called rc now
rc['text.usetex'] = False # Tex fonts
rc['font.family'] = 'sans-serif'
rc['font.serif'].insert(0,'cm') # Default font is computer modern for latex
rc['font.size'] = F1
rc['xtick.labelsize'] = 'small'
rc['ytick.labelsize'] = 'small'
rc['legend.fontsize'] = F2

phi = np.loadtxt("output.txt")
#pi = np.loadtxt("pi.txt")
x = np.loadtxt("grid.txt")
t = np.loadtxt("time.txt")
N_time = len(t)
N_space = len(x)

print (N_time, N_space, np.shape(phi))

images = []

fname_template = "Figures/phi_at_time.{i}.txt"

for i in range(0, np.shape(phi)[0], 1):
    fig = plt.figure(figsize = (12,8))
    fig.suptitle('Time step {} of {}, time = {}'.format(i, np.shape(phi)[0], t[-1]*i/N_time), fontsize=20)
    plt.ylabel('Field')
    plt.xlabel('x')
    plt.xlim([x[0], x[-1]])
    #plt.xlim([75, 125])
    plt.ylim([-0.015, 2.6])
    plt.plot(x,phi[i,:], linewidth = line, color = 'black')
    plt.savefig("Figures/phi%04i.png" %float(i),bbox_inches = 'tight')
    plt.close()
    frame = Image.open("Figures/phi%04i.png" %float(i))
    images.append(frame)
    print("Saved figure ", i, "of", np.shape(phi)[0])
#    plt.savefig('Figures/figname_template.{i}.png'.format(i = i), bbox_inches = 'tight')
images[0].save('Figures/Domain_Wall.gif', save_all=True, append_images=images[1:], duration=75, loop=0)