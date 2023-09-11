from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation as animation
import numpy as np
import time 
import math

st = time.time()

residues = np.array([[0, 0, 1, 1, 2, 2, 3, 3, 3, 2, 1, 0, 0, 1, 2, 3], [1, 0, 0, 1, 1, 0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 3]])
res_prob = residues.copy()

ES = -1.25

E_native = 9 * ES  #Total energy in native state

turns = 2000

Ei = E_native


fig = plt.figure()
ax = fig.add_subplot(111)
protein, = ax.plot(residues[0], residues[1], color = 'black', linewidth = 2, marker = "o", markerfacecolor = 'red')

ax.text(0.5, 1.100, "Protein Unfolding at E = %f units" % ES, bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 5}, transform=ax.transAxes, ha="center")

plt.xlim(-8, 12)
plt.ylim(-8, 12)

E_max = E_native
E_max_turns = 0

#for i in range(turns):
def update(frame):
    global res_prob, E, Ei, residues, E_max, E_max_turns

    cm = np.random.randint(16); #Randomly choosing one residue
    x_shift = [-1, 1, -1, 1]
    y_shift = [-1, -1, 1, 1]
    if cm in (0, 15):
        #End move (Available only for residues 1 and 16)
        r = np.random.randint(4)
        xpt = residues[0][cm] + x_shift[r]
        ypt = residues[1][cm] + y_shift[r]

        if cm == 0:
            m = 0
            if math.dist([xpt, ypt], [residues[0][1], residues[1][1]]) == 1:
                for k in range(0, 16):
                    if math.dist([xpt, ypt], [residues[0][k], residues[1][k]]) == 0:
                        m = 1
                        break

                if m == 0:
                    res_prob[0][cm] = xpt
                    res_prob[1][cm] = ypt
                
        elif cm == 15:
            m = 0
            if math.dist([xpt, ypt], [residues[0][14], residues[1][14]]) == 1:
                for k in range(0, 16):
                    if math.dist([xpt, ypt], [residues[0][k], residues[1][k]]) == 0:
                        m = 1
                        break

                if m == 0:
                    res_prob[0][cm] = xpt
                    res_prob[1][cm] = ypt
                 
    else:
        #Corner Move
        xpt = 0
        ypt = 0
        if residues[0][cm] == residues[0][cm-1] and residues[1][cm] == residues[1][cm+1]:
            xpt = residues[0][cm + 1]
            ypt = residues[1][cm - 1]
        elif residues[0][cm] == residues[0][cm + 1] and residues[1][cm] == residues[1][cm - 1]:
            xpt = residues[0][cm - 1]
            ypt = residues[1][cm + 1]
        
        m = 0
        if math.dist([xpt, ypt], [residues[0][cm + 1], residues[1][cm + 1]]) == 1 and math.dist([xpt, ypt], [residues[0][cm - 1], residues[1][cm - 1]]) == 1:
            for k in range(0, 16):
                if math.dist([xpt, ypt], [residues[0][k], residues[1][k]]) == 0:
                    m = 1
                    break

            if m == 0:
                res_prob[0][cm] = xpt
                res_prob[1][cm] = ypt
            
    #Mapping the non-covalent Interactions
    n = 0
    '''if math.dist([residues[0][0], residues[1][0]], [residues[0][3], residues[1][3]]) == 1:
        n += 1
    if math.dist([residues[0][0], residues[1][0]], [residues[0][11], residues[1][11]]) == 1:
        n += 1
    if math.dist([residues[0][3], residues[1][3]], [residues[0][10], residues[1][10]]) == 1:
        n += 1
    if math.dist([residues[0][2], residues[1][2]], [residues[0][5], residues[1][5]]) == 1:
        n += 1
    if math.dist([residues[0][4], residues[1][4]], [residues[0][7], residues[1][7]]) == 1:
        n += 1
    if math.dist([residues[0][4], residues[1][4]], [residues[0][9], residues[1][9]]) == 1:
        n += 1
    if math.dist([residues[0][8], residues[1][8]], [residues[0][15], residues[1][15]]) == 1:
        n += 1
    if math.dist([residues[0][9], residues[1][9]], [residues[0][14], residues[1][14]]) == 1:
        n += 1
    if math.dist([residues[0][10], residues[1][10]], [residues[0][13], residues[1][13]]) == 1:
        n += 1'''
    
    for m in range(16):
        for m1 in range(m+2, 16):
            if math.dist([residues[0][m], residues[1][m]], [residues[0][m1], residues[1][m1]]) == 1:
                n += 1

    Et = n * ES #Energy for the current configuration (Temporary)
    dE = Et - Ei
    w = np.exp(-dE)
    if w>1:
        Ei = Et
        residues = res_prob.copy()
    else:
        tn = np.random.rand()
        if w>tn:
            Ei = Et
            residues = res_prob.copy()
        else:
            res_prob = residues.copy()

    #Updating the values for the animation
    protein.set_xdata(residues[0])
    protein.set_ydata(residues[1])

ani = animation(fig, update, frames = range(turns))
ani.save('Anim_random_E_-1.25.mp4', writer='ffmpeg', fps=30)

et = time.time()
print("Elapsed time:",et - st, " seconds")