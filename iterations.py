from matplotlib import pyplot as plt
import numpy as np
import time
import math
import seaborn as sns
import pandas as pd

st = time.time()

ES = float(input("Enter Energy: "))

E_max_array = []
E_max_turn_array = []

Rg_max_array = []
Rg_max_turn_array = []

e2e_max_array = []
e2e_max_turn_array = []

E_matrix = np.empty((100000, 0), float)
Rg_matrix = np.empty((100000, 0), float)
e2e_matrix = np.empty((100000, 0), float)

microstates_50 = []

n_sum = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
n_mean = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for molecule in range(50):
    residues = np.array([[0, 0, 1, 1, 2, 2, 3, 3, 3, 2, 1, 0, 0, 1, 2, 3], [1, 0, 0, 1, 1, 0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 3]])
    res_prob = residues.copy()
    res_max_E = residues.copy()
    res_max_Rg = residues.copy()
    res_max_e2e = residues.copy()

    n_states = np.empty((100000, 10), int)

    n = 9
    #print("Initial number of interactions: ", n)
    E_native = n * ES  #Total energy in native state
    #print("Initial interaction energy of the protein: ", E_native)
    turns = 100000


    Ei = E_native
    E = []


    i = 0

    E_max = E_native
    E_max_turns = 0
    yd = [0]

    cmx = sum(residues[0]) / 16
    cmy = sum(residues[1]) / 16
    srg = 0

    #Radius of Gyration
    for k in range(16):
        srg = srg + (math.dist([residues[0][k], residues[1][k]], [cmx, cmy]) ** 2)
    Rg = []
    #Rg = [(srg / 16) ** 0.5]


    Rg_max = (srg / 16) ** 0.5
    Rg_max_turn = 0

    #End-to-End Distance
    e2e = []
    #e2e = [math.dist([residues[0][0], residues[1][0]], [residues[0][15], residues[1][15]])]
    e2e_max = math.dist([residues[0][0], residues[1][0]], [residues[0][15], residues[1][15]])
    e2e_max_turn = 0


    xpt = 0
    ypt = 0

    microstates = []

    n0 = []
    n1 = []
    n2 = []
    n3 = []
    n4 = []
    n5 = []
    n6 = []
    n7 = []
    n8 = []
    n9 = []

    for i in range(turns):
        yd.append(i)
        #IE = ES
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
        nr = 0
        if math.dist([residues[0][0], residues[1][0]], [residues[0][3], residues[1][3]]) == 1:
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
            n += 1

        #If all residues can interact with each other instead of just native interaction
        for m in range(16):
            for m1 in range(m+2, 16):
                if math.dist([residues[0][m], residues[1][m]], [residues[0][m1], residues[1][m1]]) == 1:
                    nr += 1



        Et = nr * ES #Energy for the current configuration (Temporary)
        
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
        E.append(Ei)

        #Finding the turn where all native interactions are lost
        '''if n == 0 and ut != 0:
            n_turn = turns
            ut = 0'''

        #Finding interaction energy with the least number of interactions
        if E_max < Ei:
            E_max = Ei
            E_max_turns = i


        #Calculating centre of mass
        cmx = sum(residues[0]) / 16
        cmy = sum(residues[1]) / 16

        #Calculating Radius of Gyration
        srg = 0
        for k in range(16):
            srg = srg + (math.dist([residues[0][k], residues[1][k]], [cmx, cmy]) ** 2)
        Rg.append((srg / 16) ** 0.5)


        #Checking for maximum Radius of Gyration
        if Rg[i] > Rg_max:
            Rg_max = Rg[i]
            Rg_max_turn = i
            res_max_Rg = residues.copy()

        
        e2e.append(math.dist([residues[0][0], residues[1][0]], [residues[0][15], residues[1][15]]))
        
        #Checking for maximum end-to-end distance
        if e2e[i] > e2e_max:
            e2e_max = e2e[i]
            e2e_max_turn = i
            res_max_e2e = residues.copy()

        if E != 0:
            Score = E[i] * e2e[i] * Rg[i]
        else:
            Score = e2e[i] * Rg[i]
        microstates.append(Score)

        if n == 0:
            n0.append(Score)
        elif n == 1:
            n1.append(Score)
        elif n == 2:
            n2.append(Score)
        elif n == 3:
            n3.append(Score)
        elif n == 4:
            n4.append(Score)
        elif n == 5:
            n5.append(Score)
        elif n == 6:
            n6.append(Score)
        elif n == 7:
            n7.append(Score)
        elif n == 8:
            n8.append(Score)
        elif n == 9:
            n9.append(Score)


    
    n_sum[0] = n_sum[0] + len(set(n0))
    n_sum[1] = n_sum[1] + len(set(n1))
    n_sum[2] = n_sum[2] + len(set(n2))
    n_sum[3] = n_sum[3] + len(set(n3))
    n_sum[4] = n_sum[4] + len(set(n4))
    n_sum[5] = n_sum[5] + len(set(n5))
    n_sum[6] = n_sum[6] + len(set(n6))
    n_sum[7] = n_sum[7] + len(set(n7))
    n_sum[8] = n_sum[8] + len(set(n8))
    n_sum[9] = n_sum[9] + len(set(n9))

    #Considering that the protein is unfolded when all native interactions are lost, storing the unfolding turns in an array
    E_max_array.append(E_max)    
    E_max_turn_array.append(E_max_turns)

    microstates_50.append(len(set(microstates)))

    E_matrix = np.append(E_matrix, np.array([E]).T, axis=1)
    Rg_matrix = np.append(Rg_matrix, np.array([Rg]).T, axis=1)
    e2e_matrix = np.append(e2e_matrix, np.array([e2e]).T, axis=1)


    '''Rg_max_array.append(Rg_max)
    Rg_max_turn_array.append(Rg_max_turn)

    e2e_max_array.append(e2e_max)
    e2e_max_turn_array.append(e2e_max_turn)'''

'''print("Maximum Energy Turns: ", E_max_turn_array)
print("Maximum Rg Turns: ", Rg_max_turn_array)
print("Maximum E2E distance turns: ", e2e_max_turn_array)
'''
#plt.hist(E_max_turn_array)
#plt.xlabel('Turns')

'''sns.histplot(E_max_turn_array, binwidth=10, kde=True)

plt.show()'''

plt.plot(E_max_turn_array)
plt.show()
print("Mean unfolding turns: ", np.mean(E_max_turn_array))
print("Median unfolding turns: ", np.median(E_max_turn_array))


for mi in range(10):
    n_mean[mi] = n_sum[mi] / 50

delG_arr = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
smarr = sum(n_mean)
for i in range(10):
    W = n_mean[i] / smarr
    delG_arr[i] = (i * (-0.25)) - np.log(W)



et = time.time()
print("Elapsed Time: ", et - st, "seconds")