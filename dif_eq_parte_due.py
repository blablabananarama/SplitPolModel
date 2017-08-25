import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from reactionsCls.py import rxn

def rxn(Z,t):

    lenf1 = 800
    lenf2 = 900
    ktr = 50 #transcription rate in nt/s^-1
    ktl = 15 #translation rate in nt/s^-1
    ktdeg = 10 #transcript deg
    kdeg = 0.28 #Protein deg
    krapa = 1.2   #T7 ass const
    kd1 = 0.2 * 10 ** -9
    kd2 = 26 * 10 ** -6
    kd3 = 12 * 10 ** -9
    kd4 = 100 * 10 ** -12



    #r1 = k5 * Z[4] * Z[6]                  #T7 association
    #r2 = k6 * Z[2]                         #T7 dissociation
    F1tr = ktr * Z[1] * Z[2]                #transcription
    F1tdeg = ktdeg * Z[4]                   #degradation of transcript
    F1tl = ktl * Z[4]                       #translation
    aF1R = kd1 * krapa * Z[5] * Z[0]        #association of F1 with Rapa
    dF1R = kd1 * krapa * Z[6]               #dissociation of F1 with Rapa
    F1deg = kdeg * Z[5]                     #deg of protein
    F1Rdeg = kdeg * Z[6]                    #deg of F1 with Rapa
    F2tr = ktr * Z[1] * Z[2]                #transcription
    F2tdeg = ktdeg * Z[7]                   #degradation of transcript
    F2tl = ktl * Z[7]                       #translation
    aF2R = kd2 * krapa * Z[8] * Z[0]        #association of F2 with Rapa
    dF2R = kd2 * krapa * Z[9]               #dissociation of F2 with Rapa
    F2deg = kdeg * Z[8]                     #deg of Protein
    F2Rdeg = kdeg * Z[9]
    aF1F2R1 = kd3 * krapa * Z[6] * Z[8]
    dF1F2R1 = kd3 * krapa * Z[3]
    aF1F2R2 = kd4 * krapa * Z[9] * Z[5]
    dF1F2R2 = kd4 * krapa * Z[3]

    dRapdt = Z[0]
    dGendt = Z[1]
    dPoldt = Z[2]
    dF1F2Rdt = aF1F2R1 + aF1F2R2 - dF1F2R1 - dF1F2R2
    dF1trdt = F1tr - F1tdeg
    dF1dt = F1tl - aF1R - F1deg + dF1R
    dF1Rdt = aF1R - dF1R - F1Rdeg - aF1F2R1 + dF1F2R1
    dF2trdt = F2tr - F2tdeg
    dF2dt = F2tl - aF2R - F2deg + dF2R
    dF2Rdt = aF2R - dF2R - F2Rdeg - aF1F2R2 + dF1F2R2


    return(dRapdt, dGendt, dPoldt, dF1F2Rdt, dF1trdt, dF1dt, dF1Rdt, dF2trdt, dF2dt, dF2Rdt)

t = np.linspace(0,2,10000)
Z0 = [1,0.01,0.1,0,0,0,0,0,0,0]
C = odeint(rxn,Z0,t)



plt.plot(t,C[:,3],color = 'b',linewidth = 2)
plt.plot(t,C[:,4],color = 'r',linewidth = 2)
plt.plot(t,C[:,5],color = 'c',linewidth = 2)
plt.plot(t,C[:,6],color = 'y',linewidth = 2)
plt.plot(t,C[:,7],color = 'g',linewidth = 2)
plt.plot(t,C[:,8],color = 'm',linewidth = 2)
plt.plot(t,C[:,9],color = 'k',linewidth = 2)
plt.xlabel('Time [s]')
plt.ylabel('concentration')
plt.legend(['F1F2R', 'F1 transcript', 'F1', 'F1R', 'F2 transcript', 'F2', 'F2R'])
plt.savefig('/home/jonas/Documents/python/output_models/t7plot.png')
plt.show()
