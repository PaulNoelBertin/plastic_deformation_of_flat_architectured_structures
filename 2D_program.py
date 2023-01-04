import matplotlib.pyplot as plt
import scipy
import numpy as np

plt.rcParams['font.family'] = ['CMU Serif', 'Serif']
plt.rcParams['font.size'] = 35



fig,ax = plt.subplots()
ax.grid(linewidth=0.2)
motifs = [0,1,2,3]
matrix_E =[]
matrix_strain = []
for motif in motifs:
    resultE = []
    resultS = []
    for index in range(2,10):
        percent = index*5
        file = str(motif)+"_"+str(percent)+".txt"
        strain, force, extension = [], [] , []
        with open(file, "r") as f_read:
            for line in f_read:
                line = line.strip()  
                if line:
                    si, fi = [float(elt) for elt in line.split(";")] 
                    strain.append(-si), force.append(-fi), extension.append((-si/8))
        max = percent*0.1-0.01
        argmax = 0
        for index2 in range(len(extension)):
            if extension[index2] > max:
                argmax = index2
        new_extension = extension[argmax:]
        new_force = force[argmax:]
        ax.plot(new_extension,new_force,"--",lw=0.5,label = str(motif)+"_"+str(percent))
        A = []
        for index2 in range(len(new_force)):
            A.append([new_extension[index2],1])
        C, residues, rank, s = scipy.linalg.lstsq(A, new_force) #Interpolation in the relaxation phase
        [a,b] = C
        Z = []
        for index2 in range(len(new_force)):
            Z.append(new_extension[index2]*a +b)
        ax.plot(new_extension,Z,lw=2,label = str(motif)+"_"+str(percent))
        resultE.append(a)
        resultS.append(-b/a)
    matrix_E.append(resultE)
    matrix_strain.append(resultS)
plt.show()


list_percent_exp = [0.50*i for i in range(0,10)]
list_percent_fem = [4.50,4.00,3.50,3.00,2.50,2.00,1.50,1.00,0.50,0]
list_percent_timo = [index*0.5 for index in range(10)]

exp_curv = [0.0,0,9.523809524,16.12903226,15.87301587,11.76470588,8.403361345,7.936507937,6.451612903,7.575757576]

list_percent = [0.50*i for i in range(2,10)]

listE_0 = matrix_E[0]
listE_1 = matrix_E[1]
listE_2 = matrix_E[2]
listE_3 = matrix_E[3]
strain_0 = matrix_strain[0]
strain_1 = matrix_strain[1]
strain_2 = matrix_strain[2]
strain_3 = matrix_strain[3]
h1 = 0.0228
h2 = 0.0228
h = h1+h2
m =1
length = 8
delta_strain = [strain_0[i]-strain_3[i] for i in range(length)]
n = [listE_0[i]/listE_3[i] for i in range(length)]
timoshenko = [24*(delta_strain[i])/(h*(14+n[i]+1/n[i])) for i in range(length)]
timoshenko = [0,0] + timoshenko
print("la veuleuru rau aysugx ",timoshenko[-1])
fig, ax1= plt.subplots(1,1, constrained_layout=True, sharex=False, sharey=False)
ax1.grid(linewidth=0.2)
ax1.set_xlim(0,4.5)
ax1.set_ylim(0,25)
ax1.set_xlabel("$ F_{11}-1$ (-)")
ax1.set_ylabel("Curvature ($m^{-1}$)")
ax1.plot(list_percent_exp,exp_curv, '.-', ms=25, lw=4, alpha=1,color="black",label="Experimental values")
ax1.plot(list_percent_timo,timoshenko, '.-', ms=25, lw=4, alpha=1,color="dimgray",label="Timoshenko")
plt.legend(title =  "AD")
plt.show()

fig,ax = plt.subplots()
ax.grid(linewidth=0.2)
ax.set_xlim(0,4.5)
ax.set_ylim(0,6)
ax.set_xlabel("Applied strain(-)")
ax.set_ylabel("Force (N)")
motifs = [0,1,2,3]
matrix_E =[]
matrix_strain = []
for motif in motifs:
    resultE = []
    resultS = []
    index_liste = [9]
    for index in range(10):
        percent = index*5
        file = str(motif)+"_"+str(percent)+".txt"
        strain, force, extension = [], [] , []
        with open(file, "r") as f_read:
            for line in f_read:
                line = line.strip()  
                if line:
                    si, fi = [float(elt) for elt in line.split(";")] 
                    strain.append(-si), force.append(-fi), extension.append((-si/8))
        max = percent*0.1-0.01
        argmax = 0
        for index2 in range(len(extension)):
            if extension[index2] > max:
                argmax = index2
        new_extension = extension[argmax:]
        new_force = force[argmax:]
        ax.plot(extension,force,"--",lw=1,alpha=1,ms=25,label ="D cell", color = "black")
        A = []
        for index2 in range(len(new_force)):
            A.append([new_extension[index2],1])
        C, residues, rank, s = scipy.linalg.lstsq(A, new_force)
        [a,b] = C
        Z = []
        for index2 in range(len(new_force)):
            Z.append(new_extension[index2]*a +b)
        ax.plot(new_extension,Z,ms=25, lw=4, alpha=1,label = "$F=(\epsilon-\epsilon_{pD})C $",color="black")
        resultE.append(a)
        resultS.append(-b/a)
    matrix_E.append(resultE)
    matrix_strain.append(resultS)
plt.legend()
plt.show()

print(matrix_E)
print(matrix_strain)


list_percent_exp = [0.50*i for i in range(0,10)]
list_percent_fem = [4.50,4.00,3.50,3.00,2.50,2.00,1.50,1.00,0.50,0]
list_percent_timo = [index*0.5 for index in range(10)]

exp_curv = [0,0,22.68190891,27.611246,26.972059,22.818083,18.907684,17.395282,17.42243,17.800316]
fem_curv = [11.670524, 10.690384, 10.717929, 11.367359, 12.605646, 14.495322, 16.891645, 17.248292, 1.98586,0]

list_percent = [0.50*i for i in range(2,10)]
len = len(matrix_E[0])
listE_0 = [matrix_E[0][i]/0.0228*0.0012 for i in range(len)]
listE_1 = matrix_E[1]
listE_2 = matrix_E[2]
listE_3 = [matrix_E[3][i]/0.0228*0.0012 for i in range(len)]
strain_0 = matrix_strain[0]
strain_1 = matrix_strain[1]
strain_2 = matrix_strain[2]
strain_3 = matrix_strain[3]
h1 = 0.0228
h2 = 0.0228
h = h1+h2
m =1
length = 8
delta_strain = [strain_0[i]-strain_3[i] for i in range(length)]
n = [listE_0[i]/listE_3[i] for i in range(length)]
timoshenko = [24*(delta_strain[i])/(h*(14+n[i]+1/n[i])) for i in range(length)]
timoshenko = [0,0] + timoshenko
h1 = 0.0228/2
h2 = 0.0228/2
h = h1+h2
m =1
length = 8
delta_strain = [strain_0[i]-strain_3[i] for i in range(length)]
n = [listE_0[i]/listE_3[i] for i in range(length)]
timoshenko2 = [24*(delta_strain[i])/(h*(14+n[i]+1/n[i])) for i in range(length)]
timoshenko2 = [0,0] + timoshenko2

fig, ax1= plt.subplots(1,1, constrained_layout=True, sharex=False, sharey=False)
ax1.grid(linewidth=0.2)
ax1.set_xlim(0,4.5)
ax1.set_ylim(0,50)
ax1.set_xlabel("$ F_{11}-1$ (-)")
ax1.set_ylabel("Curvature ($m^{-1}$)")
ax1.plot(list_percent_exp,exp_curv, '.-', ms=25, lw=4, alpha=1,color="black",label="Experimental values")
#ax1.plot(list_percent_fem,fem_curv, '.-', ms=25, lw=4, alpha=1,color="lightgray",label="FEM")
ax1.plot(list_percent_timo,timoshenko, '.-', ms=25, lw=4, alpha=1,color="dimgray",label="Timoshenko : full length")
ax1.plot(list_percent_timo,timoshenko2, '.-', ms=25, lw=4, alpha=1,color="darkgray",label="Timoshenko : half length")
plt.legend(title = "Hybrid AD")
plt.show()


fig, ax1= plt.subplots(1,1, constrained_layout=True, sharex=False, sharey=False)
ax1.grid(linewidth=0.2)
ax1.set_xlim(1,4.5)
ax1.set_ylim(0,5)
ax1.set_xlabel("Applied strain (-)")
ax1.set_ylabel("Plastic strain (-)")
ax1.plot(list_percent,strain_0, '.-', ms=25, lw=4, alpha=1,color="darkgrey",label="A")
ax1.plot(list_percent,strain_1, '.-', ms=25, lw=4, alpha=1,color="grey",label="B")
ax1.plot(list_percent,strain_2, '.-', ms=25, lw=4, alpha=1,color="dimgrey",label="C")
ax1.plot(list_percent,strain_3, '.-', ms=25, lw=4, alpha=1,color="black",label="D")
plt.legend(title = "Cells")
plt.show()

fig, ax1= plt.subplots(1,1, constrained_layout=True, sharex=False, sharey=False)
ax1.grid(linewidth=0.2)
ax1.set_xlim(1,4.5)
ax1.set_ylim(0,2.1)
ax1.set_xlabel("$F_{11}-1$ (-)")
ax1.set_ylabel("Young's modulus $E$ $(Pa)$")
ax1.plot(list_percent,listE_0, '.-', ms=25, lw=4, alpha=1,color="darkgrey",label="A")
ax1.plot(list_percent,listE_3, '.-', ms=25, lw=4, alpha=1,color="black",label="D")
plt.legend(title = "Cells")
plt.show()

liste = [0,0.10762591479702276,42.32055475622099,195.64303186679388,222.45257217007253,166.79807918908756,158.5460911287254,126.56331615500947,93.61316918791205,104.1803209928674]
list_percent = [index*0.5 for index in range(10)]
exp_curv1 = [0,0,22.68190891,27.611246,26.972059,22.818083,18.907684,17.395282,17.42243,17.800316]
exp_curv2 = [0.0,0,9.523809524,16.12903226,15.87301587,11.76470588,8.403361345,7.936507937,6.451612903,7.575757576]

fig, ax1= plt.subplots(1,1, constrained_layout=True, sharex=False, sharey=False)
ax1.grid(linewidth=0.2)
ax1.set_xlim(1,4.5)
ax1.set_ylim(0,2.1)
ax1.set_xlabel("Applied strain(-)")
ax1.set_ylabel("Young's modulus $E$ $(Pa)$")
ax1.plot(list_percent,liste, '.-', ms=25, lw=4, alpha=1,color="black",label="ADDDA")
ax1.plot(list_percent,exp_curv1, '.-', ms=25, lw=4, alpha=1,color="grey",label="AD hybrid")
ax1.plot(list_percent,exp_curv2, '.-', ms=25, lw=4, alpha=1,color="dimgrey",label="AD")
plt.legend(title = "Cells")
plt.show()
