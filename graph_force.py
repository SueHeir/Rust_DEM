import matplotlib.pyplot as plt
import numpy as np


openfile = open("force.txt", 'r')

line = openfile.readline()


displacement = []
force = []
tangential_displacement = []
tangential_force = []


un_displacement = []
un_force = []
un_tangential_displacement = []
un_tangential_force = []

while line:
    stringavalues = line.split()
    k_loading = float(stringavalues[0])
    if k_loading == 0:
        displacement.append(float(stringavalues[1]))
        force.append(float(stringavalues[2]) * 1000)
        tangential_displacement.append(float(stringavalues[3])* 1000000 )
        tangential_force.append(float(stringavalues[5])* 1000)
        line = openfile.readline()
    if k_loading == 1:
        un_displacement.append(float(stringavalues[1]))
        un_force.append(float(stringavalues[2]) * 1000)
        un_tangential_displacement.append(float(stringavalues[3])* 1000000)
        un_tangential_force.append(float(stringavalues[5])* 1000)
        line = openfile.readline()

plt.figure(1)
plt.plot(displacement,force)
plt.plot(un_displacement,un_force)
plt.show()


plt.figure(2)
plt.plot(force,tangential_force)
plt.plot(un_force,un_tangential_force)
plt.show()



plt.figure(3)
plt.plot(tangential_displacement ,tangential_force)
plt.plot(un_tangential_displacement ,un_tangential_force)
plt.show()