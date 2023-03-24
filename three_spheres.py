import math
import matplotlib.pyplot as plt
import numpy as np
test = 0.025E-04 * math.sqrt(3)

print( test + 0.475E-04)



position_1_y = []
position_2_y = []
position_3_y = []

position_1_x = []
position_2_x = []
position_3_x = []


filename = "target/release/three_sphere_data.txt"

file = open(filename, 'r')


count = 0
for line in file.readlines():
    stringvalues = line.split()

    if count % 3 == 0:
        position_1_y.append(float(stringvalues[1]))
        position_1_x.append(float(stringvalues[0]))
    if count % 3 == 1:
        position_2_y.append(float(stringvalues[1]))
        position_2_x.append(float(stringvalues[0]))
    if count % 3 == 2:
        position_3_y.append(float(stringvalues[1]))
        position_3_x.append(float(stringvalues[0]))

    count+= 1


position_1_y = np.array(position_1_y)
position_2_y = np.array(position_2_y)
position_3_y = np.array(position_3_y)

position_1_x = np.array(position_1_x)
position_2_x = np.array(position_2_x)
position_3_x = np.array(position_3_x)



delta_bot = (position_2_x - position_1_x - 0.05E-04) * 1000
delta_right_y = (position_3_y - position_1_y)
delta_right_x = (position_3_x - position_1_x)
delta_left_y = (position_2_y - position_3_y)
delta_left_x = (position_2_x - position_3_x)


delta_right = []
delta_left = []
for i in range(len(delta_left_y)):
    delta_right.append((math.sqrt(delta_right_x[i]**2 + delta_right_y[i]**2) - 0.05E-04) * 1000)
    delta_left.append((math.sqrt(delta_left_x[i]**2 + delta_left_y[i]**2) - 0.05E-04) * 1000)

delta_right = np.array(delta_right)
delta_left = np.array(delta_left)

plt.figure(1)
plt.plot(position_1_y)
plt.plot(position_2_y)
plt.plot(position_3_y)
plt.plot(delta_left)
plt.plot(delta_right)
plt.plot(delta_bot)
plt.legend(["p1 y_pos","p2 y_pos","p3 y_pos","left bond overlap x1000","right bond overlap  x1000","bot bond overlap  x1000"])
plt.show()



plt.figure(2)
plt.plot(delta_right-delta_left)
plt.legend(["right left dif in overlap"])
plt.show()