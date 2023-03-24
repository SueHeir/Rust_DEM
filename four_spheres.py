import math
import matplotlib.pyplot as plt
import numpy as np
# test = 0.025E-04 * math.sqrt(3)

# print( test + 0.5E-04)



position_1_y = []
position_2_y = []
position_3_y = []
position_4_y = []

position_1_x = []
position_2_x = []
position_3_x = []
position_4_x = []


filename = "target/release/data.txt"

file = open(filename, 'r')


count = 0
for line in file.readlines():
    stringvalues = line.split()

    if count % 4 == 0:
        position_1_y.append(float(stringvalues[1]))
        position_1_x.append(float(stringvalues[0]))
    if count % 4 == 1:
        position_2_y.append(float(stringvalues[1]))
        position_2_x.append(float(stringvalues[0]))
    if count % 4 == 2:
        position_3_y.append(float(stringvalues[1]))
        position_3_x.append(float(stringvalues[0]))
    if count % 4 == 3:
        position_4_y.append(float(stringvalues[1]))
        position_4_x.append(float(stringvalues[0]))

    count+= 1


position_1_y = np.array(position_1_y)
position_2_y = np.array(position_2_y)
position_3_y = np.array(position_3_y)
position_4_y = np.array(position_4_y)

position_1_x = np.array(position_1_x)
position_2_x = np.array(position_2_x)
position_3_x = np.array(position_3_x)
position_4_x = np.array(position_4_x)



delta_left = (position_2_y - position_1_y - 0.05E-04) * 1000
delta_right = (position_3_y - position_4_y - 0.05E-04) * 1000


delta_top = (position_3_x - position_2_x - 0.05E-04) * 1000
delta_bottom = (position_4_x - position_1_x - 0.05E-04) * 1000


plt.figure(1)
plt.plot(position_1_y)
plt.plot(position_2_y)
plt.plot(position_3_y)
plt.plot(position_4_y)
plt.plot(delta_left)
plt.plot(delta_right)
plt.plot(delta_top)
plt.plot(delta_bottom)
plt.legend(["py1","py2","py3","py4","left","right","top","bot"])
plt.show()


plt.figure(2)
plt.plot(position_1_y)
plt.plot(position_2_y)
plt.plot(delta_left)
plt.plot(delta_top)
plt.legend(["py1","py2","left","top"])
plt.show()

plt.figure(3)
plt.plot(position_3_y)
plt.plot(position_4_y)
plt.plot(delta_right)
plt.plot(delta_bottom)
plt.legend(["py1","py2","right","bot"])
plt.show()


plt.figure(4)
plt.plot(delta_right-delta_left)
plt.plot(delta_bottom-delta_top)
plt.legend(["right left dif","top bot dif"])
plt.show()