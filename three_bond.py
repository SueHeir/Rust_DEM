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


filename = "target/release/bond_data.txt"

file = open(filename, 'r')


count = 0
for line in file.readlines():
    stringvalues = line.split()
    i = int(stringvalues[0])
    j = int(stringvalues[1])

    if i == 0 and j == 1:
        position_1_y.append(float(stringvalues[2]))
        position_1_x.append(float(stringvalues[3]))
    if i == 1 and j == 2:
        position_2_y.append(float(stringvalues[2]))
        position_2_x.append(float(stringvalues[3]))
    if i == 0 and j == 2:
        position_3_y.append(float(stringvalues[2]))
        position_3_x.append(float(stringvalues[3]))

    count+= 1


position_1_y = np.array(position_1_y)
position_2_y = np.array(position_2_y)
position_3_y = np.array(position_3_y)

position_1_x = np.array(position_1_x)
position_2_x = np.array(position_2_x)
position_3_x = np.array(position_3_x)



# delta_bot = (position_2_x - position_1_x - 0.05E-04) * 1000
# delta_bot_y = (position_2_y - position_1_y) * 100000
# delta_right_y = (position_3_y - position_1_y)
# delta_right_x = (position_3_x - position_1_x)
# delta_left_y = (position_2_y - position_3_y)
# delta_left_x = (position_2_x - position_3_x)


# delta_right = []
# delta_left = []
# for i in range(len(delta_left_y)):
#     delta_right.append((math.sqrt(delta_right_x[i]**2 + delta_right_y[i]**2) - 0.05E-04) * 1000)
#     delta_left.append((math.sqrt(delta_left_x[i]**2 + delta_left_y[i]**2) - 0.05E-04) * 1000)

# delta_right = np.array(delta_right)
# delta_left = np.array(delta_left)

plt.figure(1)
plt.plot(position_1_y)
plt.plot(position_1_x)
plt.plot(position_2_y)
plt.plot(position_2_x)
plt.plot(position_3_y)
plt.plot(position_3_x)

plt.legend(["bot sigma","bot tau","right sigma","right tau","left sigma","left tau"])
plt.show()

