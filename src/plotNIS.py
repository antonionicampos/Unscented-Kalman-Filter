import matplotlib.pyplot as plt

nisfile = open("NIS.txt", "r")
NIS_radar = []
NIS_laser = []

for line in nisfile:
	NIS = line.split(" ");
	NIS_radar.append(NIS[0])
	NIS_laser.append(NIS[1][:-1])

time = range(len(NIS_radar))

f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

ax1.plot(time, NIS_radar)
ax1.axhline(y=7.815, color='r', linestyle='-')
ax1.set_title('NIS for Radar')
ax1.set_ylabel('NIS')

ax2.plot(time, NIS_laser)
ax2.axhline(y=5.991, color='r', linestyle='-')
ax2.set_title('NIS for Lidar')
ax2.set_xlabel('Time Step')
ax2.set_ylabel('NIS')

plt.show()
