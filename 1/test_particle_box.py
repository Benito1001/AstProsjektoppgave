"egen kode"

import time
import numpy as np
from particle_box import EngineBox
import ast2000tools.constants as cn
import matplotlib.pyplot as plt

box = EngineBox(100000, cn.m_H2, 1e-6, 3e3)
change_box = box # EngineBox(100000, cn.m_H2, 1e-6, 3e3, change_velocities=True)


timesteps = 1000
dt = (1e-6/1e4)/20
total_time = dt*timesteps


"plot histogram"
# normal = lambda x, sigma, mu: (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*((x - mu)/sigma)**2)
# x_ray = np.linspace(-4*box.stande, 4*box.stande, 1000)
# normal_ray = normal(x_ray, box.stande, 0)
# plt.figure(figsize=(14, 8))
# plt.hist(box.vel.x, 100, density=True, color="#00dd0044", label="$v_x$ histogram")
# plt.hist(box.vel.y, 100, density=True, color="#dd000044", label="$v_y$ histogram")
# plt.hist(box.vel.z, 100, density=True, color="#0000dd44", label="$v_z$ histogram")
# plt.plot(x_ray, normal_ray, label="expected normal distribution")
# plt.legend()
# plt.savefig("box_histogram.png", dpi=200)


"run simulation and save values"
posx0 = np.zeros(timesteps)
posy0 = np.zeros(timesteps)
posz0 = np.zeros(timesteps)

times = np.zeros(timesteps)
vel_change_temps = np.zeros(timesteps)
vel_same_temps = np.zeros(timesteps)

pressure = np.zeros(timesteps)

incy = timesteps/100
start_time = time.time()
for i in range(timesteps):
	if i % incy == 0:
		print(f"\r{i/incy:.0f}/{timesteps/incy:.0f} %", end="")
	posx0[i] = box.pos.x[0]
	posy0[i] = box.pos.y[0]
	posz0[i] = box.pos.z[0]

	times[i] = i*dt
	vel_same_temps[i] = (box.atom_mass*np.sum(box.vel.x**2 + box.vel.y**2 + box.vel.z**2)/box.atom_count)/(3*cn.k_B)
	vel_change_temps[i] = (change_box.atom_mass*np.sum(change_box.vel.x**2 + change_box.vel.y**2 + change_box.vel.z**2)/change_box.atom_count)/(3*cn.k_B)

	pressure[i] = box.atom_count*cn.k_B*vel_same_temps[i]

	box.step(dt)
	# change_box.step(dt)
print(f"\ntime: {time.time() - start_time} s")

"""plot temperature"""
# fig = plt.figure(figsize=(10, 10))
# plt.plot([times[0], times[-1]], [3000, 3000], "r", label="expected temperature")
# plt.plot(times, vel_change_temps, label="actual temperature, change velocities")
# plt.plot(times, vel_same_temps, markersize=3, label="actual temperature, same velocities")
# plt.legend()
# plt.savefig("box_temp.png", dpi=200)

# """plot pressure"""
# fig = plt.figure(figsize=(10, 10))
# plt.plot([times[0], times[-1]], [box.atom_count*cn.k_B*3000, box.atom_count*cn.k_B*3000], "r", label="expected pressure")
# plt.plot(times, pressure, label="pressure in the rocket motor over time")
# plt.legend()
# plt.savefig("box_pressure.png", dpi=200)

"""plot position"""
# fig = plt.figure(figsize=(6, 6))
# plt.plot(posx0[0], posz0[0], "go", label="start")
# sub_posx0 = []
# sub_posz0 = []
# for i in range(len(posx0) - 1):
# 	x, z = posx0[i], posz0[i]
# 	distance = (x - posx0[i+1])**2 + (z - posz0[i+1])**2
# 	if distance > (box.box_size/2)**2:
# 		plt.plot(sub_posx0, sub_posz0, ".", markersize=3)
# 		plt.plot(x, z, "ko", label="escaped box")
# 		plt.plot(posx0[i+1], posz0[i+1], "co", label="new particle")
# 		sub_posx0 = []
# 		sub_posz0 = []
# 	else:
# 		sub_posx0.append(x)
# 		sub_posz0.append(z)
#
# plt.plot(sub_posx0, sub_posz0, ".", markersize=3)
# plt.plot(posx0[-1], posz0[-1], "ro", label="end")
# plt.legend()
# plt.xlabel("x [m]")
# plt.ylabel("z [m]")
# plt.xlim(0, box.box_size)
# plt.ylim(0, box.box_size)
# for axticks in [plt.xticks, plt.yticks]:
# 	locs, labels = axticks()
# 	[label.set_text(loc) for loc, label in zip(locs, labels)]
# 	axticks(locs, labels)
# plt.savefig("box_position.png", dpi=200)
