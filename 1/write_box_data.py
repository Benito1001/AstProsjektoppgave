from particle_box import EngineBox
import ast2000tools.constants as cn

box = EngineBox(100000, cn.m_H2, 1e-6, 3e3)

timesteps = 1000
dt = (1e-6/1e4)/20
total_time = dt*timesteps

incy = timesteps/100
for i in range(timesteps):
	if i % incy == 0:
		print(f"\r{i/incy:.0f}/{timesteps/incy:.0f} %", end="")
	box.step(dt)
print()

with open("box_data.dat", "w") as datafile:
	datafile.write(f"{(box.total_atom_loss*cn.m_H2)/total_time}, {-box.total_momentum_loss/total_time}")
