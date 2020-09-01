"Egen kode"

import time
import numpy as np
import ast2000tools.constants as cn

generator = np.random.default_rng(69)

class Vec3Matrix:
	"""
	this class is basicly just an abstraction og an 3 x n numpy array to
	allow you to use .x/.y/.z to get the values
	"""

	def __init__(self, value_matrix):
		self.values = value_matrix

	@property
	def x(self):
		return self.values[0]

	@property
	def y(self):
		return self.values[1]

	@property
	def z(self):
		return self.values[2]

	def __add__(self, other):
		return Vec3Matrix(self.values + other.values)

	def __mul__(self, scalar):
		return Vec3Matrix(scalar*self.values)

class EngineBox:
	def __init__(self, atom_count, atom_mass, box_size, box_temperature, change_velocities=False):
		self.change_velocities = change_velocities

		self.atom_count = atom_count
		self.atom_mass = atom_mass
		self.box_size = box_size
		self.stande = np.sqrt((cn.k_B*box_temperature)/atom_mass)
		self.pos, self.vel = self.generate_atoms()

		self.density = atom_count/box_size**3

		self.total_atom_loss = 0
		self.total_momentum_loss = 0

	def generate_atoms(self):
		pos = Vec3Matrix(generator.uniform(0, self.box_size, size=(3, self.atom_count)))
		vel = Vec3Matrix(generator.normal(0, self.stande, size=(3, self.atom_count)))
		return pos, vel

	def release_atoms(self):
		"""
		the escape hole is located in the middle of the bottom of the box,
		new particles appear from a hole in the middle of the top of the box
		"""
		hole_size = self.box_size/2
		hole_left = self.box_size/2 - hole_size/2
		hole_right = self.box_size/2 + hole_size/2

		x_vals = (self.pos.x > hole_left) & (self.pos.x < hole_right)
		y_vals = (self.pos.y > hole_left) & (self.pos.y < hole_right)
		indices = (self.pos.z < 0) & x_vals & y_vals

		escaped_count = np.sum(indices)
		lost_momentum = self.atom_mass*np.sum(self.vel.z)

		# this would look bettes as self.vel.values[:, indices] = ... , but that is actualy noticeably slower
		self.pos.x[indices], self.pos.y[indices], self.pos.z[indices] = *generator.uniform(hole_left, hole_right, size=(2, escaped_count)), np.full(escaped_count, self.box_size)
		if self.change_velocities:
			# changing the velocity makes the temperature decrease over time
			self.vel.x[indices], self.vel.y[indices], self.vel.z[indices] = generator.uniform(0, self.box_size, size=(3, escaped_count))

		return escaped_count, lost_momentum

	def step(self, dt):
		self.pos += self.vel*dt

		lost_atoms, lost_momentum = self.release_atoms()
		self.total_atom_loss += lost_atoms
		self.total_momentum_loss += lost_momentum

		self.vel.values[(self.pos.values < 0) | (self.pos.values > self.box_size)] *= -1

"preformance profiling"
if __name__ == "__main__":
	box = EngineBox(100000, cn.m_H2, 1e-6, 3e3)

	timesteps = 10000
	dt = (1e-6/1e4)/20
	total_time = dt*timesteps

	incy = timesteps/100
	start_time = time.time()
	for i in range(timesteps):
		if i % incy == 0:
			print(f"\r{i/incy:.0f}/{timesteps/incy:.0f} %", end="")
		box.step(dt)
	print(f"\ntime: {time.time() - start_time} s")

	# 10000 timesteps of 100000 atoms on desktop:
	#   start:        12.77 s, 750549, -2.55e-16 Ns
	#   fancy vector: 11.60 s, 750549, -2.55e-16 Ns   (somewhat faster, somewhat prettier)
	#   thing
