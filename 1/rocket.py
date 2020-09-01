"Egen kode"

from ast2000tools.space_mission import SpaceMission
from ast2000tools.constants import G
import matplotlib.pyplot as plt
import numpy as np

seed = 98638
mission = SpaceMission(seed)

with open("box_data.dat", "r") as datafile:
	box_fuel_consumption, box_propulsion = [float(data) for data in datafile.readline().split(",")]

print(f"bfc: {box_fuel_consumption} kg/s, bpr: {box_propulsion} N")

def get_consumed_fuel(thrust_force, fuel_consumption, initial_rocket_mass, delta_v, dt=1e-2):
	# get consumed fuel using numericly using eulers method
	rocket_mass = initial_rocket_mass
	vel = 0
	while vel < delta_v:
		acc = thrust_force/rocket_mass
		vel += acc*dt
		rocket_mass -= fuel_consumption*dt
		if rocket_mass < mission.spacecraft_mass:
			raise ValueError(f"rocket ran out of fuel at {100*vel/delta_v:.3g} %")
	return initial_rocket_mass - rocket_mass

def get_gravity(from_surface, mass):
	planet_mass = mission.system.masses[0]*1.98847e30  # mass  in  kg
	planet_radius = mission.system.radii[0]*1000       # radius in m
	return G*(mass*planet_mass)/(planet_radius + from_surface)**2

rocket_box_count = 1.5e12
fuel_consumption = box_fuel_consumption*rocket_box_count
propulsion = box_propulsion*rocket_box_count

extra_fuel = mission.spacecraft_mass

print(f"\nfc: {fuel_consumption} kg/s, pr: {propulsion} N")
print(f"gravity: {get_gravity(0, mission.spacecraft_mass + extra_fuel)}")

print(f"\n0 - 1000 m/s: {get_consumed_fuel(propulsion, fuel_consumption, mission.spacecraft_mass+extra_fuel, 1000, 0.01):.4g} kg fuel")

accs = []
vels = []
poss = []
tims = []
v_0 = np.array([611, 0, 0])
def launch_rocket(thrust_force, fuel_consumption, rocket_mass, fuel_mass, dt=1e-4):
	vel = v_0
	pos = np.array([0, 0, 0])
	total_mass = rocket_mass + fuel_mass
	total_time = 0

	while np.linalg.norm(vel) < 15.93*1000:
		force = thrust_force - get_gravity(pos[2], total_mass)
		acc = force/rocket_mass
		vel += np.array([0, 0, acc*dt])
		pos += vel*dt
		total_mass -= fuel_consumption*dt
		if total_mass < rocket_mass:
			raise ValueError(f"rocket ran out of fuel at r={pos:.3g}")
		total_time += dt

	return pos, rocket_mass + fuel_mass - total_mass, total_time

end_pos, consumed_mass, elapsed_time = launch_rocket(propulsion, fuel_consumption, mission.spacecraft_mass, extra_fuel, 1e-3)
print(f"\n\n0 - 1000 m/s: {end_pos} m, {consumed_mass} kg, {elapsed_time/60:.3g} min")

# plt.plot(tims, accs)
# plt.plot(tims, vels)
# plt.plot(tims, poss)
# plt.savefig("temp.png", dpi=100)
