import numpy as np
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

seed = utils.get_seed('bendikdn')
print(f"seed: {seed}\n")

system = SolarSystem(seed)

print(f'My system has a {system.star_mass:g} solar mass star with a radius of {system.star_radius:g} kilometers\n')

for planet_id in range(system.number_of_planets):
	print(f"""
		Planet {planet_id:d} is a {system.types[planet_id]} planet with a mass of
		{system.masses[planet_id]*1.98847e30/5.972e24:.3g} earth masses and a radius of
		{system.radii[planet_id]:.5g} km rotating at {(2*np.pi)/(system.rotational_periods[planet_id]*24*3600):.3g} 1/s
		""")


mission = SpaceMission(seed)

print(f'\n\nMy spacecraft has a mass of {mission.spacecraft_mass:g} kg and a cross-sectional area of {mission.spacecraft_area:g} m^2.')
