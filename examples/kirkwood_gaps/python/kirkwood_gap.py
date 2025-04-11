"""
Demonstration on using the gravity simulator API to simulate the formation of Kirkwood gaps.
"""

import numpy as np

from grav_sim import GravitySimulatorAPI

gs = GravitySimulatorAPI()

NUM_PARTICLES = 100
DT = 180.0
TF = 5000000.0 * 365.24
OUTPUT_INTERVAL = 2500 * 365.24 # Stores 2000 snapshots
OUTPUT_METHOD = "hdf5"
OUTPUT_DIR = "../snapshots/"

def main():
    # ---------- Initialization ---------- #
    gs = GravitySimulatorAPI()
    system = gs.get_built_in_system("solar_system")

    # Remove Mercury, Venus, Earth, Uranus, and Neptune
    system.remove([1, 2, 3, 7, 8])

    rng = np.random.default_rng()
    a = rng.uniform(2.0, 3.35, size=NUM_PARTICLES)
    ecc = np.abs(rng.normal(loc=0.0, scale=0.12, size=NUM_PARTICLES))
    inc = np.abs(rng.normal(loc=0.0, scale=0.3, size=NUM_PARTICLES))
    argument_of_periapsis = rng.uniform(0, 2 * np.pi, size=NUM_PARTICLES)
    long_asc_node = rng.uniform(0, 2 * np.pi, size=NUM_PARTICLES)
    true_anomaly = rng.uniform(0, 2 * np.pi, size=NUM_PARTICLES)

    for i in range(NUM_PARTICLES):
        system.add_keplerian(
            semi_major_axis=a[i],
            eccentricity=ecc[i],
            inclination=inc[i],
            argument_of_periapsis=argument_of_periapsis[i],
            longitude_of_ascending_node=long_asc_node[i],
            true_anomaly=true_anomaly[i],
            m=0.0,
            primary_particle_id=0,
        )
    system.center_of_mass_correction()

    # ---------- Simulation ---------- #
    acc_param, integrator_param, output_param, settings = gs.get_new_parameters()

    acc_param.method = "massless"

    integrator_param.integrator = "whfast"
    integrator_param.dt = DT
    integrator_param.whfast_remove_invalid_particles = True

    output_param.method = OUTPUT_METHOD
    output_param.output_interval = OUTPUT_INTERVAL
    output_param.output_dir = OUTPUT_DIR
    output_param.coordinate_output_dtype = "float"
    output_param.velocity_output_dtype = "float"
    output_param.mass_output_dtype = "float"

    gs.launch_simulation(system, acc_param, integrator_param, output_param, settings, TF)
    
if __name__ == "__main__":
    main()
