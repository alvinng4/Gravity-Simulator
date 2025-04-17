from pathlib import Path

import h5py
import numpy as np
from grav_sim import GravitySimulatorAPI

IC_path = Path(__file__).parent.parent / "galaxy_collision_IC.hdf5"

G_cgs = 6.67430e-8  # cm^3 g^-1 s^-2

ACC_METHOD = "Barnes_Hut"
OPENING_ANGLE = 0.5
SOFTENING_LENGTH = 0.0

TF = 4000000000 * 365.24 * 24 * 3600 # 4 billion years

INTEGRATOR = "Leapfrog"
DT = TF / 2000

OUTPUT_METHOD = "hdf5"
OUTPUT_INTERVAL = TF / 500 # 500 snapshots
OUTPUT_DIR = Path(__file__).parent / "snapshots/"

def main():
    # Load the simulation data
    with h5py.File(IC_path, "r") as file:
        type_1_particle_IDs = file["PartType1"]["ParticleIDs"][:]
        type_2_particle_IDs = file["PartType2"]["ParticleIDs"][:]

        type_1_mass = file["PartType1"]["Masses"][:]
        type_1_position = file["PartType1"]["Coordinates"][:]
        type_1_velocity = file["PartType1"]["Velocities"][:]

        type_2_mass = file["PartType2"]["Masses"][:]
        type_2_position = file["PartType2"]["Coordinates"][:]
        type_2_velocity = file["PartType2"]["Velocities"][:]

        unit_length_in_cgs = file["Units"].attrs["Unit length in cgs (U_L)"]
        unit_mass_in_cgs = file["Units"].attrs["Unit mass in cgs (U_M)"]
        unit_time_in_cgs = file["Units"].attrs["Unit time in cgs (U_t)"]

    # Combine the data into a single array
    particle_IDs = np.concatenate((type_1_particle_IDs, type_2_particle_IDs))
    masses = np.concatenate((type_1_mass, type_2_mass))
    positions = np.concatenate((type_1_position, type_2_position))
    velocities = np.concatenate((type_1_velocity, type_2_velocity))

    # Create the simulator
    gs = GravitySimulatorAPI()

    system = gs.get_new_system()
    system.add(
        x = positions,
        v = velocities,
        m = masses,
    )

    system.G = G_cgs * unit_mass_in_cgs * unit_time_in_cgs**2 / unit_length_in_cgs**3
    # Convert the units
    global TF, DT, OUTPUT_INTERVAL
    TF /= unit_time_in_cgs
    DT /= unit_time_in_cgs
    OUTPUT_INTERVAL /= unit_time_in_cgs

    # Parameters
    acc_param, integrator_param, output_param, settings = gs.get_new_parameters()

    acc_param.method = ACC_METHOD
    acc_param.opening_angle = OPENING_ANGLE
    acc_param.softening_length = SOFTENING_LENGTH

    integrator_param.integrator = INTEGRATOR
    integrator_param.dt = DT

    output_param.method = OUTPUT_METHOD
    output_param.output_interval = OUTPUT_INTERVAL
    output_param.output_dir = OUTPUT_DIR
    output_param.coordinate_output_dtype = "float"
    output_param.velocity_output_dtype = "float"
    output_param.mass_output_dtype = "float"

    gs.launch_simulation(system, acc_param, integrator_param, output_param, settings, TF)

if __name__ == "__main__":
    main()
