from pathlib import Path
import sys
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)

import numpy as np
import matplotlib.pyplot as plt
import pytest

from gravity_sim.simulator import Simulator
from gravity_sim import simulator

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559

@pytest.mark.skip(reason="For testing manually")
def test():
    #integrator = "euler"
    #integrator = "euler_cromer"
    #integrator = "rk2"
    integrator = "rk4"
    #integrator = "leapfrog"
    #test_two_vectors(integrator)
    test_solar_system(integrator)
    #test_figure_8(integrator)



@pytest.mark.skip(reason="For testing manually")
def test_two_vectors():
    # Initialize
    R1 = np.array([1.0, 0.0, 0.0])
    R2 = np.array([-1.0, 0.0, 0.0])
    V1 = np.array([0.0, 0.5, 0.0])
    V2 = np.array([0.0, -0.5, 0.0])
    x = np.zeros((2, 3))
    v = np.zeros((2, 3))
    x[0] = R1
    x[1] = R2
    v[0] = V1
    v[1] = V2
    m = [1.0 / G, 1.0 / G]
    t0 = 0.0
    tf = 10000.0
    dt = 0.001

    # Simulation
    npts = int(np.floor((tf - t0) / dt)) + 1
    sol_time = np.linspace(t0, t0 + dt * (npts - 1), npts)
    energy = np.zeros(npts)
    for count, t in enumerate(sol_time):
        a = simulator.ode_n_body_first_order(2, x, m)
        x, v = simulator.leapfrog(2, x, v, a, m, dt)
        energy[count] = simulator.total_energy(2, x, v, m)

    # Plotting
    plt.figure()
    plt.semilogy(sol_time, np.abs((energy - energy[0]) / energy[0]))
    plt.xlabel("Time")
    plt.ylabel("|(E(t)-E0)/E0|")
    plt.show()

    plt.figure()
    plt.semilogy(sol_time, np.abs(energy))
    plt.xlabel("Time")
    plt.ylabel("E(t)")
    plt.show()


@pytest.mark.skip(reason="For testing manually")
def test_solar_system(integrator):
    # Initialize
    x = np.zeros((9, 3))
    v = np.zeros((9, 3))

    x[0] = np.array([-0.007967955691534, -0.002906227441573, 0.000210305430155])
    x[1] = np.array([-0.282598326953863, 0.197455979595808, 0.041774335580637])
    x[2] = np.array([-0.723210370166638, -0.079483020263124, 0.040428714281743])
    x[3] = np.array([-0.173819201725705, 0.966324555023514, 0.000155390185490])
    x[4] = np.array([-0.301326239258265, -1.454029331393295, -0.023005314339914])
    x[5] = np.array([3.485202469657675, 3.552136904413157, -0.092710354427984])
    x[6] = np.array([8.988104223143450, -3.719064854634689, -0.293193777732359])
    x[7] = np.array([12.263024178975050, 15.297387924805450, -0.102054902688356])
    x[8] = np.array([29.835014609847410, -1.793812957956852, -0.650640113225459])

    v[0] = np.array([0.000004875094764, -0.000007057133214, -0.000000045734537])
    v[1] = np.array([-0.022321659001897, -0.021572071031763, 0.000285519341050])
    v[2] = np.array([0.002034068201002, -0.020208286265930, -0.000394563984386])
    v[3] = np.array([-0.017230012325382, -0.002967721342619, 0.000000638212538])
    v[4] = np.array([0.014248322593453, -0.001579236181581, -0.000382372279616])
    v[5] = np.array([-0.005470970658852, 0.005642487338479, 0.000098961906021])
    v[6] = np.array([0.001822013845554, 0.005143470425888, -0.000161723590489])
    v[7] = np.array([-0.003097615358317, 0.002276781932346, 0.000048604332222])
    v[8] = np.array([0.000167653661182, 0.003152098732862, -0.000068775010957])

    m = [
        1.0,
        1.66051140935277e-07,
        2.44827371182131e-06,
        3.00329789031573e-06,
        3.22773848604808e-07,
        0.000954532562518104,
        0.00028579654259599,
        4.3655207025844e-05,
        5.1499991953912e-05,
    ]
    t0 = 0.0
    tf = 10000.0
    dt = 0.01

    # Simulation
    npts = int(np.floor((tf - t0) / dt)) + 1
    #sol_state = np.zeros ((npts ,len(x)))
    sol_time = np.linspace(t0, t0 + dt * (npts - 1), npts)
    energy = np.zeros(npts)

    
    match integrator:
        case "euler":
            for count, t in enumerate(sol_time):
                a = simulator.ode_n_body_first_order(9, x, m)
                x, v = simulator.euler(x, v, a, dt)
                #sol_state[count,:] = x
                energy[count] = simulator.total_energy(9, x, v, m)
        case "euler_cromer":
            for count, t in enumerate(sol_time):
                a = simulator.ode_n_body_first_order(9, x, m)
                x, v = simulator.euler_cromer(x, v, a, dt)
                #sol_state[count,:] = x
                energy[count] = simulator.total_energy(9, x, v, m)
        case "rk2":
            for count, t in enumerate(sol_time):
                a = simulator.ode_n_body_first_order(9, x, m)
                x, v = simulator.rk2(9, x, v, a, m, dt)
                #sol_state[count,:] = x
                energy[count] = simulator.total_energy(9, x, v, m)            
        case "rk4":
              for count, t in enumerate(sol_time):
                a = simulator.ode_n_body_first_order(9, x, m)
                x, v = simulator.rk4(9, x, v, a, m, dt)
                #sol_state[count,:] = x
                energy[count] = simulator.total_energy(9, x, v, m)   
        case "leapfrog":
            a = simulator.ode_n_body_first_order(9, x, m)
            for count, t in enumerate(sol_time):
                x, v, a = simulator.leapfrog(9, x, v, a, m, dt)
                #sol_state[count,:] = x
                energy[count] = simulator.total_energy(9, x, v, m)

    # Plotting
    plt.figure()
    plt.semilogy(sol_time, np.abs((energy - energy[0]) / energy[0]))
    plt.xlabel("Time")
    plt.ylabel("|(E(t)-E0)/E0|")
    plt.show()

    plt.figure()
    plt.semilogy(sol_time, np.abs(energy))
    plt.xlabel("Time")
    plt.ylabel("E(t)")
    plt.show()

    #fig = plt.figure()
    #ax = fig.add_subplot (111, aspect='equal ')
    #ax.plot(sol_state [:,0], sol_state [:,1],"b-")
    #ax.plot(sol_state [:,0 + 3], sol_state [:,1 + 3],"g-")
    #plt.show()


@pytest.mark.skip(reason="For testing manually")
def test_figure_8():
    # Initialize
    R1 = np.array([1.0, 0.0, 0.0])
    R2 = np.array([-1.0, 0.0, 0.0])
    V1 = np.array([0.0, 0.5, 0.0])
    V2 = np.array([0.0, -0.5, 0.0])
    x = np.zeros((2, 3))
    v = np.zeros((2, 3))
    x[0] = R1
    x[1] = R2
    v[0] = V1
    v[1] = V2
    m = [1.0 / G, 1.0 / G]
    t0 = 0.0
    tf = 1000.0
    dt = 0.001

    # Simulation
    npts = int(np.floor((tf - t0) / dt)) + 1
    sol_time = np.linspace(t0, t0 + dt * (npts - 1), npts)
    energy = np.zeros(npts)
    for count, t in enumerate(sol_time):
        a = Simulator.ode_n_body_first_order(2, x, m)
        x, v = Simulator.leapfrog(2, x, v, a, m, dt)
        energy[count] = Simulator.total_energy(2, x, v, m)

    # Plotting
    plt.figure()
    plt.semilogy(sol_time, np.abs((energy - energy[0]) / energy[0]))
    plt.xlabel("Time")
    plt.ylabel("|(E(t)-E0)/E0|")
    plt.show()

    plt.figure()
    plt.semilogy(sol_time, np.abs(energy))
    plt.xlabel("Time")
    plt.ylabel("E(t)")
    plt.show()






if __name__ == "__main__":
    test()