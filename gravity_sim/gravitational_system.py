"""
Class to represent a system of celestial bodies
"""

import csv
import math
from pathlib import Path
import plotting
import typing
import warnings

import numpy as np

import common


class GravitationalSystem:
    """
    Initial state of a gravitational system
    """

    # Conversion factor from km^3 s^-2 to AU^3 d^-2
    CONVERSION_FACTOR = (86400**2) / (149597870.7**3)
    # GM values (km^3 s^-2)
    # ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
    GM_SI = {
        "Sun": 132712440041.279419,
        "Mercury": 22031.868551,
        "Venus": 324858.592000,
        "Earth": 398600.435507,
        "Mars": 42828.375816,
        "Jupiter": 126712764.100000,
        "Saturn": 37940584.841800,
        "Uranus": 5794556.400000,
        "Neptune": 6836527.100580,
        "Moon": 4902.800118,
        "Pluto": 975.500000,
        "Ceres": 62.62890,
        "Vesta": 17.288245,
    }
    # GM values (AU^3 d^-2)
    GM = {
        "Sun": 132712440041.279419 * CONVERSION_FACTOR,
        "Mercury": 22031.868551 * CONVERSION_FACTOR,
        "Venus": 324858.592000 * CONVERSION_FACTOR,
        "Earth": 398600.435507 * CONVERSION_FACTOR,
        "Mars": 42828.375816 * CONVERSION_FACTOR,
        "Jupiter": 126712764.100000 * CONVERSION_FACTOR,
        "Saturn": 37940584.841800 * CONVERSION_FACTOR,
        "Uranus": 5794556.400000 * CONVERSION_FACTOR,
        "Neptune": 6836527.100580 * CONVERSION_FACTOR,
        "Moon": 4902.800118 * CONVERSION_FACTOR,
        "Pluto": 975.500000 * CONVERSION_FACTOR,
        "Ceres": 62.62890 * CONVERSION_FACTOR,
        "Vesta": 17.288245 * CONVERSION_FACTOR,
    }
    # Solar system masses (M_sun^-1)
    SOLAR_SYSTEM_MASSES = {
        "Sun": 1.0,
        "Mercury": GM_SI["Mercury"] / GM_SI["Sun"],
        "Venus": GM_SI["Venus"] / GM_SI["Sun"],
        "Earth": GM_SI["Earth"] / GM_SI["Sun"],
        "Mars": GM_SI["Mars"] / GM_SI["Sun"],
        "Jupiter": GM_SI["Jupiter"] / GM_SI["Sun"],
        "Saturn": GM_SI["Saturn"] / GM_SI["Sun"],
        "Uranus": GM_SI["Uranus"] / GM_SI["Sun"],
        "Neptune": GM_SI["Neptune"] / GM_SI["Sun"],
        "Moon": GM_SI["Moon"] / GM_SI["Sun"],
        "Pluto": GM_SI["Pluto"] / GM_SI["Sun"],
        "Ceres": GM_SI["Ceres"] / GM_SI["Sun"],
        "Vesta": GM_SI["Vesta"] / GM_SI["Sun"],
    }
    # Gravitational constant (kg^-1 m^3 s^-2):
    CONSTANT_G_SI = 6.67430e-11
    # Gravitational constant (M_sun^-1 AU^3 d^-2):
    CONSTANT_G = GM["Sun"]

    # Solar system position and velocities data
    # Units: AU-D
    # Coordinate center: Solar System Barycenter
    # Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
    # Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
    SOLAR_SYSTEM_POS = {
        "Sun": [-7.967955691533730e-03, -2.906227441573178e-03, 2.103054301547123e-04],
        "Mercury": [
            -2.825983269538632e-01,
            1.974559795958082e-01,
            4.177433558063677e-02,
        ],
        "Venus": [
            -7.232103701666379e-01,
            -7.948302026312400e-02,
            4.042871428174315e-02,
        ],
        "Earth": [-1.738192017257054e-01, 9.663245550235138e-01, 1.553901854897183e-04],
        "Mars": [-3.013262392582653e-01, -1.454029331393295e00, -2.300531433991428e-02],
        "Jupiter": [3.485202469657674e00, 3.552136904413157e00, -9.271035442798399e-02],
        "Saturn": [8.988104223143450e00, -3.719064854634689e00, -2.931937777323593e-01],
        "Uranus": [1.226302417897505e01, 1.529738792480545e01, -1.020549026883563e-01],
        "Neptune": [
            2.983501460984741e01,
            -1.793812957956852e00,
            -6.506401132254588e-01,
        ],
        "Moon": [-1.762788124769829e-01, 9.674377513177153e-01, 3.236901585768862e-04],
        "Pluto": [1.720200478843485e01, -3.034155683573043e01, -1.729127607100611e00],
        "Ceres": [-1.103880510367569e00, -2.533340440444230e00, 1.220283937721780e-01],
        "Vesta": [-8.092549658731499e-02, 2.558381434460076e00, -6.695836142398572e-02],
    }
    SOLAR_SYSTEM_VEL = {
        "Sun": [4.875094764261564e-06, -7.057133213976680e-06, -4.573453713094512e-08],
        "Mercury": [
            -2.232165900189702e-02,
            -2.157207103176252e-02,
            2.855193410495743e-04,
        ],
        "Venus": [
            2.034068201002341e-03,
            -2.020828626592994e-02,
            -3.945639843855159e-04,
        ],
        "Earth": [
            -1.723001232538228e-02,
            -2.967721342618870e-03,
            6.382125383116755e-07,
        ],
        "Mars": [1.424832259345280e-02, -1.579236181580905e-03, -3.823722796161561e-04],
        "Jupiter": [
            -5.470970658852281e-03,
            5.642487338479145e-03,
            9.896190602066252e-05,
        ],
        "Saturn": [
            1.822013845554067e-03,
            5.143470425888054e-03,
            -1.617235904887937e-04,
        ],
        "Uranus": [
            -3.097615358317413e-03,
            2.276781932345769e-03,
            4.860433222241686e-05,
        ],
        "Neptune": [
            1.676536611817232e-04,
            3.152098732861913e-03,
            -6.877501095688201e-05,
        ],
        "Moon": [
            -1.746667306153906e-02,
            -3.473438277358121e-03,
            -3.359028758606074e-05,
        ],
        "Pluto": [2.802810313667557e-03, 8.492056438614633e-04, -9.060790113327894e-04],
        "Ceres": [
            8.978653480111301e-03,
            -4.873256528198994e-03,
            -1.807162046049230e-03,
        ],
        "Vesta": [
            -1.017876585480054e-02,
            -5.452367109338154e-04,
            1.255870551153315e-03,
        ],
    }

    DEFAULT_SYSTEMS = [
        "circular_binary_orbit",
        "eccentric_binary_orbit",
        "3d_helix",
        "sun_earth_moon",
        "figure-8",
        "pyth-3-body",
        "solar_system",
        "solar_system_plus",
        "custom",
    ]

    def __init__(self, name: str = None) -> None:
        self.name = name
        self.objects_names = []
        self.x = None
        self.v = None
        self.m = None
        self.objects_count = 0
        self.G = self.CONSTANT_G

    def __eq__(self, other: str) -> bool:
        return self.name == other

    def __str__(self) -> str:
        return self.name

    def add(
        self,
        x: typing.Union[list, np.ndarray],
        v: typing.Union[list, np.ndarray],
        m: float,
        object_name: str = None,
    ) -> None:
        """
        Add a celestial body to the system

        Parameters
        ----------
        x : list or np.ndarray
            Position vector
        v : list or np.ndarray
            Velocity vector
        m : float
            Mass
        """
        if self.x is None:
            self.x = np.array(x)
            self.v = np.array(v)
            self.m = np.array([m,])
        else:
            self.x = np.vstack((self.x, x))
            self.v = np.vstack((self.v, v))
            self.m = np.append(self.m, m)

        self.objects_count = self.m.size
        self.objects_names.append(object_name)

    def add_keplerian(
        self,
        semi_major_axis: float,
        eccentricity: float,
        inclination: float,
        argument_of_periapsis: float,
        longitude_of_ascending_node: float,
        true_anomaly: float,
        m: float,
        primary_object_name: str=None,
        primary_object_x: np.ndarray=None,
        primary_object_v: np.ndarray=None,
        primary_object_m: float=None,
        object_name: str = None,
    ):
        """
        Add a celestial body to the system using Keplerian elements

        Warning: This method use the G value from the system. Make sure
                 to set the correct G value before using this method.

        Parameters
        ----------
        semi_major_axis : float
            Semi-major axis
        eccentricity : float
            Eccentricity
        inclination : float
            Inclination
        argument_of_periapsis : float
            Argument of periapsis
        longitude_of_ascending_node : float
            Longitude of ascending node
        true_anomaly : float
            True anomaly
        m : float
            Mass
        primary_object : float
            Name of the primary object
        primary_x : np.ndarray
            Position vector of the primary object
        primary_v : np.ndarray
            Velocity vector of the primary object
        primary_m : float
            Mass of the primary object
        object_name : str
            Name of the new object
        """
        if primary_object_name is not None:
            primary_object_idx = self.objects_names.index(primary_object_name)

            if primary_object_x is None:
                primary_object_x = self.x[primary_object_idx]

            if primary_object_v is None:
                primary_object_v = self.v[primary_object_idx]

            if primary_object_m is None:
                primary_object_m = self.m[primary_object_idx]

        x, v = common.keplerian_to_cartesian(
            semi_major_axis=semi_major_axis,
            eccentricity=eccentricity,
            inclination=inclination,
            argument_of_periapsis=argument_of_periapsis,
            longitude_of_ascending_node=longitude_of_ascending_node,
            true_anomaly=true_anomaly,
            total_mass=(primary_object_m + m),
            G=self.G,
        )
        self.add(x + primary_object_x, v + primary_object_v, m, object_name)

    def remove(
        self,
        name: str = None,
        index: int = None,
        indices: typing.Union[list, np.ndarray] = None,
    ) -> None:
        """
        Remove a celestial body from the system,
        either by index or by name. An array of
        indices can also be used, in which case
        the index and name will be ignored.

        Parameters
        ----------
        name : str (optional)
            Name of the celestial body to be removed
        index : int (optional)
            Index of the celestial body to be removed
        indices : list, np.ndarray (optional)
            Indices of the celestial bodies to be removed

        Raises
        ------
        ValueError
            If no parameter or more than one parameter is provided
        ValueError
            If name is not found in system
        ValueError
            If index is out of range
        """
        parameters_count = 0
        if name is not None:
            parameters_count += 1
        if index is not None:
            parameters_count += 1
        if indices is not None:
            parameters_count += 1

        if parameters_count == 0:
            raise ValueError("Error: at least one parameter must be provided.")
        if parameters_count > 1:
            raise ValueError("Error: only one parameter is allowed.")

        if indices is not None:
            if isinstance(indices, list):
                indices = np.array(indices)
            mask = np.ones(self.objects_count, dtype=bool)
            mask[indices] = False
            self.x = self.x[mask]
            self.v = self.v[mask]
            self.m = self.m[mask]
            self.objects_count = self.m.size
            self.objects_names = np.array(self.objects_names)[mask].tolist()

        else:
            if name is not None:
                if name not in self.objects_names:
                    raise ValueError("Error: name not found in system.")
                else:
                    index = self.objects_names.index(name)

            if index >= self.m.size or index < 0:
                raise ValueError("Error: index out of range.")

            if index is not None:
                self.x = np.delete(self.x, index, axis=0)
                self.v = np.delete(self.v, index, axis=0)
                self.m = np.delete(self.m, index)
                self.objects_count = self.m.size
                del self.objects_names[index]

    def center_of_mass_correction(self) -> None:
        """
        Set center of mass of position and V_CM to zero
        """
        r_cm = np.sum(self.m[:, np.newaxis] * self.x, axis=0) / np.sum(self.m)
        v_cm = np.sum(self.m[:, np.newaxis] * self.v, axis=0) / np.sum(self.m)

        self.x -= r_cm
        self.v -= v_cm

    def save(
        self,
        system_name: str = None,
        file_path: str = None,
    ) -> None:
        warnings.warn("Warning: currently the objects names are not saved.", UserWarning)
        if system_name is not None:
            self.name = system_name

        if file_path is None:
            file_path = Path(__file__).parent / "customized_systems.csv"

        with open(file_path, "a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                [self.name, self.G, self.objects_count]
                + self.m.tolist()
                + self.x.flatten().tolist()
                + self.v.flatten().tolist()
            )

        print(f"System \"{self.name}\" saved to \"{file_path}\"")

    def load(self, system: str) -> None:
        """
        Load a pre-defined system

        Parameters
        ----------
        system : str
            Name of the system
        G : float, optional
            Gravitational constant, by default CONSTANT_G

        Raises
        ------
        FileNotFoundError
            If system is not recognized and customized_systems.csv is not found
        ValueError
            If system is not recognized in both pre-defined system and
            customized_systems.csv
        """

        loaded_system_flag = False
        match system:
            # Default systems
            case "circular_binary_orbit":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                R1 = np.array([1.0, 0.0, 0.0])
                R2 = np.array([-1.0, 0.0, 0.0])
                V1 = np.array([0.0, 0.5, 0.0])
                V2 = np.array([0.0, -0.5, 0.0])
                self.add(R1, V1, 1.0 / self.G, None)
                self.add(R2, V2, 1.0 / self.G, None)

            case "eccentric_binary_orbit":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                R1 = np.array([1.0, 0.0, 0.0])
                R2 = np.array([-1.25, 0.0, 0.0])
                V1 = np.array([0.0, 0.5, 0.0])
                V2 = np.array([0.0, -0.625, 0.0])
                self.add(R1, V1, 1.0 / self.G, None)
                self.add(R2, V2, 0.8 / self.G, None)

            case "3d_helix":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                R1 = np.array([0.0, 0.0, -1.0])
                R2 = np.array([-math.sqrt(3.0) / 2.0, 0.0, 0.5])
                R3 = np.array([math.sqrt(3.0) / 2.0, 0.0, 0.5])
                v0 = math.sqrt(1.0 / math.sqrt(3))
                V1 = np.array([-v0, 0.5, 0.0])
                V2 = np.array([0.5 * v0, 0.5, (math.sqrt(3.0) / 2.0) * v0])
                V3 = np.array([0.5 * v0, 0.5, -(math.sqrt(3.0) / 2.0) * v0])
                self.add(R1, V1, 1.0 / self.G, None)
                self.add(R2, V2, 1.0 / self.G, None)
                self.add(R3, V3, 1.0 / self.G, None)

            case "sun_earth_moon":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                m = np.array(
                    [
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Sun"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Earth"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Moon"],
                    ]
                )

                R1 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Sun"])
                R2 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Earth"])
                R3 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Moon"])
                V1 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Sun"])
                V2 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Earth"])
                V3 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Moon"])
                self.add(R1, V1, m[0], "Sun")
                self.add(R2, V2, m[1], "Earth")
                self.add(R3, V3, m[2], "Moon")

                self.center_of_mass_correction()

            case "figure-8":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                R1 = np.array([0.970043, -0.24308753, 0.0])
                R2 = np.array([-0.970043, 0.24308753, 0.0])
                R3 = np.array([0.0, 0.0, 0.0])
                V1 = np.array([0.466203685, 0.43236573, 0.0])
                V2 = np.array([0.466203685, 0.43236573, 0.0])
                V3 = np.array([-0.93240737, -0.86473146, 0.0])
                x = np.array([R1, R2, R3])
                v = np.array([V1, V2, V3])
                m = np.array([1.0 / self.G, 1.0 / self.G, 1.0 / self.G])
                self.add(R1, V1, 1.0 / self.G)
                self.add(R2, V2, 1.0 / self.G)
                self.add(R3, V3, 1.0 / self.G)

            case "pyth-3-body":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                R1 = np.array([1.0, 3.0, 0.0])
                R2 = np.array([-2.0, -1.0, 0.0])
                R3 = np.array([1.0, -1.0, 0.0])
                V1 = np.array([0.0, 0.0, 0.0])
                V2 = np.array([0.0, 0.0, 0.0])
                V3 = np.array([0.0, 0.0, 0.0])
                self.add(R1, V1, 3.0 / self.G)
                self.add(R2, V2, 4.0 / self.G)
                self.add(R3, V3, 5.0 / self.G)

            case "solar_system":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                m = np.array(
                    [
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Sun"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Mercury"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Venus"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Earth"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Mars"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Jupiter"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Saturn"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Uranus"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Neptune"],
                    ]
                )

                R1 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Sun"])
                R2 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Mercury"])
                R3 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Venus"])
                R4 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Earth"])
                R5 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Mars"])
                R6 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Jupiter"])
                R7 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Saturn"])
                R8 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Uranus"])
                R9 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Neptune"])

                V1 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Sun"])
                V2 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Mercury"])
                V3 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Venus"])
                V4 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Earth"])
                V5 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Mars"])
                V6 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Jupiter"])
                V7 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Saturn"])
                V8 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Uranus"])
                V9 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Neptune"])

                self.add(R1, V1, m[0], "Sun")
                self.add(R2, V2, m[1], "Mercury")
                self.add(R3, V3, m[2], "Venus")
                self.add(R4, V4, m[3], "Earth")
                self.add(R5, V5, m[4], "Mars")
                self.add(R6, V6, m[5], "Jupiter")
                self.add(R7, V7, m[6], "Saturn")
                self.add(R8, V8, m[7], "Uranus")
                self.add(R9, V9, m[8], "Neptune")

                self.center_of_mass_correction()

            case "solar_system_plus":
                loaded_system_flag = True
                self.G = self.CONSTANT_G
                m = np.array(
                    [
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Sun"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Mercury"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Venus"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Earth"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Mars"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Jupiter"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Saturn"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Uranus"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Neptune"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Pluto"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Ceres"],
                        GravitationalSystem.SOLAR_SYSTEM_MASSES["Vesta"],
                    ]
                )

                R1 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Sun"])
                R2 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Mercury"])
                R3 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Venus"])
                R4 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Earth"])
                R5 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Mars"])
                R6 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Jupiter"])
                R7 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Saturn"])
                R8 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Uranus"])
                R9 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Neptune"])
                R10 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Pluto"])
                R11 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Ceres"])
                R12 = np.array(GravitationalSystem.SOLAR_SYSTEM_POS["Vesta"])

                V1 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Sun"])
                V2 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Mercury"])
                V3 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Venus"])
                V4 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Earth"])
                V5 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Mars"])
                V6 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Jupiter"])
                V7 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Saturn"])
                V8 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Uranus"])
                V9 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Neptune"])
                V10 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Pluto"])
                V11 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Ceres"])
                V12 = np.array(GravitationalSystem.SOLAR_SYSTEM_VEL["Vesta"])

                self.add(R1, V1, m[0], "Sun")
                self.add(R2, V2, m[1], "Mercury")
                self.add(R3, V3, m[2], "Venus")
                self.add(R4, V4, m[3], "Earth")
                self.add(R5, V5, m[4], "Mars")
                self.add(R6, V6, m[5], "Jupiter")
                self.add(R7, V7, m[6], "Saturn")
                self.add(R8, V8, m[7], "Uranus")
                self.add(R9, V9, m[8], "Neptune")
                self.add(R10, V10, m[9], "Pluto")
                self.add(R11, V11, m[10], "Ceres")
                self.add(R12, V12, m[11], "Vesta")

                self.center_of_mass_correction()

            # Customized system
            case _:
                file_path = Path(__file__).parent / "customized_systems.csv"
                try:
                    with open(file_path, "r") as file:
                        reader = csv.reader(file)
                        for row in reader:
                            if system == row[0]:
                                loaded_system_flag = True
                                self.G = float(row[1])
                                objects_count = int(row[2])
                                m = np.zeros(objects_count)
                                for i in range(objects_count):
                                    m[i] = row[3 + i]

                                x = np.zeros(3)
                                v = np.zeros(3)
                                for i in range(objects_count):
                                    for j in range(3):
                                        x[j] = row[3 + objects_count + i * 3 + j]
                                        v[j] = row[
                                            3
                                            + objects_count
                                            + objects_count * 3
                                            + i * 3
                                            + j
                                        ]
                                    self.add(x, v, m[i])

                except FileNotFoundError:
                    warnings.warn(
                        "customized_systems.csv not found in the "
                        + f"{Path(__file__).parent} folder."
                    )

        if loaded_system_flag:
            self.name = system
        else:
            raise ValueError("Error: system name is not recognized.")

    def plot_2d_system(
        self,
        colors=None,
        labels=None,
        legend=False,
        xlabel="$x$ (AU)",
        ylabel="$y$ (AU)",
        marker="o",
        markersize=6,
    ) -> None:
        initial_state = np.concatenate(
            [
                self.x.flatten(),
                self.v.flatten(),
            ]
        )[np.newaxis, :]

        plotting.plot_2d_trajectory(
            self.objects_count,
            initial_state,
            colors,
            labels,
            legend,
            xlabel,
            ylabel,
            marker,
            markersize,
        )

    def plot_3d_system(
        self,
        colors=None,
        labels=None,
        legend=False,
        xlabel="$x$ (AU)",
        ylabel="$y$ (AU)",
        zlabel="$z$ (AU)",
        marker="o",
        markersize=6,
    ) -> None:
        initial_state = np.concatenate(
            [
                self.x.flatten(),
                self.v.flatten(),
            ]
        )[np.newaxis, :]

        plotting.plot_3d_trajectory(
            self.objects_count,
            initial_state,
            colors,
            labels,
            legend,
            xlabel,
            ylabel,
            zlabel,
            marker,
            markersize,
        )

    def sort_by_distance(
        self, primary_object_name: str = None, primary_object_index: int = None
    ):
        """
        Sort objects by distance from the origin
        """
        if primary_object_name is not None and primary_object_index is not None:
            warnings.warn(
                "Both primary_object_name and primary_object_idx are provided. "
                + "primary_object_idx will be used."
            )
        elif primary_object_name is not None:
            primary_object_index = self.objects_names.index(primary_object_name)

        x_diff = self.x - self.x[primary_object_index]
        distances = np.linalg.norm(x_diff, axis=1)
        sorted_indices = np.argsort(distances)
        self.x = self.x[sorted_indices]
        self.v = self.v[sorted_indices]
        self.m = self.m[sorted_indices]
        self.objects_names = [self.objects_names[i] for i in sorted_indices]
