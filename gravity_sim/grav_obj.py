import sys
import os

import pygame
import pygame.gfxdraw
from pygame.sprite import Sprite


class Grav_obj(Sprite):
    def __init__(
        self,
        grav_sim,
        params: dict,
        img_path: str = None,
        name: str = None,
    ):
        super().__init__()
        self.screen = grav_sim.screen
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.params = params
        self.diameter = 2 * self.params["R"]
        if name == "Sun":
            self.img_diameter = (
                    self.diameter * 0.25 * self.settings.sun_img_scale * self.settings.screen_height
                )
        else:
            self.img_diameter = (
                self.diameter * 0.25 * self.settings.img_scale * self.settings.screen_height
            )

        # Note: apparent_scale = real_scale * img_scale
        if img_path:
            try:
                load_image = pygame.image.load(img_path).convert_alpha()
                self.image = pygame.transform.scale(
                    load_image, (self.img_diameter, self.img_diameter)
                )
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Make sure the image path provided for Grav_obj is correct."
                )
        else:
            pass
            # Create new object by holding left click.

        self.rect.centerx = self.params["r1"]
        self.rect.centery = self.params["r2"]

    def update(self):
        """Update the apparent position of all grav_objs with camera"""
        self.rect.center = (
            self.params["r1"] * 0.25 * self.settings.screen_height + self.screen.get_rect().centerx - self.camera.pos_x,
            - self.params["r2"] * 0.25 * self.settings.screen_height + self.screen.get_rect().centery - self.camera.pos_y,
        )

    @classmethod
    def create_grav_obj(self, grav_sim):
        pass
        # grav_obj = Grav_obj(grav_sim, params, 1)
        # grav_sim.grav_objs.add(grav_obj)

    @classmethod
    def create_solor_system(self, grav_sim):
        main_dir_path = os.path.dirname(__file__)
        path_sun = os.path.join(main_dir_path, "images/sun.png")
        path_mercury = os.path.join(main_dir_path, "images/mercury.png")
        path_venus = os.path.join(main_dir_path, "images/venus.png")
        path_earth = os.path.join(main_dir_path, "images/earth.png")
        path_mars = os.path.join(main_dir_path, "images/mars.png")
        path_jupiter = os.path.join(main_dir_path, "images/jupiter.png")
        path_saturn = os.path.join(main_dir_path, "images/saturn.png")
        path_uranus = os.path.join(main_dir_path, "images/uranus.png")
        path_neptune = os.path.join(main_dir_path, "images/neptune.png")
        # Data dated on 2017 - Jun - 22
        # r1 - r3: Positions (AU), v1 - v3: Velocities (AU/d), m: Mass (Solar masses)
        sun = Grav_obj(
            grav_sim,  # Wrong data, fix later
            {
                "r1": 0,
                "r2": 0,
                "r3": 0,
                "v1": 0,
                "v2": 0,
                "v3": 0,
                "m": 1,
                "R": 0.004650467261,
            },
            path_sun,
            name = "Sun",
        )
        mercury = Grav_obj(
            grav_sim,
            {
                "r1": -1.202266214811173e-02,
                "r2": 3.129738317339329e-01,
                "r3": 2.637608427362013e-02,
                "v1": -3.374999571504927e-02,
                "v2": -3.223139133573000e-04,
                "v3": 3.069097338939302e-03,
                "m": 1.66051140935277e-07,
                "R": 1.63083872e-05,
            },
            path_mercury,
        )
        venus = Grav_obj(
            grav_sim,
            {
                "r1": 6.036656393784035e-01,
                "r2": -4.041148416177896e-01,
                "r3": -4.043024664738540e-02,
                "v1": 1.125511509001221e-02,
                "v2": 1.664344623604423e-02,
                "v3": -4.214323452340928e-04,
                "m": 2.44827371182131e-06,
                "R": 4.04537843e-05,
            },
            path_venus,
        )
        earth = Grav_obj(
            grav_sim,
            {
                "r1": 1.241238601620544e-02,
                "r2": -1.011219129247450e00,
                "r3": -9.819526634120970e-05,
                "v1": 1.692619974976791e-02,
                "v2": 1.014868174374918e-04,
                "v3": -6.047160617924358e-08,
                "m": 3.00329789031573e-06,
                "R": 4.25875046e-05,
            },
            path_earth,
        )
        mars = Grav_obj(
            grav_sim,
            {
                "r1": -4.926800380968982e-01,
                "r2": 1.537007440322637e00,
                "r3": 4.411949077995760e-02,
                "v1": -1.278895103122624e-02,
                "v2": -3.109871472361932e-03,
                "v3": 2.485576543075108e-04,
                "m": 3.22773848604808e-07,
                "R": 2.26574081e-05,
            },
            path_mars,
        )
        jupiter = Grav_obj(
            grav_sim,
            {
                "r1": -4.989446787630805e00,
                "r2": -2.184763688656508e00,
                "r3": 1.206579440062781e-01,
                "v1": 2.938669025252137e-03,
                "v2": -6.553859846840747e-03,
                "v3": -3.848907559519295e-05,
                "m": 0.000954532562518104,
                "R": 4.6732617e-04,
            },
            path_jupiter,
        )
        saturn = Grav_obj(
            grav_sim,
            {
                "r1": -9.681157061545530e-01,
                "r2": -1.000423489517810e01,
                "r3": 2.124749402062057e-01,
                "v1": 5.246396731312688e-03,
                "v2": -5.546463665223218e-04,
                "v3": -1.994279405146075e-04,
                "m": 0.00028579654259599,
                "R": 3.89256877e-04,
            },
            path_saturn,
        )
        uranus = Grav_obj(
            grav_sim,
            {
                "r1": 1.806198902787260e01,
                "r2": 8.416356280190394e00,
                "r3": -2.027374282309190e-01,
                "v1": -1.689894219169289e-03,
                "v2": 3.381692838015134e-03,
                "v3": 3.446267721267768e-05,
                "m": 4.3655207025844e-05,
                "R": 1.69534499e-04,
            },
            path_uranus,
        )
        neptune = Grav_obj(
            grav_sim,
            {
                "r1": 2.850592355224314e01,
                "r2": -9.173827312094703e00,
                "r3": -4.680305005676293e-01,
                "v1": 9.407596025584859e-04,
                "v2": 3.006460319698939e-03,
                "v3": -8.395033744165947e-05,
                "m": 5.1499991953912e-05,
                "R": 1.64587904e-04,
            },
            path_neptune,
        )
        grav_sim.grav_objs.add(sun)
        grav_sim.grav_objs.add(mercury)
        grav_sim.grav_objs.add(venus)
        grav_sim.grav_objs.add(earth)
        grav_sim.grav_objs.add(mars)
        grav_sim.grav_objs.add(jupiter)
        grav_sim.grav_objs.add(saturn)
        grav_sim.grav_objs.add(uranus)
        grav_sim.grav_objs.add(neptune)

    @classmethod
    def create_figure_8(self, grav_sim):
        pass
        # Create a figure 8 orbit
