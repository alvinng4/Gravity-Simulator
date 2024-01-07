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
        apparent_scale,
        img_path: str = None,
    ):
        super().__init__()
        self.screen = grav_sim.screen
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.params = params
        self.diameter = apparent_scale * self.settings.screen_height
        # Note: apparent_scale = real_scale * img_scale
        if img_path:
            try:
                load_image = pygame.image.load(img_path).convert_alpha()
                self.image = pygame.transform.scale(
                    load_image, (self.diameter, self.diameter)
                )
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Make sure the image path provided for Grav_obj is correct."
                )
        else:
            pass
            # Create new object by holding left click.

        self.rect.centerx = self.params["x"]
        self.rect.centery = self.params["y"]

    def update(self):
        """Update the apparent position of all grav_objs with camera"""
        self.rect.center = (
            self.params["x"] - self.camera.pos_x,
            self.params["y"] - self.camera.pos_y,
        )

    @classmethod
    def create_solor_system(self, grav_sim):
        main_dir_path = os.path.dirname(__file__)
        centerx = grav_sim.screen.get_rect().centerx
        centery = grav_sim.screen.get_rect().centery
        path_sun = os.path.join(main_dir_path, "images/sun.png")
        path_mercury = os.path.join(main_dir_path, "images/mercury.png")
        path_venus = os.path.join(main_dir_path, "images/venus.png")
        path_earth = os.path.join(main_dir_path, "images/earth.png")
        path_mars = os.path.join(main_dir_path, "images/mars.png")
        path_jupiter = os.path.join(main_dir_path, "images/jupiter.png")
        path_saturn = os.path.join(main_dir_path, "images/saturn.png")
        path_uranus = os.path.join(main_dir_path, "images/uranus.png")
        path_neptune = os.path.join(main_dir_path, "images/neptune.png")
        sun = Grav_obj(
            grav_sim,
            {"x": centerx, "y": centery},
            0.05 * grav_sim.settings.sun_img_scale,
            path_sun,
        )
        mercury = Grav_obj(
            grav_sim,
            {"x": 1.1 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_mercury,
        )
        venus = Grav_obj(
            grav_sim,
            {"x": 1.2 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_venus,
        )
        earth = Grav_obj(
            grav_sim,
            {"x": 1.3 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_earth,
        )
        mars = Grav_obj(
            grav_sim,
            {"x": 1.5 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_mars,
        )
        jupiter = Grav_obj(
            grav_sim,
            {"x": 1.8 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_jupiter,
        )
        saturn = Grav_obj(
            grav_sim,
            {"x": 2 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_saturn,
        )
        uranus = Grav_obj(
            grav_sim,
            {"x": 4 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
            path_uranus,
        )
        neptune = Grav_obj(
            grav_sim,
            {"x": 5 * centerx, "y": centery},
            0.0005 * grav_sim.settings.img_scale,
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
