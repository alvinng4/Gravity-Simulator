import sys

import pygame
import pygame.gfxdraw
from pygame.sprite import Sprite
import numpy


class Grav_obj(Sprite):
    def __init__(self, grav_sim, x0, y0, img_path: str = None):
        super().__init__()
        self.screen = grav_sim.screen
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.color = self.settings.GRAV_OBJ_COLOR

        if img_path:
            try:
                self.image = pygame.image.load(img_path).convert_alpha()
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Check that you are in the correct directory, and the image path provided for Grav_obj is correct."
                )

        self.rect.centerx = x0
        self.rect.centery = y0

    def draw(self):
        """Draw the object at its current location."""
        draw_rect = self.rect.copy()
        draw_rect.center = (
            draw_rect.centerx - self.camera.pos_x,
            draw_rect.centery - self.camera.pos_y,
        )
        self.screen.blit(self.image, draw_rect)


class solar_system(Grav_obj):
    def __init__(self, grav_sim):
        """Initialize the solar system"""
        super().__init__()
        self.sun = Grav_obj(grav_sim, "images/sun.png")
        self.mercury = Grav_obj(grav_sim, "images/mercury.png")
        self.venus = Grav_obj(grav_sim, "images/venus.png")
        self.earth = Grav_obj(grav_sim, "images/earth.png")
        self.mars = Grav_obj(grav_sim, "images/mars.png")
        self.jupiter = Grav_obj(grav_sim, "images/jupiter.png")
        self.saturn = Grav_obj(grav_sim, "images/saturn.png")
        self.uranus = Grav_obj(grav_sim, "images/uranus.png")
        self.neptune = Grav_obj(grav_sim, "images/neptune.png")
