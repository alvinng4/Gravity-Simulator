import sys

import pygame
import pygame.gfxdraw
from pygame.sprite import Sprite

from textbox import Textbox


class Grav_obj(Sprite):
    def __init__(
        self, grav_sim, pos_x, pos_y, real_scale, name: str = None, img_path: str = None
    ):
        super().__init__()
        self.screen = grav_sim.screen
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.color = self.settings.grav_obj_color
        if not name.strip().title() == "Sun":
            self.img_scale = grav_sim.camera.img_scale
        else:
            self.img_scale = 1

        self.pos_x = pos_x
        self.pos_y = pos_y
        self.diameter = real_scale * self.img_scale * self.settings.screen_height
        if img_path:
            try:
                load_image = pygame.image.load(img_path).convert_alpha()
                self.image = pygame.transform.scale(
                    load_image, (self.diameter, self.diameter)
                )
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Check that you are in the correct directory, and the image path provided for Grav_obj is correct."
                )

        self.rect.centerx = pos_x
        self.rect.centery = pos_y

    def update(self):
        """Update the apparent position of all grav_objs with camera"""
        self.rect.center = (
            self.pos_x - self.camera.pos_x,
            self.pos_y - self.camera.pos_y,
        )
