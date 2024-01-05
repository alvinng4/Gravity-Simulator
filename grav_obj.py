import sys

import pygame
import pygame.gfxdraw
from pygame.sprite import Sprite


class Grav_obj(Sprite):
    def __init__(self, grav_sim, pos_x, pos_y, real_scale, name: str = None, img_path: str = None):
        super().__init__()
        self.screen = grav_sim.screen
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.color = self.settings.GRAV_OBJ_COLOR
        if not name.strip().title() == "Sun":
            self.img_scale = grav_sim.camera.img_scale
        else:
            self.img_scale = 1

        self.diameter = real_scale * self.img_scale * self.settings.SCREEN_HEIGHT
        if img_path:
            try:
                load_image = pygame.image.load(img_path).convert_alpha()
                self.image = pygame.transform.scale(load_image, (self.diameter, self.diameter))
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Check that you are in the correct directory, and the image path provided for Grav_obj is correct."
                )

        self.rect.centerx = pos_x
        self.rect.centery = pos_y

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
