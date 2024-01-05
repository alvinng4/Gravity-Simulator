import sys

import pygame.font

from textbox import Textbox


class Menu:
    """A class to build the menu"""

    def __init__(self, grav_sim):
        """Initialize button attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = grav_sim.screen.get_rect()
        self.settings = grav_sim.settings

        self.menu_active == False

        self.resume_button = Textbox(
            grav_sim,
            0.25,
            0.1,
            "Resume",
            48,
            (245, 245, 245),
            (0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - 0.3 * self.settings.screen_height,
        )
        self.void_button = Textbox(
            grav_sim,
            0.25,
            0.1,
            "Void",
            48,
            (245, 245, 245),
            (0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - 0.1 * self.settings.screen_height,
        )
        self.solar_system_button = Textbox(
            grav_sim,
            0.25,
            0.1,
            "Solar System",
            48,
            (245, 245, 245),
            (0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - (-0.1) * self.settings.screen_height,
        )
        self.exit_button = Textbox(
            grav_sim,
            0.25,
            0.1,
            "Exit",
            48,
            (245, 245, 245),
            (0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - (-0.3) * self.settings.screen_height,
        )

    def menu_active(self):
        return self.menu_active

    def draw(self):
        self.resume_button.draw()
        self.void_button.draw()
        self.solar_system_button.draw()
        self.exit_button.draw()

    def _check_button(self, mouse_pos, grav_sim):
        if self.resume_button.rect.collidepoint(mouse_pos):
            self.menu_active = False
        if self.void_button.rect.collidepoint(mouse_pos):
            grav_sim.grav_objs.empty()
            self.menu_active = False
        if self.solar_system_button.rect.collidepoint(mouse_pos):
            grav_sim.grav_objs.empty()
            grav_sim._create_solor_system()
            self.menu_active = False
        if self.exit_button.rect.collidepoint(mouse_pos):
            sys.exit()
