import sys

import pygame.font

from text_box import Text_box
from grav_obj import Grav_obj


class Menu:
    """A class to build the menu"""

    def __init__(self, grav_sim):
        """Initialize button attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = grav_sim.screen.get_rect()
        self.settings = grav_sim.settings
        self.start_menu_active = True
        self.menu_active = True

        self.start_menu_caption = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="2D N-body Gravity Simulator",
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - 0.3 * self.settings.screen_height,
        )
        self.resume_button = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="Resume",
            text_box_color=(245, 245, 245),
            text_color=(0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - 0.3 * self.settings.screen_height,
        )
        self.void_button = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="Void",
            text_box_color=(245, 245, 245),
            text_color=(0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - 0.15 * self.settings.screen_height,
        )
        self.solar_system_button = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="Solar System",
            text_box_color=(245, 245, 245),
            text_color=(0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery,
        )
        self.figure_8_button = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="Figure 8",
            text_box_color=(245, 245, 245),
            text_color=(0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - (-0.15) * self.settings.screen_height,
        )
        self.exit_button = Text_box(
            grav_sim,
            0.25,
            0.08,
            48,
            msg="Exit",
            text_box_color=(245, 245, 245),
            text_color=(0, 0, 0),
            centerx=self.screen_rect.centerx,
            centery=self.screen_rect.centery - (-0.3) * self.settings.screen_height,
        )

    def menu_active(self):
        return self.menu_active

    def draw(self):
        if self.start_menu_active == True:
            self.start_menu_caption.draw()
        else:
            self.resume_button.draw()
        self.void_button.draw()
        self.solar_system_button.draw()
        self.figure_8_button.draw()
        self.exit_button.draw()

    def _check_button(self, mouse_pos, grav_sim):
        if self.resume_button.rect.collidepoint(mouse_pos):
            self.menu_active = False
        if self.void_button.rect.collidepoint(mouse_pos):
            grav_sim.grav_objs.empty()
            self.menu_active = False
            self.start_menu_active = False
        if self.solar_system_button.rect.collidepoint(mouse_pos):
            grav_sim.grav_objs.empty()
            Grav_obj.create_solor_system(grav_sim)
            self.menu_active = False
            self.start_menu_active = False
        if self.figure_8_button.rect.collidepoint(mouse_pos):
            grav_sim.grav_objs.empty()
            Grav_obj.create_figure_8(grav_sim)
            self.menu_active = False
            self.start_menu_active = False
        if self.exit_button.rect.collidepoint(mouse_pos):
            sys.exit()
