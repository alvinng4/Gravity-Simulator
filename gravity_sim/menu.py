import sys

from text_box import Text_box
from grav_obj import Grav_obj


class Menu:
    """A class to build the menu"""

    def __init__(self, grav_sim):
        """Initialize button attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = self.screen.get_rect()
        self.settings = grav_sim.settings
        self.main_menu_active = True
        self.menu_active = True

        self.main_menu_caption = Text_box(
            grav_sim,
            72,
            0,
            0,
            msg="2D N-body Gravity Simulator",
            font="Manrope",
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - 0.3 * self.settings.screen_height,
            ),
        )
        self.resume_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Resume",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - 0.22 * self.settings.screen_height,
            ),
        )
        self.void_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Void",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - 0.14 * self.settings.screen_height,
            ),
        )
        self.solar_system_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Solar System",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - 0.06 * self.settings.screen_height,
            ),
        )
        self.figure_8_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Figure 8 orbit",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - (-0.02) * self.settings.screen_height,
            ),
        )
        self.pyth_3_body_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Pythagorean three-body",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - (-0.1) * self.settings.screen_height,
            ),
        )
        self.exit_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Exit",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - (-0.18) * self.settings.screen_height,
            ),
        )
        self.main_menu_button = Text_box(
            grav_sim,
            48,
            0.25,
            0.05,
            msg="Main Menu",
            text_box_color=(220, 220, 220),
            text_color=(0, 0, 0),
            center=(
                self.screen_rect.centerx,
                self.screen_rect.centery - (-0.18) * self.settings.screen_height,
            ),
        )

    def draw(self):
        """Draw the menu buttons"""
        if self.main_menu_active == True:
            self.main_menu_caption.draw()
            self.exit_button.draw()
        else:
            self.resume_button.draw()
            self.main_menu_button.draw()

        self.void_button.draw()
        self.solar_system_button.draw()
        self.figure_8_button.draw()
        self.pyth_3_body_button.draw()

    def check_button(self, grav_sim, mouse_pos):
        """Check if there is any click on the buttons"""
        if self.main_menu_active == False:
            if self.resume_button.rect.collidepoint(mouse_pos):
                self.menu_active = False
            if self.main_menu_button.rect.collidepoint(mouse_pos):
                grav_sim.grav_objs.empty()
                grav_sim.stats.reset(grav_sim)
                self.main_menu_active = True
        else:
            if self.exit_button.rect.collidepoint(mouse_pos):
                sys.exit()
        if self.void_button.rect.collidepoint(mouse_pos):
            self._menu_common_actions(grav_sim)
        if self.solar_system_button.rect.collidepoint(mouse_pos):
            self._menu_common_actions(grav_sim)
            Grav_obj.create_solor_system(grav_sim)
        if self.figure_8_button.rect.collidepoint(mouse_pos):
            self._menu_common_actions(grav_sim)
            Grav_obj.create_figure_8(grav_sim)
        if self.pyth_3_body_button.rect.collidepoint(mouse_pos):
            self._menu_common_actions(grav_sim)
            Grav_obj.create_pyth_3_body(grav_sim)

    def _menu_common_actions(self, grav_sim):
        grav_sim.grav_objs.empty()
        grav_sim.stats.reset(grav_sim)
        self.menu_active = False
        self.main_menu_active = False
