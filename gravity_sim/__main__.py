import sys
import argparse
from os import environ

# Remove the "Hello from the pygame community." message when starting the program.
environ["PYGAME_HIDE_SUPPORT_PROMPT"] = "1"

import pygame

from settings import Settings
from camera import Camera
from menu import Menu
from stats import Stats
from grav_obj import Grav_obj
from simulator import Simulator


class GravitySimulator:
    """Overall class to manage the main program."""

    def __init__(self):
        """
        Initialize the main program.
        Initialization dependencies:
            settings: none
            menu: none
            camera: none
            stats: settings
            grav_objs: camera, settings
            simulator: stats, settings
        """
        self._read_command_line_arg()
        pygame.init()
        self.settings = Settings(
            screen_width=self.args.resolution[0],
            screen_height=self.args.resolution[1],
        )
        self.screen = pygame.display.set_mode(
            (self.settings.screen_width, self.settings.screen_height),
            pygame.SCALED,
            vsync=1,
        )
        pygame.display.set_caption("Gravity Simulator")
        self.clock = pygame.time.Clock()
        self.menu = Menu(self)
        self.camera = Camera()
        self.stats = Stats(self)
        self.grav_objs = pygame.sprite.Group()
        self.simulator = Simulator(self)

    def run_prog(self):
        """The main loop for the program"""
        while True:
            self._check_events()
            self._update_events()
            self._simulation()
            self._update_screen()
            self.clock.tick(self.settings.MAX_FPS)

    def _check_events(self):
        self.simulator.check_current_integrator()
        self.settings.check_current_changing_parameter()
        for event in pygame.event.get():
            match event.type:
                case pygame.KEYDOWN:
                    self._check_key_down_events(event)
                case pygame.KEYUP:
                    self._check_key_up_events(event)
                case pygame.MOUSEBUTTONDOWN:
                    self._check_mouse_button_down_events(event)
                case pygame.MOUSEBUTTONUP:
                    self._check_mouse_button_up_events(event)
                case pygame.MOUSEWHEEL:
                    self.settings.scroll_change_parameters(event.y)
                case pygame.QUIT:
                    sys.exit()

    def _update_events(self):
        self.camera.update_movement()
        self.grav_objs.update()
        self.stats.update(self)

    def _simulation(self):
        if self.grav_objs and not self.stats.is_paused:
            self.simulator.run_simulation(self)
            self.simulator.unload_value(self)

    def _update_screen(self):
        self.screen.fill(Settings.BG_COLOR)
        self.grav_objs.draw(self.screen)
        if self.settings.is_hide_gui == False:
            self.stats.draw(self)
        if self.stats.is_holding_rclick == True:
            self._new_star_draw_line_circle()
        if self.menu.menu_active == True:
            self.menu.draw()
        pygame.display.flip()

    def _check_key_up_events(self, event):
        match event.key:
            case pygame.K_d:
                self.camera.moving_right = False
            case pygame.K_a:
                self.camera.moving_left = False
            case pygame.K_w:
                self.camera.moving_up = False
            case pygame.K_s:
                self.camera.moving_down = False

    def _check_key_down_events(self, event):
        match event.key:
            case pygame.K_d:
                self.camera.moving_right = True
            case pygame.K_a:
                self.camera.moving_left = True
            case pygame.K_w:
                self.camera.moving_up = True
            case pygame.K_s:
                self.camera.moving_down = True
            case pygame.K_p:
                if self.stats.is_paused == False:
                    self.stats.start_pause()
                elif self.stats.is_paused == True:
                    self.stats.end_pause()
            case pygame.K_f:
                pygame.display.toggle_fullscreen()
            case pygame.K_h:
                self.settings.is_hide_gui = not self.settings.is_hide_gui
            case pygame.K_r:
                self.settings.reset_parameters()
            case pygame.K_ESCAPE:
                if self.menu.main_menu_active == False:
                    self.menu.menu_active = not self.menu.menu_active

    def _check_mouse_button_down_events(self, event):
        if event.button == 1:  # left click
            mouse_pos = pygame.mouse.get_pos()
            self.stats.check_button(self, mouse_pos)
            if self.menu.menu_active == True:
                self.menu.check_button(self, mouse_pos)
        elif event.button == 3:  # right click
            if self.menu.menu_active == False:
                mouse_pos = pygame.mouse.get_pos()
                self.stats.start_holding_rclick()
                self.new_star_mouse_pos = mouse_pos
                self.new_star_camera_pos = self.camera.pos

    def _check_mouse_button_up_events(self, event):
        if event.button == 3:  # right click up
            if self.stats.is_holding_rclick == True:
                self.new_star_drag_mouse_pos = (
                    pygame.mouse.get_pos()
                )  # for object's velocity
                self.stats.end_holding_rclick()
                Grav_obj.create_star(
                    self,
                    self.new_star_mouse_pos,
                    self.new_star_camera_pos,
                    self.new_star_drag_mouse_pos,
                    self.camera.pos,
                )

    def _read_command_line_arg(self):
        parser = argparse.ArgumentParser(description="N-body gravity simulator")
        parser.add_argument(
            "--resolution",
            "-r",
            nargs=2,
            default=Settings.DEFAULT_RESOLUTION,
            type=float,
            help="Usage: --resolution <width>, <height>",
        )
        self.args = parser.parse_args()
        if not (self.args.resolution[0] > 0 and self.args.resolution[1] > 0):
            sys.exit("Invalid resolution")

    def _new_star_draw_line_circle(self):
        pygame.draw.line(
            self.screen,
            "white",
            (
                self.new_star_mouse_pos[0]
                + (self.new_star_camera_pos[0] - self.camera.pos[0]),
                self.new_star_mouse_pos[1]
                + (self.new_star_camera_pos[1] - self.camera.pos[1]),
            ),
            pygame.mouse.get_pos(),
        )
        m = 1 * 0.5 * self.stats.holding_rclick_time * self.settings.new_star_mass_scale
        R = Grav_obj.SOLAR_RADIUS * (m ** (1.0 / 3.0))
        img_R = (
            R
            * (699.0 / 894.0)  # Actual Sun size in images/sun.png with size (894 x 894)
            * self.settings.star_img_scale
        )
        new_star_circle_pos = [
            self.new_star_mouse_pos[0]
            + (self.new_star_camera_pos[0] - self.camera.pos[0]),
            self.new_star_mouse_pos[1]
            + (self.new_star_camera_pos[1] - self.camera.pos[1]),
        ]
        pygame.draw.circle(self.screen, "orange", new_star_circle_pos, img_R, width=1)


if __name__ == "__main__":
    grav_sim = GravitySimulator()
    grav_sim.run_prog()
