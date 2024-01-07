import sys
import argparse

import pygame

from settings import Settings
from grav_obj import Grav_obj
from camera import Camera
from menu import Menu
from stats import Stats


class GravitySimulator:
    """Overall class to manage the main program."""

    def __init__(self):
        self._read_command_line_arg()
        pygame.init()
        self.settings = Settings(
            screen_width=self.args.resolution[0],
            screen_height=self.args.resolution[1],
            sun_img_scale=self.args.img_scale[0],
            img_scale=self.args.img_scale[1],
        )
        self.screen = pygame.display.set_mode(
            (self.settings.screen_width, self.settings.screen_height),
            pygame.SCALED,
            vsync=1,
        )
        pygame.display.set_caption("Gravity Simulator")
        self.clock = pygame.time.Clock()
        self.menu = Menu(self)
        self.stats = Stats(self)
        self.camera = Camera()
        self.grav_objs = pygame.sprite.Group()

    def run_prog(self):
        # Start the main loop for the program.
        while True:
            self._check_events()
            self._update_screen()
            self.clock.tick(60)

    def _check_events(self):
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN:
                self._check_key_down_events(event)
            elif event.type == pygame.KEYUP:
                self._check_key_up_events(event)
            elif event.type == pygame.MOUSEBUTTONDOWN:
                mouse_pos = pygame.mouse.get_pos()
                if self.menu.menu_active == True:
                    self.menu._check_button(mouse_pos, self)
            elif event.type == pygame.QUIT:
                sys.exit()

    def _check_key_up_events(self, event):
        if event.key == pygame.K_d:
            self.camera.moving_right = False
        elif event.key == pygame.K_a:
            self.camera.moving_left = False
        elif event.key == pygame.K_w:
            self.camera.moving_up = False
        elif event.key == pygame.K_s:
            self.camera.moving_down = False

    def _check_key_down_events(self, event):
        if event.key == pygame.K_d:
            self.camera.moving_right = True
        elif event.key == pygame.K_a:
            self.camera.moving_left = True
        elif event.key == pygame.K_w:
            self.camera.moving_up = True
        elif event.key == pygame.K_s:
            self.camera.moving_down = True
        elif event.key == pygame.K_f:
            pygame.display.toggle_fullscreen()
        elif event.key == pygame.K_ESCAPE:
            if self.menu.start_menu_active == False:
                self.menu.menu_active = not self.menu.menu_active

    def _update_screen(self):
        self.camera.update_movement()
        self.grav_objs.update()
        self.screen.fill(self.settings.bg_color)
        self.grav_objs.draw(self.screen)
        self.stats.update(self)
        if self.menu.menu_active == True:
            self.menu.draw()
        pygame.display.flip()

    def _read_command_line_arg(self):
        parser = argparse.ArgumentParser(description="2D N-body gravity simulator")
        parser.add_argument(
            "--resolution",
            "-r",
            nargs=2,
            default=[1920, 1080],
            type=float,
            help="Usage: --resolution width, height",
        )
        parser.add_argument(
            "--img_scale",
            "-i",
            nargs=2,
            default=[10, 30],
            type=float,
            help="Usage: --img_scale Solar image scale, obj image scale",
        )
        self.args = parser.parse_args()
        if self.args.resolution[0] != 0 and self.args.resolution[1] != 0:
            pass
        else:
            sys.exit("Invalid resolution")


if __name__ == "__main__":
    grav_sim = GravitySimulator()
    grav_sim.run_prog()
