import sys

import pygame

from settings import Settings
from grav_obj import Grav_obj
from camera import Camera


class GravitySimulator:
    """Overall class to manage the main program."""

    def __init__(self):
        pygame.init()
        self.clock = pygame.time.Clock()
        self.settings = Settings()

        self.screen = pygame.display.set_mode((Settings.SCREEN_WIDTH, Settings.SCREEN_HEIGHT), pygame.SCALED, vsync = 1)
        pygame.display.set_caption("Gravity Simulator")
        self.camera = Camera()

        ### For testing
        self.sun = Grav_obj(self, self.screen.get_rect().centerx, self.screen.get_rect().centery, "images/sun.png")
        ###

    def run_prog(self):
        # Start the main loop for the program.
        while True:
            self._check_events()
            self._update_screen()
            self.clock.tick(60)

    def _check_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
            elif event.type == pygame.KEYDOWN:
                self._check_key_down_events(event)
            elif event.type == pygame.KEYUP:
                self._check_key_up_events(event)


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
        elif event.key == pygame.K_ESCAPE:  # Temporary
            sys.exit()

    def _update_screen(self):
        self.camera.update()
        self.screen.fill(Settings.BG_COLOR)
        ### For testing
        self.sun.draw()
        ###
        pygame.display.flip()


if __name__ == "__main__":
    grav_sim = GravitySimulator()
    grav_sim.run_prog()
