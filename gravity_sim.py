import sys

import pygame

from settings import Settings
from grav_obj import Grav_obj
from camera import Camera
from menu import Menu


class GravitySimulator:
    """Overall class to manage the main program."""

    def __init__(self):
        pygame.init()
        self.clock = pygame.time.Clock()
        self.settings = Settings()

        self.screen = pygame.display.set_mode(
            (Settings.SCREEN_WIDTH, Settings.SCREEN_HEIGHT), pygame.SCALED, vsync=1
        )
        pygame.display.set_caption("Gravity Simulator")
        self.camera = Camera(img_scale=10)

        self.menu = Menu(self)
        self.grav_objs = pygame.sprite.Group()


    def run_prog(self):
        # Start the main loop for the program.
        while True:
            self._check_events()
            self._update_grav_objs()
            self._update_screen()
            self.clock.tick(self.settings.FPS)

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
        elif event.key == pygame.K_ESCAPE:
            self.menu.menu_active = True


    def _update_screen(self):
        self.camera.update_movement()
        self.screen.fill(Settings.BG_COLOR)

        self.grav_objs.draw(self.screen)

        if self.menu.menu_active == True:
            self.menu.draw()

        pygame.display.flip()

    def _create_grav_obj(self):
        grav_obj = Grav_obj(self, 0, 0, 1)
        self.grav_objs.add(grav_obj)
    
    def _create_solor_system(self):
        centerx = self.screen.get_rect().centerx
        centery = self.screen.get_rect().centery
        sun = Grav_obj(self, centerx, centery, 0.05, "Sun", "images/sun.png")
        mercury = Grav_obj(self, 1.1 *centerx, centery, 0.0005, "Mercury", "images/mercury.png")
        venus = Grav_obj(self, 1.2 * centerx, centery, 0.0005, "Venus", "images/venus.png")
        earth = Grav_obj(self, 1.3 * centerx, centery, 0.0005, "Earth", "images/earth.png")
        mars = Grav_obj(self, 1.5 * centerx, centery, 0.0005, "Mars", "images/mars.png")
        jupiter = Grav_obj(self, 1.8 * centerx, centery, 0.0005, "Jupiter", "images/jupiter.png")
        saturn = Grav_obj(self, 2 * centerx, centery, 0.0005, "Saturn", "images/saturn.png")
        uranus = Grav_obj(self, 4 * centerx, centery, 0.0005, "Uranus", "images/uranus.png")
        neptune = Grav_obj(self, 5 * centerx, centery, 0.0005, "Neptune", "images/neptune.png")
        self.grav_objs.add(sun)
        self.grav_objs.add(mercury)
        self.grav_objs.add(venus)
        self.grav_objs.add(earth)
        self.grav_objs.add(mars)
        self.grav_objs.add(jupiter)
        self.grav_objs.add(saturn)
        self.grav_objs.add(uranus)
        self.grav_objs.add(neptune)

    def _update_grav_objs(self):
        """Update the apparent position of all grav_objs with camera"""
        self.grav_objs.update()


if __name__ == "__main__":
    grav_sim = GravitySimulator()
    grav_sim.run_prog()
