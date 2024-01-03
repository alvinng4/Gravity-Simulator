import sys

import pygame

from settings import Settings
from grav_obj import Grav_obj
import program_functions as pf


def main():
    pygame.init()
    settings = Settings()

    screen = pygame.display.set_mode((Settings.SCREEN_WIDTH, Settings.SCREEN_HEIGHT))
    pygame.display.set_caption("Gravity Simulator")

    sun = Grav_obj(screen, "images/sun.png", 0, 0)

    # Start the main loop for the program.
    while True:
        pf.check_events()
        pf.update_screen(screen, sun)


if __name__ == "__main__":
    main()
