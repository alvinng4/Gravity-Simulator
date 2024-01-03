import sys

import pygame

from settings import Settings
from grav_obj import Grav_obj

def check_events():
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            sys.exit()


def update_screen(screen, *grav_objs):
    screen.fill(Settings.BG_COLOR)
    for grav_obj in grav_objs:
        grav_obj.blitme()


    # Make the most recently drawn screen visible.
    pygame.display.flip()
