import sys

import pygame

from settings import Settings
from grav_obj import Grav_obj

# Not finished. May look for sprite.
class Camera():
    def __init__(self):
        self.pos_x = 0
        self.pos_y = 0

        self.speed_x = 1
        self.speed_y = 1

        # Movement flag
        self.moving_right = False
        self.moving_left = False
        self.moving_up = False
        self.moving_down = False

    def speed_x(self):
        return self.speed_x 
    
    def speed_y(self):
        return self.speed_y
    
    def update(self):
        if self.moving_right == True:
            self.pos_x += self.speed_x
        if self.moving_left == True:
            self.pos_x -= self.speed_x
        if self.moving_up == True:
            self.pos_y += self.speed_y
        if self.moving_down == True:
            self.pos_y -= self.speed_y

                 
def check_key_up_events(event, camera):
    if event.type == pygame.K_d or event.type == pygame.K_RIGHT:
        camera.moving_right == False
    elif event.type == pygame.K_a or event.type == pygame.K_LEFT:
        camera.moving_left == False       
    elif event.type == pygame.K_w or event.type == pygame.K_UP:
        camera.moving_up == False
    elif event.type == pygame.K_d or event.type == pygame.K_DOWN:
        camera.moving_down == False


def check_key_down_events(event, camera):
    if event.type == pygame.K_d or event.type == pygame.K_RIGHT:
        camera.moving_right == True
    elif event.type == pygame.K_a or event.type == pygame.K_LEFT:
        camera.moving_left == True        
    elif event.type == pygame.K_w or event.type == pygame.K_UP:
        camera.moving_up == True
    elif event.type == pygame.K_d or event.type == pygame.K_DOWN:
        camera.moving_down == True
    elif event.key == pygame.K_ESCAPE:
        sys.exit()


def check_events(camera):
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            sys.exit()
        elif event.type == pygame.KEYDOWN:   
            check_key_down_events(event, camera)
        elif event.type == pygame.KEYUP:
            check_key_up_events(event, camera)


def update_screen(screen, *grav_objs):
    screen.fill(Settings.BG_COLOR)
    for grav_obj in grav_objs:
        grav_obj.draw()


    # Make the most recently drawn screen visible.
    pygame.display.flip()




