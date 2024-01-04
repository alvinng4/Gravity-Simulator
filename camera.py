import sys

import pygame

from settings import Settings

class Camera():
    def __init__(self, grav_sim):
        self.screen_rect = grav_sim.screen.get_rect()

        self.pos_x = self.screen_rect.centerx
        self.pos_y = self.screen_rect.centery

        self.speed_x = 10
        self.speed_y = 10

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

                 






    




