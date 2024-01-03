import pygame

class Grav_obj():
    def __init__(self, screen, img_path: str, x0, y0):
        self.screen = screen

        self.image = pygame.image.load(img_path)
        self.rect = self.image.get_rect()
        #self.screen_rect = screen.get_rect()

        self.rect.centerx = x0
        self.rect.centery = y0


    def blitme(self):
        """Draw the object at its current location."""
        self.screen.blit(self.image, self.rect)





def solar_system(screen):
    """Initialize the solar system when user hit a key"""
    sun = Grav_obj(screen, "images/sun.png")
    mercury = Grav_obj()
    venus = Grav_obj()
    earth = Grav_obj()
    mars = Grav_obj()
    jupiter = Grav_obj()
    saturn = Grav_obj()
    uranus = Grav_obj()
    neptune = Grav_obj()
