import pygame

class Grav_obj():
    def __init__(self, screen, img_path: str, x0, y0):
        self.screen = screen

        self.image = pygame.image.load(img_path)
        self.rect = self.image.get_rect()
        #self.screen_rect = screen.get_rect()

        self.rect.centerx = x0
        self.rect.centery = y0


    def draw(self):
        """Draw the object at its current location."""
        self.screen.blit(self.image, self.rect)










def solar_system(screen):
    """Initialize the solar system when user hit a key"""
    sun = Grav_obj(screen, "images/sun.png", 0, 0)
    mercury = Grav_obj(screen, "images/mercury.png")
    venus = Grav_obj(screen, "images/venus.png")
    earth = Grav_obj(screen, "images/earth.png")
    mars = Grav_obj(screen, "images/mars.png")
    jupiter = Grav_obj(screen, "images/jupiter.png")
    saturn = Grav_obj(screen, "images/saturn.png")
    uranus = Grav_obj(screen, "images/uranus.png")
    neptune = Grav_obj(screen, "images/neptune.png")
