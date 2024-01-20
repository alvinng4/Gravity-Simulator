import os

import pygame.font


class Text_box:
    """A class to build text boxes"""

    def __init__(
        self,
        grav_sim,
        font_size: int,
        size_factor_x: float = None,
        size_factor_y: float = None,
        size_x: int = None,
        size_y: int = None,
        font: str = None,
        msg: str = None,
        text_box_color: tuple = None,
        text_color: tuple = (255, 255, 255),
        center: tuple = None,
        text_box_left_top: tuple = (0, 0),
    ):
        """Initialize text box attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = self.screen.get_rect()

        # Set the dimensions and properties of the text box.
        if size_factor_x != None:
            self.width = size_factor_x * grav_sim.settings.screen_width
        else:
            self.width = size_x
        if size_factor_y != None:
            self.height = size_factor_y * grav_sim.settings.screen_height
        else:
            self.height = size_y

        self.textbox_color = text_box_color
        self.text_color = text_color
        main_dir_path = os.path.dirname(__file__)
        path_manrope = os.path.join(main_dir_path, "assets/fonts/Manrope-Regular.ttf")
        if font == "Manrope":
            self.font = pygame.font.Font(path_manrope, font_size)
        else:
            self.font = pygame.font.SysFont(font, font_size)

        # Build the text box's rect object and center it.
        self.center = center
        self.text_box_left_top = text_box_left_top
        self.rect = pygame.Rect(
            self.text_box_left_top[0],
            self.text_box_left_top[1],
            self.width,
            self.height,
        )
        if self.center:
            self.rect.center = self.center

        # The message needs to be printed only once.
        if msg:
            self.print_msg(msg)

    def print_msg(self, msg):
        """Turn msg into a rendered image and center text on the text box."""
        self.msg_image = self.font.render(
            msg, True, self.text_color, self.textbox_color
        )
        self.msg_image_rect = self.msg_image.get_rect()
        if self.center:
            self.msg_image_rect.center = self.rect.center
        else:
            self.msg_image_rect.left = self.text_box_left_top[0]
            self.msg_image_rect.top = self.text_box_left_top[1]

    def draw(self):
        """Draw blank text box and then draw message."""
        if self.textbox_color:
            self.screen.fill(self.textbox_color, self.rect)

        self.screen.blit(self.msg_image, self.msg_image_rect)
