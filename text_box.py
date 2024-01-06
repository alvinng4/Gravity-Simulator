import pygame.font


class Text_box:
    """A class to build text boxes"""

    def __init__(
        self,
        grav_sim,
        size_factor_x: float,
        size_factor_y: float,
        font_size: int,
        font: str = None,
        msg: str = None,
        text_box_color: tuple = None,
        text_color: tuple = (255, 255, 255),
        centerx: float = None,
        centery: float = None,
        text_box_left: float = 0,
        text_box_top: float = 0,
    ):
        """Initialize text box attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = self.screen.get_rect()

        # Set the dimensions and properties of the text box.
        self.width = size_factor_x * grav_sim.settings.screen_width
        self.height = size_factor_y * grav_sim.settings.screen_height

        self.textbox_color = text_box_color
        self.text_color = text_color
        self.font = pygame.font.SysFont(font, font_size)

        # Build the text box's rect object and center it.
        self.rect = pygame.Rect(text_box_left, text_box_top, self.width, self.height)
        if centerx:
            self.rect.centerx = centerx
        if centery:
            self.rect.centery = centery

        # The message needs to be printed only once.
        if msg:
            self.print_msg(msg)

    def print_msg(self, msg):
        """Turn msg into a rendered image and center text on the text box."""
        self.msg_image = self.font.render(
            msg, True, self.text_color, self.textbox_color
        )
        self.msg_image_rect = self.msg_image.get_rect()
        self.msg_image_rect.center = self.rect.center

    def draw(self):
        """Draw blank text box and then draw message."""
        if self.textbox_color:
            self.screen.fill(self.textbox_color, self.rect)

        self.screen.blit(self.msg_image, self.msg_image_rect)
