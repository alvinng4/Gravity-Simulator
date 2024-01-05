import pygame.font


class Textbox:
    """A class to build text boxes"""

    def __init__(
        self,
        grav_sim,
        size_factor_x: float,
        size_factor_y: float,
        msg: str,
        font_size: int,
        textbox_color: tuple = None,
        text_color: tuple = (0, 0, 0),
        centerx: float = None,
        centery: float = None,
        textbox_left: float = 0,
        textbox_top: float = 0,
    ):
        """Initialize button attributes."""
        self.screen = grav_sim.screen
        self.screen_rect = self.screen.get_rect()

        # Set the dimensions and properties of the button.
        self.width = size_factor_x * grav_sim.settings.screen_width
        self.height = size_factor_y * grav_sim.settings.screen_height

        self.textbox_color = textbox_color
        self.text_color = text_color
        self.font = pygame.font.SysFont(None, font_size)

        # Build the button's rect object and center it.
        self.rect = pygame.Rect(0, 0, self.width, self.height)
        if centerx:
            self.rect.centerx = centerx
        if centery:
            self.rect.centery = centery

        # The button message needs to be printed only once.
        self._print_msg(msg)

    def _print_msg(self, msg):
        """Turn msg into a rendered image and center text on the button."""
        self.msg_image = self.font.render(
            msg, True, self.text_color, self.textbox_color
        )
        self.msg_image_rect = self.msg_image.get_rect()
        self.msg_image_rect.center = self.rect.center

    def draw(self):
        """Draw blank button and then draw message."""
        if self.textbox_color:
            self.screen.fill(self.textbox_color, self.rect)

        self.screen.blit(self.msg_image, self.msg_image_rect)
