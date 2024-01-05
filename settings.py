class Settings:
    """A class to store all settings for gravity simulator."""

    def __init__(self):
        self.screen_width = 1366
        self.screen_height = 768
        self.bg_color = (0, 0, 0)  # Background color

        self.fps = 60

        self.grav_obj_color = (255, 255, 255)
