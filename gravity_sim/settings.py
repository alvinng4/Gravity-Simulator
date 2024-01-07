class Settings:
    """A class to store all settings for gravity simulator."""

    def __init__(self, screen_width, screen_height, sun_img_scale, img_scale):
        # To change the default settings of screen_width, screen_height, sun_img_scale, img_scale, go to _read_command_line_arg function in __main__.py
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.bg_color = (0, 0, 0)  # Background color
        self.fps = 60

        self.grav_obj_color = (255, 255, 255)
        self.sun_img_scale = sun_img_scale
        self.img_scale = img_scale

    @property
    def sun_img_scale(self):
        return self._sun_img_scale

    # Img may corrupt if the scale is too large.
    @sun_img_scale.setter
    def sun_img_scale(self, value):
        if value > 100:
            self._sun_img_scale = 100
        elif value < 0:
            self._sun_img_scale = 0
        else:
            self._sun_img_scale = value

    @property
    def img_scale(self):
        return self._img_scale

    # Img may corrupt if the scale is too large.
    @img_scale.setter
    def img_scale(self, value):
        if value > 100:
            self._img_scale = 100
        elif value < 0:
            self._img_scale = 0
        else:
            self._img_scale = value
