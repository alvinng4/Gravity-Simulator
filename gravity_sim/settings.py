class Settings:
    """A class to store all settings for gravity simulator."""

    def __init__(
        self,
        screen_width,
        screen_height,
        sun_img_scale,
        img_scale,
        time_speed: int = 512,
    ):
        # To change the default settings of screen_width, screen_height, sun_img_scale, img_scale, go to _read_command_line_arg function in __main__.py
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.bg_color = (0, 0, 0)  # Background color
        self.fps = 60

        self.grav_obj_color = (255, 255, 255)
        self.sun_img_scale = sun_img_scale
        self.img_scale = img_scale
        self.distance_scale = 1
        self.time_speed = time_speed

    @property
    def screen_width(self):
        return self._screen_width

    @property
    def screen_height(self):
        return self._screen_height

    @property
    def sun_img_scale(self):
        return self._sun_img_scale

    @property
    def img_scale(self):
        return self._img_scale

    @property
    def distance_scale(self):
        return self._distance_scale

    @property
    def time_speed(self):
        return self._time_speed

    @screen_width.setter
    def screen_width(self, value):
        if value < 0:
            self._screen_width = 0
        else:
            self._screen_width = value

    @screen_height.setter
    def screen_height(self, value):
        if value < 0:
            self._screen_height = 0
        else:
            self._screen_height = value

    # Img may corrupt if the scale is too large.
    @sun_img_scale.setter
    def sun_img_scale(self, value):
        if value > 100:
            self._sun_img_scale = 100
        elif value < 0:
            self._sun_img_scale = 0
        else:
            self._sun_img_scale = value

    # Img may corrupt if the scale is too large.
    @img_scale.setter
    def img_scale(self, value):
        if value > 1000:
            self._img_scale = 1000
        elif value < 0:
            self._img_scale = 0
        else:
            self._img_scale = value

    @distance_scale.setter
    def distance_scale(self, value):
        if value > 1000:
            self._distance_scale = 1000
        elif value <= 0.1:
            self._distance_scale = 0.1
        else:
            self._distance_scale = value

    @time_speed.setter
    def time_speed(self, value):
        if value > 512:
            self._time_speed = 512
        elif value < 1:
            self._time_speed = 1
        else:
            self._time_speed = int(value)
