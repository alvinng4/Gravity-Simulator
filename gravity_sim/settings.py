class Settings:
    """A class to store all settings for gravity simulator."""

    DEFAULT_RESOLUTION = (1920, 1080)
    DEFAULT_STAR_IMG_SCALE = 5000
    DEFAULT_PLANET_IMG_SCALE = 100000
    DEFAULT_DISTANCE_SCALE = 200
    DEFAULT_TIME_SPEED = 1
    DEFAULT_dt = 0.1
    MAX_FPS = 60
    BG_COLOR = (0, 0, 0)  # Background color

    CHANGE_STAR_IMG_SCALE_SPEED = 500
    CHANGE_IMG_SCALE_SPEED = 10000
    CHANGE_DISTANCE_SCALE_SPEED = 10
    CHANGE_dt_SPEED = 0.1
    CHANGE_TIME_SPEED_SPEED = 10

    DEFAULT_NEW_OBJECT_VELOCITY_SCALE = 0.0001

    def __init__(
        self,
        screen_width: int,
        screen_height: int,
        time_speed: int = DEFAULT_TIME_SPEED,
        dt: float = DEFAULT_dt,
    ):
        # To change the default settings of screen_width and screen_height, go to _read_command_line_arg function in __main__.py
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.star_img_scale = self.DEFAULT_STAR_IMG_SCALE
        self.img_scale = self.DEFAULT_PLANET_IMG_SCALE
        self.time_speed = time_speed
        self.dt = dt
        self.distance_scale = self.DEFAULT_DISTANCE_SCALE

    @property
    def screen_width(self):
        return self._screen_width

    @property
    def screen_height(self):
        return self._screen_height

    @property
    def star_img_scale(self):
        return self._star_img_scale

    @property
    def img_scale(self):
        return self._img_scale

    @property
    def distance_scale(self):
        return self._distance_scale

    @property
    def dt(self):
        return self._dt

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
    @star_img_scale.setter
    def star_img_scale(self, value):
        if value > 100000:
            self._star_img_scale = 100000
        elif value < 0:
            self._star_img_scale = 0
        else:
            self._star_img_scale = value

    # Img may corrupt if the scale is too large.
    @img_scale.setter
    def img_scale(self, value):
        if value > 5000000:
            self._img_scale = 5000000
        elif value < 0:
            self._img_scale = 0
        else:
            self._img_scale = value

    @distance_scale.setter
    def distance_scale(self, value):
        if value > 1000:
            self._distance_scale = 1000
        elif value <= 10:
            self._distance_scale = 10
        else:
            self._distance_scale = value

    @dt.setter
    def dt(self, value):
        if value < 0:
            self._dt = 0
        else:
            self._dt = value

    @time_speed.setter
    def time_speed(self, value):
        if value > 10000:
            self._time_speed = 10000
        elif value < 1:
            self._time_speed = 1
        else:
            self._time_speed = int(value)
