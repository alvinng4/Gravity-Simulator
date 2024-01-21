import math


class Settings:
    """A class to store all settings for gravity simulator."""

    DEFAULT_RESOLUTION = (1920, 1080)

    MAX_FPS = 60
    BG_COLOR = (0, 0, 0)  # Background color

    DEFAULT_STAR_IMG_SCALE = 5000
    DEFAULT_PLANET_IMG_SCALE = 100000
    DEFAULT_DISTANCE_SCALE = 200
    DEFAULT_dt = 0.1
    DEFAULT_TIME_SPEED = 1
    DEFAULT_EPSILON = 1e-2
    DEFAULT_NEW_STAR_MASS_SCALE = 1

    MAX_STAR_IMG_SCALE = 100000
    MIN_STAR_IMG_SCALE = 0
    MAX_PLANET_IMG_SCALE = 500000
    MIN_PLANET_IMG_SCALE = 0
    MAX_DISTANCE_SCALE = 1000
    MIN_DISTANCE_SCALE = 0
    MAX_DT = 100
    MIN_DT = 1e-10
    MAX_TIME_SPEED = 10000
    MIN_TIME_SPEED = 1
    MAX_EPSILON = 1e2
    MIN_EPSILON = 1e-15
    MAX_NEW_STAR_MASS_SCALE = 1000
    MIN_NEW_STAR_MASS_SCALE = 1e-3

    DEFAULT_CHANGE_STAR_IMG_SCALE_SPEED = 500
    DEFAULT_CHANGE_PLANET_IMG_SCALE_SPEED = 10000
    DEFAULT_CHANGE_DISTANCE_SCALE_SPEED = 10

    DEFAULT_NEW_STAR_VELOCITY_SCALE = 0.0001

    def __init__(
        self,
        screen_width: int,
        screen_height: int,
    ):
        # To change the default settings of screen_width and screen_height, go to _read_command_line_arg function in __main__.py
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.star_img_scale = self.DEFAULT_STAR_IMG_SCALE
        self.planet_img_scale = self.DEFAULT_PLANET_IMG_SCALE
        self.time_speed = self.DEFAULT_TIME_SPEED
        self.dt = self.DEFAULT_dt
        self.distance_scale = self.DEFAULT_DISTANCE_SCALE
        self.epsilon = self.DEFAULT_EPSILON
        self.new_star_mass_scale = self.DEFAULT_NEW_STAR_MASS_SCALE
        self.set_all_parameters_changing_false()
        self.current_changing_parameter = None
        self.is_hide_gui = False

    def scroll_change_parameters(self, magnitude):
        match self.current_changing_parameter:
            case "star_img_scale":
                self.star_img_scale += (
                    self.DEFAULT_CHANGE_STAR_IMG_SCALE_SPEED * magnitude
                )
            case "planet_img_scale":
                self.planet_img_scale += (
                    self.DEFAULT_CHANGE_PLANET_IMG_SCALE_SPEED * magnitude
                )
            case "distance_scale":
                self.distance_scale += (
                    self.DEFAULT_CHANGE_DISTANCE_SCALE_SPEED * magnitude
                )
            case "new_star_mass_scale":
                if magnitude > 0:
                    for _ in range(magnitude):
                        self.new_star_mass_scale += self._inc_speed(
                            self.new_star_mass_scale
                        )
                elif magnitude < 0:
                    for _ in range(-magnitude):
                        self.new_star_mass_scale -= self._dec_speed(
                            self.new_star_mass_scale
                        )
            case "dt":
                if magnitude > 0:
                    for _ in range(magnitude):
                        self.dt += self._inc_speed(self.dt)
                elif magnitude < 0:
                    for _ in range(-magnitude):
                        self.dt -= self._dec_speed(self.dt)
            case "time_speed":
                if magnitude > 0:
                    for _ in range(magnitude):
                        self.time_speed += self._inc_speed(self.time_speed)
                elif magnitude < 0:
                    for _ in range(-magnitude):
                        self.time_speed -= self._dec_speed(self.time_speed)
            case "epsilon":
                if magnitude > 0:
                    for _ in range(magnitude):
                        self.epsilon += self._inc_speed(self.epsilon)
                elif magnitude < 0:
                    for _ in range(-magnitude):
                        self.epsilon -= self._dec_speed(self.epsilon)

    @staticmethod
    def _inc_speed(x: float) -> float:
        return 10 ** math.floor(math.log10(x))

    @staticmethod
    def _dec_speed(x: float) -> float:
        if f"{x:e}"[0] == "1":
            return x * 0.1
        else:
            return 10 ** math.floor(math.log10(x))

    def check_current_changing_parameter(self):
        if self.is_changing_star_img_scale == True:
            self.current_changing_parameter = "star_img_scale"
        elif self.is_changing_planet_img_scale == True:
            self.current_changing_parameter = "planet_img_scale"
        elif self.is_changing_distance_scale == True:
            self.current_changing_parameter = "distance_scale"
        elif self.is_changing_dt == True:
            self.current_changing_parameter = "dt"
        elif self.is_changing_time_speed == True:
            self.current_changing_parameter = "time_speed"
        elif self.is_changing_epsilon == True:
            self.current_changing_parameter = "epsilon"
        elif self.is_changing_new_star_mass_scale == True:
            self.current_changing_parameter = "new_star_mass_scale"

    def set_all_parameters_changing_false(self):
        self.is_changing_star_img_scale = False
        self.is_changing_planet_img_scale = False
        self.is_changing_distance_scale = False
        self.is_changing_dt = False
        self.is_changing_time_speed = False
        self.is_changing_epsilon = False
        self.is_changing_new_star_mass_scale = False

    def reset_parameters(self):
        self.star_img_scale = self.DEFAULT_STAR_IMG_SCALE
        self.planet_img_scale = self.DEFAULT_PLANET_IMG_SCALE
        self.time_speed = self.DEFAULT_TIME_SPEED
        self.dt = self.DEFAULT_dt
        self.distance_scale = self.DEFAULT_DISTANCE_SCALE
        self.epsilon = self.DEFAULT_EPSILON
        self.new_star_mass_scale = self.DEFAULT_NEW_STAR_MASS_SCALE

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
    def planet_img_scale(self):
        return self._planet_img_scale

    @property
    def distance_scale(self):
        return self._distance_scale

    @property
    def dt(self):
        return self._dt

    @property
    def time_speed(self):
        return self._time_speed

    @property
    def epsilon(self):
        return self._epsilon

    @property
    def new_star_mass_scale(self):
        return self._new_star_mass_scale

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
        if value > self.MAX_STAR_IMG_SCALE:
            self._star_img_scale = self.MAX_STAR_IMG_SCALE
        elif value < self.MIN_STAR_IMG_SCALE:
            self._star_img_scale = self.MIN_STAR_IMG_SCALE
        else:
            self._star_img_scale = value

    # Img may corrupt if the scale is too large.
    @planet_img_scale.setter
    def planet_img_scale(self, value):
        if value > self.MAX_PLANET_IMG_SCALE:
            self._planet_img_scale = self.MAX_PLANET_IMG_SCALE
        elif value < self.MIN_PLANET_IMG_SCALE:
            self._planet_img_scale = self.MIN_PLANET_IMG_SCALE
        else:
            self._planet_img_scale = value

    @distance_scale.setter
    def distance_scale(self, value):
        if value > self.MAX_DISTANCE_SCALE:
            self._distance_scale = self.MAX_DISTANCE_SCALE
        elif value <= self.MIN_DISTANCE_SCALE:
            self._distance_scale = self.MIN_DISTANCE_SCALE
        else:
            self._distance_scale = value

    @dt.setter
    def dt(self, value):
        if value > self.MAX_DT:
            self._dt = self.MAX_DT
        elif value < self.MIN_DT:
            self._dt = self.MIN_DT
        else:
            self._dt = round(value, ndigits=15)

    @time_speed.setter
    def time_speed(self, value):
        if value > self.MAX_TIME_SPEED:
            self._time_speed = self.MAX_TIME_SPEED
        elif value < self.MIN_TIME_SPEED:
            self._time_speed = self.MIN_TIME_SPEED
        else:
            self._time_speed = int(value)

    @epsilon.setter
    def epsilon(self, value):
        if value > self.MAX_EPSILON:
            self._epsilon = self.MAX_EPSILON
        elif value < self.MIN_EPSILON:
            self._epsilon = self.MIN_EPSILON
        else:
            self._epsilon = round(value, ndigits=15)

    @new_star_mass_scale.setter
    def new_star_mass_scale(self, value):
        if value > self.MAX_NEW_STAR_MASS_SCALE:
            self._new_star_mass_scale = self.MAX_NEW_STAR_MASS_SCALE
        elif value < self.MIN_NEW_STAR_MASS_SCALE:
            self._new_star_mass_scale = self.MIN_NEW_STAR_MASS_SCALE
        else:
            self._new_star_mass_scale = round(value, ndigits=15)
