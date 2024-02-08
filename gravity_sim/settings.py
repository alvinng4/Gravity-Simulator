import math

import pygame

class Settings:
    """A class to store all settings for gravity simulator."""

    MAX_FPS = 60
    BG_COLOR = (0, 0, 0)  # Background color

    DEFAULT_STAR_IMG_SCALE = 5000
    DEFAULT_PLANET_IMG_SCALE = 100000
    DEFAULT_DISTANCE_SCALE = 150
    DEFAULT_NEW_STAR_MASS_SCALE = 1
    DEFAULT_dt = 0.1
    DEFAULT_TIME_SPEED = 1
    DEFAULT_RK_MAX_ITERATION = 10
    DEFAULT_RK_MIN_ITERATION = 1
    DEFAULT_TOLERANCE = 1e-6
    DEFAULT_EXPECTED_TIME_SCALE = 1e4

    MAX_STAR_IMG_SCALE = 100000
    MIN_STAR_IMG_SCALE = 1
    MAX_PLANET_IMG_SCALE = 500000
    MIN_PLANET_IMG_SCALE = 1
    MAX_DISTANCE_SCALE = 1000
    MIN_DISTANCE_SCALE = 1
    MAX_NEW_STAR_MASS_SCALE = 1000
    MIN_NEW_STAR_MASS_SCALE = 1e-3
    MAX_DT = 100
    MIN_DT = 1e-10
    MAX_TIME_SPEED = 10000
    MIN_TIME_SPEED = 1
    MAX_RK_MAX_ITERATION = 50000
    MIN_RK_MAX_ITERATION = 1
    MAX_RK_MIN_ITERATION = 50000
    MIN_RK_MIN_ITERATION = 1
    MAX_TOLERANCE = 1e-4
    MIN_TOLERANCE = 1e-15
    MAX_EXPECTED_TIME_SCALE = 1e10
    MIN_EXPECTED_TIME_SCALE = 1

    DEFAULT_CHANGE_STAR_IMG_SCALE_SPEED = 1000
    DEFAULT_CHANGE_PLANET_IMG_SCALE_SPEED = 10000
    DEFAULT_CHANGE_DISTANCE_SCALE_SPEED = 10

    DEFAULT_NEW_STAR_VELOCITY_SCALE = 0.0001

    def __init__(
        self,
        screen_width: int,
        screen_height: int,
    ) -> None:
        # To change the default settings of screen_width and screen_height, go to _read_command_line_arg function in __main__.py
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.star_img_scale = self.DEFAULT_STAR_IMG_SCALE
        self.planet_img_scale = self.DEFAULT_PLANET_IMG_SCALE
        self.distance_scale = self.DEFAULT_DISTANCE_SCALE
        self.new_star_mass_scale = self.DEFAULT_NEW_STAR_MASS_SCALE
        self.dt = self.DEFAULT_dt
        self.time_speed = self.DEFAULT_TIME_SPEED
        # Initializing internal variable directly because rk_min_iteration.setter and rk_max_iteration.setter depends on each other
        self._rk_max_iteration = self.DEFAULT_RK_MAX_ITERATION
        self._rk_min_iteration = self.DEFAULT_RK_MIN_ITERATION
        self.tolerance = self.DEFAULT_TOLERANCE
        self.expected_time_scale = self.DEFAULT_EXPECTED_TIME_SCALE
        self.set_all_parameters_changing_false()
        self.current_changing_parameter = None
        self.is_hide_gui = False

    def scroll_change_parameters(self, magnitude):
        match self.current_changing_parameter:
            case "star_img_scale":
                for _ in range(abs(magnitude)):
                    if self.star_img_scale >= 10000:
                        if magnitude > 0:
                            self.star_img_scale += (
                                self.DEFAULT_CHANGE_STAR_IMG_SCALE_SPEED
                            )
                        elif magnitude < 0:
                            self.star_img_scale -= (
                                self.DEFAULT_CHANGE_STAR_IMG_SCALE_SPEED
                            )
                    else:
                        self.star_img_scale += self._rate_of_change(
                            self.star_img_scale, magnitude
                        )
            case "planet_img_scale":
                for _ in range(abs(magnitude)):
                    if self.planet_img_scale >= 100000:
                        if magnitude > 0:
                            self.planet_img_scale += (
                                self.DEFAULT_CHANGE_PLANET_IMG_SCALE_SPEED
                            )
                        elif magnitude < 0:
                            self.planet_img_scale -= (
                                self.DEFAULT_CHANGE_PLANET_IMG_SCALE_SPEED
                            )
                    else:
                        self.planet_img_scale += self._rate_of_change(
                            self.planet_img_scale, magnitude
                        )
            case "distance_scale":
                for _ in range(abs(magnitude)):
                    if self.distance_scale >= 100:
                        if magnitude > 0:
                            self.distance_scale += (
                                self.DEFAULT_CHANGE_DISTANCE_SCALE_SPEED
                            )
                        elif magnitude < 0:
                            self.distance_scale -= (
                                self.DEFAULT_CHANGE_DISTANCE_SCALE_SPEED
                            )
                    else:
                        self.distance_scale += self._rate_of_change(
                            self.distance_scale, magnitude
                        )
            case "new_star_mass_scale":
                for _ in range(abs(magnitude)):
                    self.new_star_mass_scale += self._rate_of_change(
                        self.new_star_mass_scale, magnitude
                    )
            case "dt":
                for _ in range(abs(magnitude)):
                    self.dt += self._rate_of_change(self.dt, magnitude)
            case "time_speed":
                for _ in range(abs(magnitude)):
                    self.time_speed += self._rate_of_change(self.time_speed, magnitude)
            case "rk_max_iteration":
                for _ in range(abs(magnitude)):
                    self.rk_max_iteration += self._rate_of_change(
                        self.rk_max_iteration, magnitude
                    )
            case "rk_min_iteration":
                for _ in range(abs(magnitude)):
                    self.rk_min_iteration += self._rate_of_change(
                        self.rk_min_iteration, magnitude
                    )
            case "tolerance":
                for _ in range(abs(magnitude)):
                    self.tolerance += self._rate_of_change(self.tolerance, magnitude)

    @staticmethod
    def _rate_of_change(x: float, magnitude: int) -> float:
        if magnitude > 0:
            return 10 ** math.floor(math.log10(x))
        elif magnitude < 0:
            if f"{x:e}"[0] == "1":
                return -x * 0.1
            else:
                return -(10 ** math.floor(math.log10(x)))

    def check_current_changing_parameter(self):
        if self.is_changing_star_img_scale == True:
            self.current_changing_parameter = "star_img_scale"
        elif self.is_changing_planet_img_scale == True:
            self.current_changing_parameter = "planet_img_scale"
        elif self.is_changing_distance_scale == True:
            self.current_changing_parameter = "distance_scale"
        elif self.is_changing_new_star_mass_scale == True:
            self.current_changing_parameter = "new_star_mass_scale"
        elif self.is_changing_dt == True:
            self.current_changing_parameter = "dt"
        elif self.is_changing_time_speed == True:
            self.current_changing_parameter = "time_speed"
        elif self.is_changing_rk_max_iteration == True:
            self.current_changing_parameter = "rk_max_iteration"
        elif self.is_changing_rk_min_iteration == True:
            self.current_changing_parameter = "rk_min_iteration"
        elif self.is_changing_tolerance == True:
            self.current_changing_parameter = "tolerance"

    def set_all_parameters_changing_false(self):
        self.is_changing_star_img_scale = False
        self.is_changing_planet_img_scale = False
        self.is_changing_distance_scale = False
        self.is_changing_new_star_mass_scale = False
        self.is_changing_dt = False
        self.is_changing_time_speed = False
        self.is_changing_rk_max_iteration = False
        self.is_changing_rk_min_iteration = False
        self.is_changing_tolerance = False

    def reset_parameters(self):
        self.star_img_scale = self.DEFAULT_STAR_IMG_SCALE
        self.planet_img_scale = self.DEFAULT_PLANET_IMG_SCALE
        self.time_speed = self.DEFAULT_TIME_SPEED
        self.new_star_mass_scale = self.DEFAULT_NEW_STAR_MASS_SCALE
        self.dt = self.DEFAULT_dt
        self.distance_scale = self.DEFAULT_DISTANCE_SCALE
        self.rk_max_iteration = self.DEFAULT_RK_MAX_ITERATION
        self.rk_min_iteration = self.DEFAULT_RK_MIN_ITERATION
        self.tolerance = self.DEFAULT_TOLERANCE

    @property
    def screen_width(self):
        return self._screen_width

    @screen_width.setter
    def screen_width(self, value):
        if value < 0:
            self._screen_width = 0
        else:
            self._screen_width = value

    @property
    def screen_height(self):
        return self._screen_height

    @screen_height.setter
    def screen_height(self, value):
        if value < 0:
            self._screen_height = 0
        else:
            self._screen_height = value

    @property
    def star_img_scale(self):
        return self._star_img_scale

    # Img may corrupt if the scale is too large.
    @star_img_scale.setter
    def star_img_scale(self, value):
        if value > self.MAX_STAR_IMG_SCALE:
            self._star_img_scale = self.MAX_STAR_IMG_SCALE
        elif value < self.MIN_STAR_IMG_SCALE:
            self._star_img_scale = self.MIN_STAR_IMG_SCALE
        else:
            self._star_img_scale = int(value)

    @property
    def planet_img_scale(self):
        return self._planet_img_scale

    # Img may corrupt if the scale is too large.
    @planet_img_scale.setter
    def planet_img_scale(self, value):
        if value > self.MAX_PLANET_IMG_SCALE:
            self._planet_img_scale = self.MAX_PLANET_IMG_SCALE
        elif value < self.MIN_PLANET_IMG_SCALE:
            self._planet_img_scale = self.MIN_PLANET_IMG_SCALE
        else:
            self._planet_img_scale = int(value)

    @property
    def distance_scale(self):
        return self._distance_scale

    @distance_scale.setter
    def distance_scale(self, value):
        if value > self.MAX_DISTANCE_SCALE:
            self._distance_scale = self.MAX_DISTANCE_SCALE
        elif value <= self.MIN_DISTANCE_SCALE:
            self._distance_scale = self.MIN_DISTANCE_SCALE
        else:
            self._distance_scale = int(value)

    @property
    def new_star_mass_scale(self):
        return self._new_star_mass_scale

    @new_star_mass_scale.setter
    def new_star_mass_scale(self, value):
        if value > self.MAX_NEW_STAR_MASS_SCALE:
            self._new_star_mass_scale = self.MAX_NEW_STAR_MASS_SCALE
        elif value < self.MIN_NEW_STAR_MASS_SCALE:
            self._new_star_mass_scale = self.MIN_NEW_STAR_MASS_SCALE
        else:
            self._new_star_mass_scale = round(value, ndigits=15)

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        if value > self.MAX_DT:
            self._dt = self.MAX_DT
        elif value < self.MIN_DT:
            self._dt = self.MIN_DT
        else:
            self._dt = round(value, ndigits=15)

    @property
    def time_speed(self):
        return self._time_speed

    @time_speed.setter
    def time_speed(self, value):
        if value > self.MAX_TIME_SPEED:
            self._time_speed = self.MAX_TIME_SPEED
        elif value < self.MIN_TIME_SPEED:
            self._time_speed = self.MIN_TIME_SPEED
        else:
            self._time_speed = int(value)

    @property
    def rk_max_iteration(self):
        return self._rk_max_iteration

    @rk_max_iteration.setter
    def rk_max_iteration(self, value):
        if value > self.MAX_RK_MAX_ITERATION:
            self._rk_max_iteration = self.MAX_RK_MAX_ITERATION
        elif value < self.MIN_RK_MAX_ITERATION:
            self._rk_max_iteration = self.MIN_RK_MAX_ITERATION
        elif value < self.rk_min_iteration:
            self._rk_max_iteration = self.rk_min_iteration
        else:
            self._rk_max_iteration = int(value)

    @property
    def rk_min_iteration(self):
        return self._rk_min_iteration

    @rk_min_iteration.setter
    def rk_min_iteration(self, value):
        if value > self.MAX_RK_MIN_ITERATION:
            self._rk_min_iteration = self.MAX_RK_MIN_ITERATION
        elif value > self.rk_max_iteration:
            self._rk_min_iteration = self.rk_max_iteration
        elif value < self.MIN_RK_MIN_ITERATION:
            self._rk_min_iteration = self.MIN_RK_MIN_ITERATION
        else:
            self._rk_min_iteration = int(value)

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value):
        if value > self.MAX_TOLERANCE:
            self._tolerance = self.MAX_TOLERANCE
        elif value < self.MIN_TOLERANCE:
            self._tolerance = self.MIN_TOLERANCE
        else:
            self._tolerance = round(value, ndigits=15)

    @property
    def expected_time_scale(self):
        return self._expected_time_scale

    @expected_time_scale.setter
    def expected_time_scale(self, value):
        if value > self.MAX_EXPECTED_TIME_SCALE:
            self._expected_time_scale = self.MAX_EXPECTED_TIME_SCALE
        elif value < self.MIN_EXPECTED_TIME_SCALE:
            self._expected_time_scale = self.MIN_EXPECTED_TIME_SCALE
        else:
            self._expected_time_scale = int(value)
