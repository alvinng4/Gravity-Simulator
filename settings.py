class Settings:
    """A class to store all settings for gravity simulator."""

    def __init__(self, img_scale=50):
        self.screen_width = 1366
        self.screen_height = 768
        self.bg_color = (0, 0, 0)  # Background color

        self.fps = 60

        self.grav_obj_color = (255, 255, 255)

        self.img_scale = img_scale


    @property
    def img_scale(self):
        return self._img_scale

    @img_scale.setter
    def img_scale(self, value):
        if value > 100:
            self._img_scale = 100
        elif value < 1:
            self._img_scale = 1
        else:
            self._img_scale = value