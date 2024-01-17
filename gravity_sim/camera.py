class Camera:
    def __init__(self):
        self._pos = [0, 0]
        self.speed = [10, 10]

        # Movement flag
        self.moving_right = False
        self.moving_left = False
        self.moving_up = False
        self.moving_down = False

    @property
    def pos(self):
        return tuple(self._pos)

    def update_movement(self):
        if self.moving_right == True:
            self._pos[0] += self.speed[0]
        if self.moving_left == True:
            self._pos[0] -= self.speed[0]
        if self.moving_up == True:
            self._pos[1] -= self.speed[1]
        if self.moving_down == True:
            self._pos[1] += self.speed[1]
