class Camera:
    def __init__(self, img_scale=1):
        self.pos_x = 0
        self.pos_y = 0
        self.speed_x = 10
        self.speed_y = 10

        self._img_scale = img_scale

        # Movement flag
        self.moving_right = False
        self.moving_left = False
        self.moving_up = False
        self.moving_down = False



    def pos(self):
        return (self.pos_x, self.pos_y)

    def speed_x(self):
        return self.speed_x

    def speed_y(self):
        return self.speed_y
    
    def update_movement(self):
        if self.moving_right == True:
            self.pos_x += self.speed_x
        if self.moving_left == True:
            self.pos_x -= self.speed_x
        if self.moving_up == True:
            self.pos_y -= self.speed_y
        if self.moving_down == True:
            self.pos_y += self.speed_y

    @property
    def img_scale(self):
        return self._img_scale   
    
    @img_scale.setter
    def img_scale(self, value):
        if self._img_scale > 100:
            self._img_scale = 100
        elif self._img_scale < 1:
            self._img_scale = 1
        else:
            self._img_scale = value


        
