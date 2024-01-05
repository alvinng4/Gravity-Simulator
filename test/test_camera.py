from camera import Camera


def test_default():
    camera = Camera()
    assert camera._img_scale == 1


def test_camera_img_scale():
    camera1 = Camera(img_scale=10.1)
    assert camera1._img_scale == 10.1
    camera2 = Camera(img_scale=1000)
    assert camera2._img_scale == 100
    camera3 = Camera(img_scale=0.1)
    assert camera3._img_scale == 1
    camera3 = Camera(img_scale=-10)
    assert camera3._img_scale == 1


def test_camera_movement():
    camera = Camera()
    camera.moving_right = True
    camera.moving_up = True
    assert camera.pos_x == 0
    assert camera.pos_y == 0
    camera.update_movement()
    assert camera.pos_x == 10
    assert camera.pos_y == -10
    camera.moving_left = True
    camera.moving_down = True
    camera.update_movement()
    assert camera.pos_x == 10
    assert camera.pos_y == -10
    camera.moving_right = False
    camera.moving_up = False
    camera.update_movement()
    camera.update_movement()
    assert camera.pos_x == -10
    assert camera.pos_y == 10

