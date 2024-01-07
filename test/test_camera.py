from gravity_sim.camera import Camera


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

