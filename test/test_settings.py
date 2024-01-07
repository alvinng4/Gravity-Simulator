from gravity_sim.settings import Settings


def test_settings_screen_width():
    settings1 = Settings(sun_img_scale=10, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings1._screen_width == 1920
    settings2 = Settings(sun_img_scale=10, img_scale=10, screen_width=200.5, screen_height=1080)
    assert settings2._screen_width == 200.5
    settings3 = Settings(sun_img_scale=10, img_scale=10, screen_width=-100, screen_height=1080)
    assert settings3._screen_width == 0
    settings4 = Settings(sun_img_scale=10, img_scale=10, screen_width=0, screen_height=1080)
    assert settings4._screen_width == 0

def test_settings_screen_height():
    settings1 = Settings(sun_img_scale=10, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings1._screen_height == 1080
    settings2 = Settings(sun_img_scale=10, img_scale=10, screen_width=1920, screen_height=200.5)
    assert settings2._screen_height == 200.5
    settings3 = Settings(sun_img_scale=10, img_scale=10, screen_width=1920, screen_height=-100)
    assert settings3._screen_height == 0
    settings4 = Settings(sun_img_scale=10, img_scale=10, screen_width=1920, screen_height=0)
    assert settings4._screen_height == 0

def test_settings_sun_img_scale():
    settings1 = Settings(sun_img_scale=10.1, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings1._sun_img_scale == 10.1
    settings2 = Settings(sun_img_scale=1000, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings2._sun_img_scale == 100
    settings3 = Settings(sun_img_scale=0.1, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings3._sun_img_scale == 0.1
    settings4 = Settings(sun_img_scale=-10, img_scale=10, screen_width=1920, screen_height=1080)
    assert settings4._sun_img_scale == 0

def test_settings_img_scale():
    settings1 = Settings(sun_img_scale=10, img_scale=10.1, screen_width=1920, screen_height=1080)
    assert settings1._img_scale == 10.1
    settings2 = Settings(sun_img_scale=10, img_scale=1000, screen_width=1920, screen_height=1080)
    assert settings2._img_scale == 100
    settings3 = Settings(sun_img_scale=10, img_scale=0.1, screen_width=1920, screen_height=1080)
    assert settings3._img_scale == 0.1
    settings4 = Settings(sun_img_scale=-10, img_scale=-10, screen_width=1920, screen_height=1080)
    assert settings4._img_scale == 0
