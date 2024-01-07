from gravity_sim.settings import Settings


def test_settings_sun_img_scale():
    settings1 = Settings(sun_img_scale=10.1)
    assert settings1._sun_img_scale == 10.1
    settings2 = Settings(sun_img_scale=1000)
    assert settings2._sun_img_scale == 100
    settings3 = Settings(sun_img_scale=0.1)
    assert settings3._sun_img_scale == 0.1
    settings4 = Settings(sun_img_scale=-10)
    assert settings4._sun_img_scale == 0

def test_settings_img_scale():
    settings1 = Settings(img_scale=10.1)
    assert settings1._img_scale == 10.1
    settings2 = Settings(img_scale=1000)
    assert settings2._img_scale == 100
    settings3 = Settings(img_scale=0.1)
    assert settings3._img_scale == 0.1
    settings4 = Settings(img_scale=-10)
    assert settings4._img_scale == 0

