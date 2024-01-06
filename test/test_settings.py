from settings import Settings


def test_settings_img_scale():
    settings1 = Settings(img_scale=10.1)
    assert settings1._img_scale == 10.1
    settings2 = Settings(img_scale=1000)
    assert settings2._img_scale == 100
    settings3 = Settings(img_scale=0.1)
    assert settings3._img_scale == 1
    settings3 = Settings(img_scale=-10)
    assert settings3._img_scale == 1

