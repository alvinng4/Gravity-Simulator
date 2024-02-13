from pathlib import Path
import sys
path = str(Path(Path(Path(__file__).parent.absolute()).parent.absolute()).parent.absolute()) + "/gravity_sim/"
sys.path.insert(0, path)

from settings import Settings


def test_settings_screen_size():
    settings1 = Settings(screen_width=1920, screen_height=1080)
    assert settings1._screen_width == 1920
    assert settings1._screen_height == 1080
    settings2 = Settings(screen_width=200.5, screen_height=201.5)
    assert settings2._screen_width == 200.5
    assert settings2._screen_height == 201.5
    settings3 = Settings(screen_width=-100, screen_height=-101)
    assert settings3._screen_width == 0
    assert settings3._screen_height == 0
    settings4 = Settings(screen_width=0, screen_height=0)
    assert settings4._screen_width == 0
    assert settings4._screen_height == 0


def test_settings_star_img_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.star_img_scale == 5000
    settings.star_img_scale = 10
    assert settings.star_img_scale == 10
    settings.star_img_scale = 0
    assert settings.star_img_scale == 1
    settings.star_img_scale = -10
    assert settings.star_img_scale == 1
    settings.star_img_scale = 150000
    assert settings.star_img_scale == 100000


def test_settings_planet_img_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.planet_img_scale == 100000
    settings.planet_img_scale = 10
    assert settings.planet_img_scale == 10
    settings.planet_img_scale = 1
    assert settings.planet_img_scale == 1
    settings.planet_img_scale = -10
    assert settings.planet_img_scale == 1
    settings.planet_img_scale = 550000
    assert settings.planet_img_scale == 500000


def test_settings_distance_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.distance_scale == 200
    settings.distance_scale = 10
    assert settings.distance_scale == 10
    settings.distance_scale = 1000
    assert settings.distance_scale == 1000
    settings.distance_scale = -10.5
    assert settings.distance_scale == 1
    settings.distance_scale = 0
    assert settings.distance_scale == 1
    settings.distance_scale = 1500
    assert settings.distance_scale == 1000

def test_new_star_mass_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.new_star_mass_scale == 1
    settings.new_star_mass_scale = 10
    assert settings.new_star_mass_scale == 10
    settings.new_star_mass_scale = 1000
    assert settings.new_star_mass_scale == 1000
    settings.new_star_mass_scale = -10.5
    assert settings.new_star_mass_scale == 1e-3
    settings.new_star_mass_scale = 0
    assert settings.new_star_mass_scale == 1e-3
    settings.new_star_mass_scale = 1500
    assert settings.new_star_mass_scale == 1000

def test_settings_dt():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.dt == 0.1
    settings.dt = 10
    assert settings.dt == 10
    settings.dt = 0.00001
    assert settings.dt == 0.00001
    settings.dt = -10.5
    assert settings.dt == 1e-10
    settings.dt = 0
    assert settings.dt == 1e-10
    settings.dt = 150
    assert settings.dt == 100


def test_settings_time_speed():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.time_speed == 1
    settings.time_speed = 10
    assert settings.time_speed == 10
    settings.time_speed = 1000
    assert settings.time_speed == 1000
    settings.time_speed = -10.5
    assert settings.time_speed == 1
    settings.time_speed = 0
    assert settings.time_speed == 1
    settings.time_speed = 15000
    assert settings.time_speed == 10000


def test_epsilon():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.epsilon == 1e-2
    settings.epsilon = 1
    assert settings.epsilon == 1
    settings.epsilon = 10
    assert settings.epsilon == 10
    settings.epsilon = -10.5
    assert settings.epsilon == 1e-15
    settings.epsilon = 0
    assert settings.epsilon == 1e-15
    settings.epsilon = 150
    assert settings.epsilon == 100



