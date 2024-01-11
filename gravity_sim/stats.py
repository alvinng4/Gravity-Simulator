import time

import pygame

from text_box import Text_box


class Stats:
    """Track statistics for Gravity Simulator."""

    def __init__(self, grav_sim):
        self.simulation_time = 0
        self.fps = grav_sim.clock.get_fps()
        self.settings = grav_sim.settings
        self.sun_img_scale = grav_sim.settings.sun_img_scale
        self.img_scale = grav_sim.settings.img_scale
        self.total_energy = 0
        self.is_paused = False
        self.create_stats_board(grav_sim)
        # Image scale is fixed in the settings
        self.sun_img_scale_board.print_msg(f"Sun Image Scale = {self.sun_img_scale}")
        self.img_scale_board.print_msg(f"Objects Image Scale = {self.img_scale}")

    def update(self, grav_sim):
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)
        self.distance_scale = round(grav_sim.settings.distance_scale, 1)
        self.dt = grav_sim.settings.dt
        self.time_speed = grav_sim.settings.time_speed
        if grav_sim.menu.main_menu_active == True:
            self.start_time = time.time()

        self.run_time = time.time() - self.start_time

    def start_pause(self):
        self.paused_start_time = time.time()
        self.is_paused = True

    def end_pause(self):
        self.start_time -= self.paused_start_time - time.time()
        self.is_paused = False

    def reset_stats(self):
        self.start_time = time.time()
        self.simulation_time = 0

    def print_msg(self):
        self.fps_board.print_msg(f"FPS = {round(self.fps, 1)}")
        self.obj_board.print_msg(f"Object = {self.objects_count}")
        self.distance_scale_board.print_msg(f"Distance Scale = {self.distance_scale}")
        self.dt_board.print_msg(f"dt = {self.dt} days / frame")
        self.time_speed_board.print_msg(f"Time Speed = {self.time_speed}x")
        self.simulation_time_board.print_msg(
            f"Simulation Time = {self.simulation_time / 365.2425:.1e} years"
        )
        self.run_time_board.print_msg(f"Run time = {int(self.run_time)} seconds")
        self.total_energy_board.print_msg(f"Total Energy = {self.total_energy:.3e}")

    def draw(self):
        self.print_msg()
        self.fps_board.draw()
        self.obj_board.draw()
        self.sun_img_scale_board.draw()
        self.img_scale_board.draw()
        self.distance_scale_board.draw()
        self.dt_board.draw()
        self.time_speed_board.draw()
        self.simulation_time_board.draw()
        self.run_time_board.draw()
        self.total_energy_board.draw()

    @classmethod
    def create_stats_board(self, grav_sim):
        self.fps_board = Text_box(
            grav_sim,
            0,
            0,
            20,
            font="Avenir",
            text_box_left_top=(10, 0),
        )
        self.obj_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 23),
        )
        self.sun_img_scale_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 46),
        )
        self.img_scale_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 69),
        )
        self.distance_scale_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 92),
        )
        self.dt_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 115),
        )
        self.time_speed_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 138),
        )
        self.simulation_time_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 161),
        )
        self.run_time_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 184),
        )
        self.total_energy_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 207),
        )

