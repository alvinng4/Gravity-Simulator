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
            self.start_time = pygame.time.get_ticks()
        self.run_time = (pygame.time.get_ticks() - self.start_time) / 1000
        self.fps_board.print_msg(f"FPS = {round(self.fps, 1)}")
        self.obj_board.print_msg(f"Object = {self.objects_count}")
        self.distance_scale_board.print_msg(f"Distance Scale = {self.distance_scale}")
        self.dt_board.print_msg(f"dt = {self.dt} days / frame")
        self.time_speed_board.print_msg(f"Time Speed = {self.time_speed}x")
        self.simulation_time_board.print_msg(f"Simulation Time = {self.simulation_time / 365.2425:.1e} years")
        self.run_time_board.print_msg(f"Run time = {int(self.run_time)} seconds")

    def reset_stats(self):
        self.start_time = pygame.time.get_ticks()
        self.simulation_time = 0

    def draw(self):
        self.fps_board.draw()
        self.obj_board.draw()
        self.sun_img_scale_board.draw()
        self.img_scale_board.draw()
        self.distance_scale_board.draw()
        self.dt_board.draw()
        self.time_speed_board.draw()
        self.simulation_time_board.draw()
        self.run_time_board.draw()

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
