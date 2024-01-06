import pygame

from text_box import Text_box


class Stats:
    """Track statistics for Gravity Simulator."""

    def __init__(self, grav_sim):
        self.fps = grav_sim.clock.get_fps()
        self.start_time = pygame.time.get_ticks()
        self.settings = grav_sim.settings

        self.create_stats_board(grav_sim)

    def update(self, grav_sim):
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)
        self.sun_img_scale = self.settings.sun_img_scale
        self.img_scale = self.settings.img_scale

        self.fps_board.print_msg(f"FPS = {round(self.fps, 1)}")
        self.obj_board.print_msg(f"Object = {self.objects_count}")
        self.sun_img_scale_board.print_msg(f"Sun Image Scale = {self.sun_img_scale}")
        self.img_scale_board.print_msg(f"Objects Image Scale = {self.img_scale}")

        self.fps_board.draw()
        self.obj_board.draw()
        self.sun_img_scale_board.draw()
        self.img_scale_board.draw()

    def reset_stats(self):
        self.start_time = pygame.time.get_ticks()

    @classmethod
    def create_stats_board(self, grav_sim):
        self.fps_board = Text_box(
            grav_sim,
            0.08,
            0.03,
            20,
            font="Avenir",
            text_box_top=0,
            text_box_left=0,
        )
        self.obj_board = Text_box(
            grav_sim,
            0.08,
            0.03,
            20,
            font="Avenir",
            text_box_top=self.fps_board.rect.bottom,
            text_box_left=0,
        )
        self.sun_img_scale_board = Text_box(
            grav_sim,
            0.145,
            0.03,
            20,
            font="Avenir",
            text_box_top=self.obj_board.rect.bottom,
            text_box_left=0,
        )
        self.img_scale_board = Text_box(
            grav_sim,
            0.18,
            0.03,
            20,
            font="Avenir",
            text_box_top=self.sun_img_scale_board.rect.bottom,
            text_box_left=0,
        )
