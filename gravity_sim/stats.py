import time

import pygame

from text_box import Text_box


class Stats:
    """Track statistics for Gravity Simulator."""

    def __init__(self, grav_sim):
        self.simulation_time = 0
        self.fps = grav_sim.clock.get_fps()
        self.settings = grav_sim.settings
        self.star_img_scale = grav_sim.settings.star_img_scale
        self.img_scale = grav_sim.settings.img_scale
        self.total_energy = 0
        self.is_paused = False
        self.is_holding_rclick = False
        self.create_stats_board(grav_sim)
        # Image scale is fixed in the settings, only need to be drawn once (Change later)
        self.star_img_scale_board.print_msg(f"Star Image Scale = {self.star_img_scale}")
        self.img_scale_board.print_msg(f"Planet Image Scale = {self.img_scale}")
        self.integrators_board.print_msg(f"Integrators: (Click to switch)")
        self.euler_board.print_msg(f"Euler")
        self.euler_cromer_board.print_msg(f"Euler-Cromer")
        self.rk2_board.print_msg(f"2nd order Runge-Kutta")
        self.rk4_board.print_msg(f"4th order Runge-Kutta")
        self.leapfrog_board.print_msg(f"Leapfrog (Verlet)")

    def update(self, grav_sim):
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)
        self.distance_scale = round(grav_sim.settings.distance_scale, 1)
        self.dt = grav_sim.settings.dt
        self.time_speed = grav_sim.settings.time_speed

        if grav_sim.menu.main_menu_active == True:
            self.start_time = time.time()

        # Pause time counting when paused
        if self.is_paused == False:
            self.run_time = time.time() - self.start_time

        if self.is_holding_rclick == True:
            self.holding_rclick_time = time.time() - self.holding_rclick_start_time

    def reset(self, grav_sim):
        self.start_time = time.time()
        self.simulation_time = 0
        grav_sim.simulator.is_initialize = True
        grav_sim.camera.pos_x = 0
        grav_sim.camera.pos_y = 0

    def start_pause(self):
        self.paused_start_time = time.time()
        self.is_paused = True

    def end_pause(self):
        self.start_time -= self.paused_start_time - time.time()
        self.is_paused = False

    def start_holding_rclick(self):
        self.holding_rclick_start_time = time.time()
        self.is_holding_rclick = True

    def end_holding_rclick(self):
        self.is_holding_rclick = False

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

    def draw(self, grav_sim):
        self.print_msg()
        self.fps_board.draw()
        self.obj_board.draw()
        self.star_img_scale_board.draw()
        self.img_scale_board.draw()
        self.distance_scale_board.draw()
        self.dt_board.draw()
        self.time_speed_board.draw()
        self.simulation_time_board.draw()
        self.run_time_board.draw()
        self.total_energy_board.draw()

        self.integrators_board.draw()
        self.euler_board.draw()
        self.euler_cromer_board.draw()
        self.rk2_board.draw()
        self.rk4_board.draw()
        self.leapfrog_board.draw()

        match grav_sim.simulator.current_integrator:
            case "euler":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.euler_board.rect.centery), 4
                )
            case "euler_cromer":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (250, self.euler_cromer_board.rect.centery),
                    4,
                )
            case "rk2":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.rk2_board.rect.centery), 4
                )
            case "rk4":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.rk4_board.rect.centery), 4
                )
            case "leapfrog":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.leapfrog_board.rect.centery), 4
                )

    def check_button(self, grav_sim, mouse_pos):
        """Check if there is any click on the buttons"""
        if self.euler_board.rect.collidepoint(mouse_pos):
            grav_sim.simulator.set_all_integrators_false()
            grav_sim.simulator.is_euler = True
        if self.euler_cromer_board.rect.collidepoint(mouse_pos):
            grav_sim.simulator.set_all_integrators_false()
            grav_sim.simulator.is_euler_cromer = True
        if self.rk2_board.rect.collidepoint(mouse_pos):
            grav_sim.simulator.set_all_integrators_false()
            grav_sim.simulator.is_rk2 = True
        if self.rk4_board.rect.collidepoint(mouse_pos):
            grav_sim.simulator.set_all_integrators_false()
            grav_sim.simulator.is_rk4 = True
        if self.leapfrog_board.rect.collidepoint(mouse_pos):
            grav_sim.simulator.set_all_integrators_false()
            grav_sim.simulator.is_leapfrog = True

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
        self.star_img_scale_board = Text_box(
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
        self.integrators_board = Text_box(
            grav_sim,
            0,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 253),
        )
        self.euler_board = Text_box(
            grav_sim,
            0.13,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 276),
        )
        self.euler_cromer_board = Text_box(
            grav_sim,
            0.13,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 299),
        )
        self.rk2_board = Text_box(
            grav_sim,
            0.13,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 322),
        )
        self.rk4_board = Text_box(
            grav_sim,
            0.13,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 345),
        )
        self.leapfrog_board = Text_box(
            grav_sim,
            0.13,
            0.03,
            20,
            font="Avenir",
            text_box_left_top=(10, 368),
        )
