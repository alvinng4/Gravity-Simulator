import time

import pygame

from text_box import Text_box


class Stats:
    """Track statistics for Gravity Simulator."""

    STATSBOARD_FONT_SIZE = 20
    STATSBOARD_SIZE_X = 260
    STATSBOARD_SIZE_Y = 23

    def __init__(self, grav_sim) -> None:
        self.simulation_time = 0
        self.fps = grav_sim.clock.get_fps()
        self.total_energy = 0
        self.settings = grav_sim.settings
        self.is_paused = False
        self.is_holding_rclick = False
        self.create_statsboard(grav_sim)
        self._statsboard_init_print_msg()

    def update(self, grav_sim) -> None:
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)

        if grav_sim.menu.main_menu_active == True:
            self.start_time = time.time()
        if self.is_paused == False:
            self.run_time = time.time() - self.start_time

        # Count time for the create new star function
        if self.is_holding_rclick == True:
            self.holding_rclick_time = time.time() - self.holding_rclick_start_time

    def reset(self, grav_sim) -> None:
        self.start_time = time.time()
        self.simulation_time = 0
        self.total_energy = 0
        grav_sim.simulator.is_initialize = True
        grav_sim.simulator.is_initialize_integrator = grav_sim.simulator.current_integrator
        grav_sim.camera._pos[0] = 0
        grav_sim.camera._pos[1] = 0

    def start_pause(self) -> None:
        self.paused_start_time = time.time()
        self.is_paused = True

    def end_pause(self) -> None:
        self.start_time -= self.paused_start_time - time.time()
        self.is_paused = False

    def start_holding_rclick(self) -> None:
        self.holding_rclick_start_time = time.time()
        self.is_holding_rclick = True

    def end_holding_rclick(self) -> None:
        self.is_holding_rclick = False

    def print_msg(self) -> None:
        self.fps_board.print_msg(f"FPS = {self.fps:2.1f}")
        self.obj_count_board.print_msg(f"Object = {self.objects_count}")
        self.simulation_time_board.print_msg(
            f"Simulation Time = {self.simulation_time / 365.2425:.1e} years"
        )
        self.run_time_board.print_msg(f"Run time = {self.run_time:.0f} seconds")
        self.total_energy_board.print_msg(f"Total Energy = {self.total_energy:.3e}")

        self.star_img_scale_board.print_msg(
            f"Star Image Scale = {self.settings.star_img_scale:d}"
        )
        self.planet_img_scale_board.print_msg(
            f"Planet Image Scale = {self.settings.planet_img_scale:d}"
        )
        self.distance_scale_board.print_msg(
            f"Distance Scale = {self.settings.distance_scale}"
        )
        self.new_star_mass_scale_board.print_msg(
            f"New star mass scale = {self.settings.new_star_mass_scale:g}x"
        )
        self.dt_board.print_msg(f"dt = {self.settings.dt:g} days / frame")
        self.time_speed_board.print_msg(f"Time Speed = {self.settings.time_speed:d}x")
        self.epsilon_board.print_msg(f"Epsilon = {self.settings.epsilon:g}")

    def draw(self, grav_sim) -> None:
        self.print_msg()
        self.fps_board.draw()
        self.obj_count_board.draw()
        self.simulation_time_board.draw()
        self.run_time_board.draw()
        self.total_energy_board.draw()

        self.parameters_board.draw()
        self.star_img_scale_board.draw()
        self.planet_img_scale_board.draw()
        self.distance_scale_board.draw()
        self.new_star_mass_scale_board.draw()
        self.dt_board.draw()
        self.time_speed_board.draw()
        self.epsilon_board.draw()

        self.integrators_board.draw()
        self.euler_board.draw()
        self.euler_cromer_board.draw()
        self.rk2_board.draw()
        self.rk4_board.draw()
        self.leapfrog_board.draw()
        self.rkf45_board.draw()

        # Visual indicator for currently changing parameter
        match self.settings.current_changing_parameter:
            case "star_img_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.star_img_scale_board.rect.centery + 5),
                    4,
                )
            case "planet_img_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.planet_img_scale_board.rect.centery + 5),
                    4,
                )
            case "distance_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.distance_scale_board.rect.centery + 5),
                    4,
                )
            case "new_star_mass_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.new_star_mass_scale_board.rect.centery + 5),
                    4,
                )
            case "dt":
                pygame.draw.circle(
                    grav_sim.screen, "yellow", (290, self.dt_board.rect.centery + 5), 4
                )
            case "time_speed":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.time_speed_board.rect.centery + 5),
                    4,
                )
            case "epsilon":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.epsilon_board.rect.centery + 5),
                    4,
                )

        # Visual indicator for currently selected integrator
        match grav_sim.simulator.current_integrator:
            case "euler":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (250, self.euler_board.rect.centery + 5),
                    4,
                )
            case "euler_cromer":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (250, self.euler_cromer_board.rect.centery + 5),
                    4,
                )
            case "rk2":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.rk2_board.rect.centery + 5), 4
                )
            case "rk4":
                pygame.draw.circle(
                    grav_sim.screen, "green", (250, self.rk4_board.rect.centery + 5), 4
                )
            case "leapfrog":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (250, self.leapfrog_board.rect.centery + 5),
                    4,
                )
            case "rkf45":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (250, self.rkf45_board.rect.centery + 5),
                    4,
                )

    def check_button(self, grav_sim, mouse_pos) -> None:
        """Check if there is any click on the buttons"""
        if self.settings.is_hide_gui == False:
            if self.star_img_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_star_img_scale = True
            if self.planet_img_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_planet_img_scale = True
            if self.distance_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_distance_scale = True
            if self.new_star_mass_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_new_star_mass_scale = True
            if self.dt_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_dt = True
            if self.time_speed_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_time_speed = True
            if self.epsilon_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_epsilon = True

            if self.euler_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_euler = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "euler"
            if self.euler_cromer_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_euler_cromer = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "euler_cromer"
            if self.rk2_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rk2 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rk2"
            if self.rk4_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rk4 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rk4"
            if self.leapfrog_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_leapfrog = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "leapfrog"
            if self.rkf45_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rkf45 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rkf45"

    def _statsboard_init_print_msg(self) -> None:
        self.parameters_board.print_msg("Parameters: (Click to select)")
        self.integrators_board.print_msg(f"Integrators: (Click to select)")
        self.euler_board.print_msg(f"Euler")
        self.euler_cromer_board.print_msg(f"Euler-Cromer")
        self.rk2_board.print_msg(f"2nd order Runge-Kutta")
        self.rk4_board.print_msg(f"4th order Runge-Kutta")
        self.leapfrog_board.print_msg(f"Leapfrog (Verlet)")
        self.rkf45_board.print_msg(f"RKF4(5)")

    @classmethod
    def create_statsboard(self, grav_sim) -> None:
        self.fps_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 0),
        )
        self.obj_count_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 23),
        )

        self.simulation_time_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 46),
        )
        self.run_time_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 69),
        )
        self.total_energy_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 92),
        )

        self.parameters_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 138),
        )
        self.star_img_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 161),
        )
        self.planet_img_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 184),
        )
        self.distance_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 207),
        )
        self.new_star_mass_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 230),
        )
        self.dt_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 253),
        )
        self.time_speed_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 276),
        )
        self.epsilon_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 299),
        )

        self.integrators_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 345),
        )
        self.euler_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 368),
        )
        self.euler_cromer_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 391),
        )
        self.rk2_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 414),
        )
        self.rk4_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 437),
        )
        self.leapfrog_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 460),
        )
        self.rkf45_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 506),
        )
