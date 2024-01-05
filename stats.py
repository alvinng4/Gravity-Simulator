import pygame


class Stats:
    """Track statistics for Gravity Simulator."""

    def __init__(self):
        self.reset_stats()


    def update(self, grav_sim):
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)




    def reset_stats(self):
        pass