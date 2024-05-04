import numpy as np
from enum import Enum
import collections

from Utils import lerp

R0_km = 6371

class NeType(Enum):
    LINEAR = 1,
    PARABOLA = 2,
    CHAPMAN = 3

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class NeProfile:
    def __init__(self):
        self.points = []

        self.type = NeType.LINEAR

        self.fc_MHz = 6
        self.h0_km = 100
        self.hm_km = 300
        self.ht_km = 500


    def linear(self, h0, hm, h_step_km = 1, k = (1/R0_km), b = 1):
        heights = np.arange(h0, hm, h_step_km)
        self.points = list(map(lambda z: Point(k*z + b, z), heights))
        self.hm_km = hm
        self.h0_km = h0
        self.ht_km = h0 + 2*(hm-h0)
        return self

    def parabola(self):
        return self

    def get_n(self, z) -> float:
        if z < self.h0_km or z > self.hm_km:
            return 0

        if self.type == NeType.LINEAR:
            closest_idx = min(range(len(self.points)), key=lambda i: abs(self.points[i].y - z))
            closest_point = self.points[closest_idx]

            if closest_point.y == z:
                return closest_point.x
            elif closest_point.y < z:
                sp = self.points[min(len(self.points) - 1, closest_idx + 1)]
            else:
                sp = self.points[max(0, closest_idx - 1)]

            return lerp(closest_point.y, closest_point.x, sp.y, sp.x, z)
