import numpy as np
from gelas import Gelas, Fallen
# Constants
K = 0.53 # Stability Coef

class Map:
    def __init__(self):
        self.list_gelas = []
        self.index_lantai = []
        self.interactions = []
        self.has_fallen = False
  
    def add_gelas(self, x_coordinate, y_coordinate, lantai = False, tahan_kiri = False, tahan_kanan = False):
        self.list_gelas.append(Gelas(x_coordinate, y_coordinate, tahan_kiri=tahan_kiri, tahan_kanan = tahan_kanan))
        if lantai:
            self.index_lantai.append(len(self.list_gelas)-1)
  
    def interact(self, gelas_atas, gelas_bawah):
        F, center_coords = gelas_atas.total_weight()
        gelas_bawah.press(F,center_coords[0])
        # Kemiringan atas mengikuti kemiringan bawah
        gelas_atas.global_tilt = gelas_bawah.local_tilt
        # Posisi terbawah gelas atas bertemu dengan garis atas gelas bawah.
        vertex_ = gelas_atas.get_vertex()
        line_function = gelas_bawah.get_upper_line()
        delta_y = vertex_[1] - line_function(vertex_[0])
        gelas_atas.y -= delta_y

    def get_total_center_of_mass(self):
        sum_x = 0
        sum_y = 0
        sum_m = 0
        for gelas in self.list_gelas:
            coords = gelas.get_center_of_mass()
            sum_x += coords[0]
            sum_y += coords[1]
            sum_m += 1
        return sum_x/sum_m, sum_y/sum_m

    def get_shapes(self):
        return [g.get_shape() for g in self.list_gelas]

    def add_interaction(self, atas, bawah):
        self.interactions.append((atas, bawah))

    def interact_all(self):
        for i in range(len(self.interactions)):
            for a, b in self.interactions[:-i]:
                try:
                    self.interact(self.list_gelas[a], self.list_gelas[b])
                except Fallen:
                    self.has_fallen = True

    def get_interval(self):
        max_x = -999999999999
        min_x = 999999999999
        for i in self.index_lantai:
            c_x, _ = self.list_gelas[i].get_vertex()
            max_x = max(max_x, c_x)
            min_x = min(min_x, c_x)
        return min_x, max_x

    def is_stable(self):
        if self.has_fallen:
            return False
        min_x, max_x = self.get_interval()
        c_x, c_y = self.get_total_center_of_mass()
        epsilon = K / c_y
        return min_x - epsilon < c_x < max_x + epsilon