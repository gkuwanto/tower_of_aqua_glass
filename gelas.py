import numpy as np
from matplotlib.path import Path

#Constants
W = 227 * 9.8 / 1000 # weight of 1 gelas (Newton)
MU = 0.035 # Friction Coef
MAX_theta = np.arctan(0.2/3.75)

def theta(force, x):
    # Calculate Tilt Physics
    if force == 0 or x == 0:
        return 0
    return np.arctan(1/(MU * ((force+W)/force) * np.sqrt((x/100)**2+(9.7/100)**2) / (x/100)))

class Fallen(Exception):
    """Error when the gelas has fallen"""
    pass


class Gelas:
    def __init__(self, x_global, y_global, tahan_kiri = False, tahan_kanan = False):
        # Global Coordinate
        self.x = x_global
        self.y = y_global
        self.global_tilt = 0
        self.F = 0
        self.local_tilt = 0
        self.tahan_kiri = tahan_kiri
        self.tahan_kanan = tahan_kanan
        # Local Coordinate
        self.axis_x = 0
        self.axis_y = 9.7 
        self.top_right_bound = (2.25,9.7)
        self.top_left_bound = (-2.25, 9.7)
        self.bot_right_bound = (3.75, 0.2)
        self.bot_left_bound = (-3.75, 0.2)
        self.bottom_vertex = (0,0)
        self.center = (0, 4.47)

    def press(self, force, x_global):
        # Update Coordinates and Tilt acording to Force acted
        self.F = force
        x = x_global - self.x
        x = (x if x<self.top_right_bound[0] else self.top_right_bound[0]) if x > self.top_left_bound[0] else self.top_left_bound[0]
        self.local_tilt = theta(force, x)
        if self.local_tilt > MAX_theta and not self.tahan_kanan:
            raise Fallen
        
        if self.local_tilt < MAX_theta and not self.tahan_kiri:
            raise Fallen

        self.local_tilt = max(-MAX_theta, min(self.local_tilt, MAX_theta))
        self.bottom_vertex = (x, 0)
        self.tilt()

    def tilt(self):
        # Update Coordinates according to local tilt
        top_delta_x = np.cos(self.local_tilt) * 2.25
        top_delta_y = np.sin(self.local_tilt) * 2.25
        bot_delta_x = np.cos(self.local_tilt) * 3.75
        bot_delta_y = np.sin(self.local_tilt) * 3.75

        self.top_right_bound = (self.axis_x + top_delta_x, self.axis_y - top_delta_y)
        self.top_left_bound = (self.axis_x - top_delta_x, self.axis_y + top_delta_y)
        self.bot_right_bound = (self.axis_x + bot_delta_x, self.axis_y - 9.5/np.cos(self.local_tilt) - bot_delta_y)
        self.bot_left_bound = (self.axis_x - bot_delta_x, self.axis_y - 9.5/np.cos(self.local_tilt) + bot_delta_y)
        self.update_center()
    
    def quadratic_solve(self):
        x1, y1 = self.bottom_vertex
        x2, y2 = self.bot_right_bound
        x3, y3 = self.bot_left_bound
        denom = (x1-x2) * (x1-x3) * (x2-x3)
        A = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
        B = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
        C = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom
        return lambda x: A*x**2 + B*x+ C

    def is_in_body(self, x,y):
        path = Path(np.array([self.top_left_bound, self.top_right_bound, self.bot_right_bound, self.bot_left_bound]))
        in_trapezium = path.contains_points([[x,y]])[0]
        g = lambda x: np.tan(self.local_tilt) * x + 9.7 - 9.5/ np.cos(self.local_tilt)
        f = self.quadratic_solve()
        in_parabola =  self.bot_left_bound[0] < x < self.bot_right_bound[0] and f(x)<y<g(x)
        return in_trapezium or in_parabola

    def update_center(self):
        shape = self.get_shape()
        x_sum = 0
        y_sum = 0
        m_sum = 0
        for i in range(len(shape)-1):
            m_sum += shape[i][0]*shape[i+1][1]-shape[i+1][0]*shape[i][1]
            x_sum += (shape[i][0]+shape[i+1][0]) * (shape[i][0]*shape[i+1][1]-shape[i+1][0]*shape[i][1])
            y_sum += (shape[i][1]+shape[i+1][1]) * (shape[i][0]*shape[i+1][1]-shape[i+1][0]*shape[i][1])
        A = m_sum / 2
        self.center = (x_sum / (6*A), y_sum / (6*A))

    def get_global_coord(self, local_coords):
        # Inverse of Rotation by global tilt composed with translation of global coords
        return (np.linalg.inv(
            np.array([[np.cos(self.global_tilt), -np.sin(self.global_tilt)], 
                    [np.sin(self.global_tilt), np.cos(self.global_tilt)]])
            ) @ np.array(local_coords) 
            )+ np.array((self.x, self.y))
        
    def total_weight(self):
        # Get Force and Position
        return W + self.F, self.get_center_of_mass()

    def get_center_of_mass(self):
        # Get global coordinates of center of mass
        return self.get_global_coord(self.center)
    def get_vertex(self):
        # Get global coordinates of bottom vertex
        return self.get_global_coord(self.bottom_vertex)
    
    def get_upper_line(self):
        # Get Line equation y = mx + c
        x1,y1 = self.get_global_coord(self.top_right_bound)
        x2,y2 = self.get_global_coord(self.top_left_bound)
        return lambda x: ((x-x1)*(y2-y1))/(x2-x1) + y1

    def get_shape(self):
        # Get list of coords pertaining to shape
        coords = [
            self.get_global_coord(self.top_left_bound), 
            self.get_global_coord(self.top_right_bound), 
            self.get_global_coord(self.bot_right_bound),
        ]
        f = self.quadratic_solve()
        for x in np.arange(self.bot_right_bound[0], self.bot_left_bound[0], -0.01):
            coords.append(self.get_global_coord((x,f(x))))
        coords.append(self.get_global_coord(self.bot_left_bound))
        coords.append(self.get_global_coord(self.top_left_bound))
        return coords