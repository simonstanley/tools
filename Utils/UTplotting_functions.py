import unittest
from plotting_functions import *

class Test_get_sensible_contour_levels(unittest.TestCase):
    def test_contour_levels(self):
        c_levels = get_sensible_contour_levels(0.2, 76.5, 5)
        self.assertEqual(list(c_levels), [0, 20, 40, 60, 80])
    
    def test_symmetry(self):
        c_levels = get_sensible_contour_levels(-3.23, 3.23, 10)
        self.assertEqual(c_levels[0], -c_levels[-1])
        

class Test_colour_map(unittest.TestCase):
    def test_cmap(self):
        cmap = colour_map(['blue', 'white', 'red'])
        self.assertEqual(cmap(0), (0.0, 0.0, 1.0, 1.0))
        self.assertEqual(cmap(256), (1.0, 0.0, 0.0, 1.0))
    
    def test_colour_postion(self):
        # Get the colour white to the middle value 2.5.
        data = [1,2,3,4]
        cmap = colour_map(['blue', 'white', 'pink', 'red'], 'white', 2.5, data=data)
        mid_colour = int((cmap.N / 2.) - 1)
        self.assertEqual(cmap(mid_colour), (1.0, 1.0, 1.0, 1.0))
        