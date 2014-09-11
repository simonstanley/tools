import unittest
from iris_functions import *

class Test_cube_time_converter(unittest.TestCase):
    def setUp(self):
        self.datetime  = datetime.datetime(2000,2,3)
        self.time_unit = iris.unit.Unit('hours since 1970-01-01 00:00:00')
        self.cube_time = 263760.
        
    def test_convert_to_cube_time(self):
        cube_time = cube_time_converter(self.datetime, self.time_unit)
        self.assertEqual(cube_time, self.cube_time)
        
    def test_convert_from_cube_time(self):
        datetime = cube_time_converter(self.cube_time, self.time_unit)
        self.assertEqual(datetime, self.datetime)
        

class Test_get_time_bound_constraints(unittest.TestCase):
    def setUp(self):
        self.time_start = 35000.
        self.time_end   = 35072.
        
    def test_dates_in_wrong_order(self):
        self.assertRaises(AssertionError, get_time_bound_constraints,
                          self.time_end, self.time_start)


class Test_get_list_of_time_bound_constraints(unittest.TestCase):
    def setUp(self):
        self.times = [[35000,35023],[35050,35061],[36000,36020]]
        
    def test_return_value(self):
        self.assertEqual(type(get_list_of_time_bound_constraints(self.times)),
                         iris.Constraint)
        

class Test_get_time_bounds_from_cube(unittest.TestCase):
    def setUp(self):
        self.cube = iris.load(iris.sample_data_path('uk_hires.pp'))[0]
        
    def test_time_bounds(self):
        bounds = get_time_bounds_from_cube(self.cube)
        self.assertSequenceEqual(bounds, [datetime.datetime(2009,11,19,10),
                                          datetime.datetime(2009,11,19,12)])
    
    def test_return_none_for_bad_cube(self):
        self.cube = iris.cube.Cube([0])
        bounds = get_time_bounds_from_cube(self.cube)
        self.assertEqual(bounds, None)
        

class Test_get_coordinate_slice_dimensions(unittest.TestCase):
    def setUp(self):
        self.cube = iris.load(iris.sample_data_path('uk_hires.pp'))[0]
        
    def test_dim_coords(self):
        coords = get_coordinate_slice_dimensions(self.cube, 
                                                 ['time', 
                                                  'model_level_number'])
        self.assertEqual(coords, ['grid_latitude', 'grid_longitude'])
        
    def test_aux_coord(self):
        coords = get_coordinate_slice_dimensions(self.cube, 'forecast_period')
        self.assertEqual(coords, ['model_level_number', 'grid_latitude', 
                                  'grid_longitude'])
        
    def test_bad_coord(self):
        self.assertRaises(ValueError, get_coordinate_slice_dimensions,
                                      self.cube, 'no_coord')
    
    
class Test_get_cubelist_data_range(unittest.TestCase):
    def setUp(self):
        cubelist = []
        for i in range(3,10):
            cubelist.append(iris.cube.Cube([i]))
        self.cubelist = iris.cube.CubeList(cubelist)
        self.min = 3
        self.max = 9
        
    def test_data_range(self):
        data_range = get_cubelist_data_range(self.cubelist)
        self.assertEqual(data_range, (self.min, self.max))
    

class Test_get_xy_coords(unittest.TestCase):
    def setUp(self):
        self.cube = iris.load(iris.sample_data_path('uk_hires.pp'))[0]
    
    def test_xy_coords(self):
        xy_coords = get_xy_coords(self.cube)
        self.assertEqual(xy_coords, [self.cube.coord('grid_longitude'),
                                     self.cube.coord('grid_latitude')])


class Test_order_cubelist(unittest.TestCase):
    def setUp(self):
        cubelist = []
        for i, t in zip(range(3,8), [4,1,2,5,3]):
            cube = iris.cube.Cube([i])
            cube.add_aux_coord(iris.coords.DimCoord([t], 'time'))
            cubelist.append(cube)
        self.cubelist = iris.cube.CubeList(cubelist)
        
    def test_order(self):
        ordered_cubes =  order_cubelist(self.cubelist)
        self.assertEqual(ordered_cubes[0].coord('time').points[0], 1)
        self.assertEqual(ordered_cubes[4].coord('time').points[0], 5)


class Test_AreaBounds(unittest.TestCase):
    def setUp(self):
        self.AreaBounds = AreaBounds(['x_min', 'x_max', 'y_min', 'y_max'], [-180,180,-90,90])

    def test_bounds_range(self):
        self.assertEqual(self.AreaBounds.bounds_range, [-180,180,-90,90])

    def test_bad_input(self):
        self.assertRaises(ValueError, AreaBounds, 
                          ['x_mn', 'x_mx', 'y_mn', 'y_mx'])
        self.assertRaises(AssertionError, AreaBounds, 
                          ['x_min', 'x_max'])


class Test_check_area_bounds(Test_AreaBounds):
    def setUp(self):
        super(Test_check_area_bounds, self).setUp()
        
    def test_good_bounds(self):
        bounds = self.AreaBounds.check_area_bounds([-100,-70,40,60])
        self.assertEqual(bounds, [-100,-70,40,60])
        
    def test_bad_bounds(self):
        # Test out of bounds
        self.assertRaises(ValueError, self.AreaBounds.check_area_bounds,
                          [-181,-70,40,60])
        # Test lon_min larger than lon_max
        self.assertRaises(ValueError, self.AreaBounds.check_area_bounds,
                          [-60,-70,40,60])
        

class Test_get_area_bound_constraints(Test_AreaBounds):
    def setUp(self):
        super(Test_get_area_bound_constraints, self).setUp()
        self.bounds = [-100,-70,40,60]
        
    def test_constrates(self):
        cons = self.AreaBounds.get_area_bound_constraints(self.bounds)
        self.assertEqual(type(cons), iris._constraints.ConstraintCombination)


class Test_get_cube_area_bounds(Test_AreaBounds):
    def setUp(self):
        super(Test_get_cube_area_bounds, self).setUp()
        self.cube = iris.load(iris.sample_data_path('uk_hires.pp'))[0]
        
    def test_area_bounds(self):
        actual_bounds = [357.4939880371094, 
                         360.0049743652344, 
                         0.14430022239685059, 
                         2.8847999572753906]
        given_bounds = self.AreaBounds.get_cube_area_bounds(self.cube)
        self.assertEqual(given_bounds, actual_bounds)

class Test_check_cube_area_bounds(Test_AreaBounds):
    def setUp(self):
        super(Test_check_cube_area_bounds, self).setUp()
        self.cube = iris.load_cube(iris.sample_data_path('pre-industrial.pp'))
        # self.AreaBounds.bounds_range are [-180,180,-90,90]
        # self.cube bounds are [0, 360.,-90., 90.]
        # This function should recognise the difference and shift the 
        # coordinate system without effecting the data.
    
    def test_shift(self):
        # Get data from specific coordinate.
        old_lon_data = self.cube.extract(iris.Constraint(longitude=345)).data
        shifted_cube = self.AreaBounds.check_cube_area_bounds(self.cube)
        # Get data from the same place using shifted coordinate.
        new_lon_data = shifted_cube.extract(iris.Constraint(longitude=(345-360))).data
        self.assertEqual(old_lon_data[10], new_lon_data[10])
        
    
    def test_no_shift(self):
        # Given bounds which still fall within cube bound range, check there is 
        # no shift.
        bounds = [20,30,40,50]
        new_lon_point = self.cube.coord('longitude').points[0]
        new_cube = self.AreaBounds.check_cube_area_bounds(self.cube, bounds=bounds)
        self.assertEqual(new_lon_point, new_cube.coord('longitude').points[0])
        

class Test_coordinate_shift(unittest.TestCase):
    def setUp(self):
        self.cube = iris.load_cube(iris.sample_data_path('pre-industrial.pp'))
        self.coordinate = 'longitude'
        self.new_range = [-180, 180]
    
    def test_shift(self):
        old_lon_data = self.cube.extract(iris.Constraint(longitude=345)).data
        shifted_cube = AreaBounds.coordinate_shift(self.cube, self.coordinate, self.new_range)
        new_lon_data = shifted_cube.extract(iris.Constraint(longitude=(345-360))).data
        self.assertEqual(old_lon_data[20], new_lon_data[20])
        
    def test_new_range_error(self):
        # Check error when new range lies within cube range.
        new_range = [10, 200]
        self.assertRaises(UserWarning, AreaBounds.coordinate_shift, self.cube, self.coordinate, new_range)
        
        