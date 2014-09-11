import unittest
from GetData import *
import inspect

this_file = inspect.getfile(inspect.currentframe())
    
def create_cube(shape, realizations_from=1,
                 fcst_date_from=datetime.datetime(2014,3,18,12),
                 inti_date_from=datetime.datetime(2014,3,10,12)):
    
    data = numpy.random.rand(*shape)
    
    real_data  = range(realizations_from, shape[0]+realizations_from)
    real_coord = iris.coords.DimCoord(real_data, standard_name='realization')
    
    time_unit = iris.unit.Unit('hours since 1970-01-01 00:00:00')
    time_data = []
    date = fcst_date_from
    for i in range(shape[1]):
        this_date = date + datetime.timedelta(days=i)
        time_data.append(cube_time_converter(this_date, time_unit))
    time_coord = iris.coords.DimCoord(time_data, standard_name='time', units=time_unit)
    
    coord_sys  = iris.coord_systems.GeogCS(6371229.0)
    long_data  = numpy.linspace(-180., 180., shape[2])
    long_coord = iris.coords.DimCoord(long_data, standard_name='longitude', coord_system=coord_sys)
    
    lat_data   = numpy.linspace(-90., 90., shape[3])
    lat_coord  = iris.coords.DimCoord(lat_data, standard_name='latitude', coord_system=coord_sys)
    
    init_data = []
    date = inti_date_from
    for i in range(shape[0]):
        this_date = date + datetime.timedelta(days=i)
        init_data.append(cube_time_converter(this_date, time_unit))
    init_coord = iris.coords.AuxCoord(init_data, standard_name='forecast_reference_time', units=time_unit)
    
    
    cube = iris.cube.Cube(data, 
                          dim_coords_and_dims=[(real_coord, 0),
                                                     (time_coord, 1),
                                                     (long_coord, 2),
                                                     (lat_coord,  3)],
                          aux_coords_and_dims=[(init_coord, 0)])
    
    return cube


class Test_remove_coords(unittest.TestCase):
    def setUp(self):
        self.cube = create_cube((4,3,10,10))
        
    def test_coord_removed(self):
        self.assertTrue(self.cube.coord('forecast_reference_time'))
        cube = remove_coords(self.cube, 'forecast_reference_time')
        self.assertNotIn('forecast_reference_time', [coord.name() for coord in cube.coords()])
        

class Test_check_dates(unittest.TestCase):
    def test_bad_dates(self):
        dates = [datetime.datetime(2000,1,2), 12345.]
        self.assertRaises(TypeError, check_dates, dates)


class Test_change_zeroth_hour(unittest.TestCase):
    def test_change(self):
        dates = [datetime.datetime(2000,1,2,0), datetime.datetime(2000,1,2,1)]
        new_dates = change_zeroth_hour(dates)
        self.assertEqual(new_dates[0].hour, 12)
        self.assertEqual(new_dates[1].hour, 1)


class Test_CubeData(unittest.TestCase):
    def setUp(self):
        self.cube = create_cube((3,4,20,20))
        self.CubeData = CubeData(self.cube)


class Test_CubeListData(unittest.TestCase):
    def setUp(self):
        cube1 = create_cube((3,4,20,20))
        cube2 = create_cube((3,4,20,20))
        cube3 = create_cube((3,4,20,20))
        self.cubelist = iris.cube.CubeList([cube1, cube2, cube3])
        self.CubeListData = CubeListData(self.cubelist)


class Test_ForecastData(unittest.TestCase):
    def setUp(self):
        self.variable    = 'surface_temperature'
        self.init_dates  = [datetime.datetime(2011,7,19), datetime.datetime(2011,7,20)]
        self.fcst_dates  = [datetime.datetime(2011,9,16), datetime.datetime(2011,11,16)]
        self.area_bounds = [-45,-25,60,80]
        self.directory   = iris.sample_data_path('GloSea4')
        self.fcst_dates_bounds    = True
        self.area_bounds_format   = ['x_min','x_max','y_min','y_max']
        self.area_bounds_range    = [-180, 180, -90, 90]
        self.init_date_str_format = '%y%m%d'
        self.unwanted_coords      = ['forecast_period']
    
    def initialise(self):
        return ForecastData(self.directory,
                            self.variable, 
                            self.init_dates, 
                            self.fcst_dates, 
                            self.area_bounds,
                            self.fcst_dates_bounds,
                            self.area_bounds_format,
                            self.area_bounds_range,
                            self.init_date_str_format,
                            self.unwanted_coords)
        
    def test_initialise(self):
        self.assertTrue(self.initialise())


class Test_HindcastData(unittest.TestCase):
    def setUp(self):
        self.variable    = 'surface_temperature'
        self.init_dates  = [datetime.datetime(2011,7,19), datetime.datetime(2011,7,20)]
        self.fcst_dates  = [datetime.datetime(2011,9,16), datetime.datetime(2011,11,16)]
        self.years       = [2011]
        self.area_bounds = [-45,-25,60,80]
        self.directory   = iris.sample_data_path('GloSea4')
        self.fcst_dates_bounds    = True
        self.years_bounds         = False
        self.area_bounds_format   = ['x_min','x_max','y_min','y_max']
        self.area_bounds_range    = [-180, 180, -90, 90]
        self.init_date_str_format = '%y%m%d'
        self.unwanted_coords      = ['forecast_period']
    
    def initialise(self):
        return HindcastData(self.directory,
                            self.variable, 
                            self.init_dates, 
                            self.fcst_dates, 
                            self.years,
                            self.area_bounds,
                            self.fcst_dates_bounds,
                            self.years_bounds,
                            self.area_bounds_format,
                            self.area_bounds_range,
                            self.init_date_str_format,
                            self.unwanted_coords)
        
    def test_initialise(self):
        self.assertTrue(self.initialise())


class Test_ObservationData(unittest.TestCase):
    def setUp(self):
        self.variable     = 'surface_temperature'
        self.dates        = [datetime.datetime(2011,9,16), datetime.datetime(2011,11,16)]
        self.years        = [2011]
        self.area_bounds  = [-45,-25,60,80]
        self.directory    = iris.sample_data_path('GloSea4')
        self.dates_bounds = True
        self.years_bounds = False
        self.area_bounds_format   = ['x_min','x_max','y_min','y_max']
        self.area_bounds_range    = [-180, 180, -90, 90]
        self.unwanted_coords      = []
    
    def initialise(self):
        return ObservationData(self.directory,
                               self.variable, 
                               self.dates,
                               self.years,
                               self.area_bounds,
                               self.dates_bounds,
                               self.years_bounds,
                               self.area_bounds_format,
                               self.area_bounds_range,
                               self.unwanted_coords)
        
    def test_initialise(self):
        self.assertTrue(self.initialise())


class Test_FourierClimatology(unittest.TestCase):
    def setUp(self):
        self.variable    = 'air_temperature'
        self.init_dates  = [datetime.datetime(2011,7,19), datetime.datetime(2011,7,20)]
        self.dates       = [datetime.datetime(2011,7,21), datetime.datetime(2011,7,22)]
        self.area_bounds = [-45,-25,60,80]
        self.directory   = '/data/local/sstanley/forecast_page/fourier_coeff/means/'
        self.dates_bounds       = True
        self.area_bounds_format = ['x_min','x_max','y_min','y_max']
        self.area_bounds_range  = [-180, 180, -90, 90]
        self.lag_str_format     = '%03d'
        self.unwanted_coords    = ['forecast_period']
    
    def initialise(self):
        return FourierClimatology(self.directory,
                                  self.variable, 
                                  self.init_dates, 
                                  self.dates, 
                                  self.area_bounds,
                                  self.dates_bounds,
                                  self.area_bounds_format,
                                  self.area_bounds_range,
                                  self.lag_str_format,
                                  self.unwanted_coords)
        
    def test_initialise(self):
        self.assertTrue(self.initialise())
    

class Test__find_files(Test_CubeData):
    def test_find(self):
        path_list = this_file.split('/')
        directory = ('/').join(path_list[:-1])
        file_name = path_list[-1]
        result = self.CubeData._find_files(directory, '', file_name, False, False)
        self.assertEqual(result[0], this_file)
    
    def test_error(self):
        self.assertRaises(UserWarning, self.CubeData._find_files, '/nonsense/1/2/3/', '', 'balony.cc', False, True)
        
        
class Test__load_all_cubes(Test_ForecastData):
    def setUp(self):
        super(Test__load_all_cubes, self).setUp()
        self.inst = self.initialise()
    
    def test_load(self):
        files_to_load = self.inst._find_files(self.directory, '', ['ensemble_000.pp', 'ensemble_001.pp'], False, False)
        cubelist = self.inst._load_all_cubes(files_to_load)
        self.assertEqual(type(cubelist), iris.cube.CubeList)
        self.assertEqual(self.inst.xy_coords, ['longitude', 'latitude'])
        self.assertEqual(self.inst.time_unit, iris.unit.Unit('hours since 1970-01-01 00:00:00'))


class Test__sort_cubelist(Test_ForecastData):
    def setUp(self):
        super(Test__sort_cubelist, self).setUp()
        self.Fdata = self.initialise()
        self.cube = create_cube((4,3,40,40))
        self.Fdata.xy_coords = [coord.name() for coord in get_xy_coords(self.cube)]
    
    def test_realisations(self):
        # Two cubes with the same realisations.
        cube2 = create_cube((4,3,40,40))
        cubelist = iris.cube.CubeList([self.cube, cube2])
        cube = self.Fdata._sort_cubelist(cubelist)
        self.assertEqual(range(1,9), list(cube.coord('realization').points))
        
    def test_time(self):
        # Two cubes with differing but overlapping times.
        cube2 = create_cube((4,2,40,40))
        cubelist = iris.cube.CubeList([self.cube, cube2])
        cube = self.Fdata._sort_cubelist(cubelist)
        self.assertEqual(3, len(cube.coord('time').points))
        
    def test_area_bounds(self):
        # Test area extraction.
        cube2 = create_cube((4,3,40,40))
        cubelist = iris.cube.CubeList([self.cube, cube2])
        cube = self.Fdata._sort_cubelist(cubelist)
        self.assertEqual(len(cube.coord('longitude').points), 2)
        self.assertEqual(len(cube.coord('latitude').points), 4)


class Test__get_model_metadata(Test_CubeData):
    def setUp(self):
        super(Test__get_model_metadata, self).setUp()
        
    def test_metadata(self):
        self.cube = create_cube((3,1,20,20))
        metadata = self.CubeData._get_model_metadata(self.cube)
        self.assertEqual(metadata['LONGITUDE_BOUNDS'], [-180.0, 180.0])
        self.assertEqual(metadata['LATITUDE_BOUNDS'], [-90.0, 90.0])
        self.assertEqual(metadata['FORECAST_DATES'], [datetime.datetime(2014, 3, 18, 12, 0)])


class Test__get_obs_metadata(Test_CubeData):
    def setUp(self):
        super(Test__get_obs_metadata, self).setUp()
        
    def test_metadata(self):
        self.cube = create_cube((3,1,20,20))
        metadata = self.CubeData._get_obs_metadata(self.cube)
        self.assertEqual(metadata['LONGITUDE_BOUNDS'], [-180.0, 180.0])
        self.assertEqual(metadata['LATITUDE_BOUNDS'], [-90.0, 90.0])
        self.assertEqual(metadata['DATES'], [datetime.datetime(2014, 3, 18, 12, 0)])


class Test_cube_coordinate_analysis(Test_CubeData):
    def setUp(self):
        super(Test_cube_coordinate_analysis, self).setUp()
        self.cube = create_cube((1,3,1,1))
    
    def test_mean(self):
        cube = self.CubeData.cube_coordinate_analysis(self.cube, 'time')
        mean_val = numpy.mean(self.cube.data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        

class Test_cube_ensemble_mean(Test_CubeData):
    def setUp(self):
        super(Test_cube_ensemble_mean, self).setUp()
        self.cube = create_cube((3,1,1,1))
    
    def test_mean(self):
        cube = self.CubeData.cube_ensemble_mean(self.cube)
        mean_val = numpy.mean(self.cube.data)
        self.assertEqual(cube.data[0][0][0], mean_val)


class Test_cube_area_analysis(Test_CubeData):
    def setUp(self):
        super(Test_cube_area_analysis, self).setUp()
        self.cube = create_cube((1,1,20,20))
    
    def test_mean(self):
        cube = self.CubeData.cube_area_analysis(self.cube)
        mean_val = numpy.mean(self.cube.data)
        # Area weighting is applied so numpy.mean is an approximation.
        self.assertAlmostEqual(cube.data[0][0], mean_val, places=1)


class Test_extract_area_bounds(Test_ForecastData):
    def setUp(self):
        super(Test_extract_area_bounds, self).setUp()
        self.Fdata = self.initialise()
        self.cube = create_cube((4,3,73,37))
        
    def test_extraction(self):
        cube = self.Fdata.extract_area_bounds(cubes=self.cube)
        self.assertEqual(min(cube.coord('longitude').points), self.Fdata.area_bounds[0])
        self.assertEqual(max(cube.coord('longitude').points), self.Fdata.area_bounds[1])
        self.assertEqual(min(cube.coord('latitude').points), self.Fdata.area_bounds[2])
        self.assertEqual(max(cube.coord('latitude').points), self.Fdata.area_bounds[3])
    
    def test_bad_bounds(self):
        self.assertRaises(ValueError, self.Fdata.extract_area_bounds, [30, 20, 20, 30], self.cube)
    
    def test_no_loaded_data(self):
        self.assertRaises(AssertionError, self.Fdata.extract_area_bounds)
    

class Test_extract_dates(Test_ForecastData):
    def setUp(self):
        super(Test_extract_dates, self).setUp()
        self.Fdata = self.initialise()
        self.cube = create_cube((4,5,20,20))
    
    def test_extraction(self):
        dates = [datetime.datetime(2014,3,19,12), datetime.datetime(2014,3,20,12)]
        cube = self.Fdata.extract_dates(dates, self.cube)
        time_unit = cube.coord('time').units
        self.assertEqual(min(cube.coord('time').points), cube_time_converter(dates[0], time_unit))
        self.assertEqual(max(cube.coord('time').points), cube_time_converter(dates[1], time_unit))

    def test_bad_dates(self):
        self.assertRaises(UserWarning, self.Fdata.extract_dates, [datetime.datetime(2013,1,1,12)], self.cube)
        
    def test_no_loaded_data(self):
        self.assertRaises(AssertionError, self.Fdata.extract_dates)


class Test_extract_coordinate_levels(Test_ForecastData):
    def setUp(self):
        super(Test_extract_coordinate_levels, self).setUp()
        self.Fdata = self.initialise()
        self.cube = create_cube((4,3,20,20))
        
    def test_extraction(self):
        cube = self.Fdata.extract_coordinate_levels('realization', [2,3], self.cube)
        self.assertEqual(min(cube.coord('realization').points), 2)
        self.assertEqual(max(cube.coord('realization').points), 3)
        
    def test_no_coordinate(self):
        self.assertRaises(UserWarning, self.Fdata.extract_coordinate_levels, 'bad_coord', [2,3], self.cube)

    def test_no_loaded_data(self):
        self.assertRaises(AssertionError, self.Fdata.extract_coordinate_levels, 'bad_coord', [2,3])
    

class Test_configure_coord_names(Test_CubeData):
    def setUp(self):
        super(Test_configure_coord_names, self).setUp()
        
    def test_coord_name_change(self):
        self.assertEqual(self.CubeData.time_coord, 'time')
        self.CubeData.configure_coord_names(period_name='t')
        self.assertEqual(self.CubeData.time_coord, 't')


class Test_MetaData(Test_CubeData):
    def setUp(self):
        super(Test_MetaData, self).setUp()
        
    def test_display(self):
        ex_dict = {'DATES':'test'}
        metadata = self.CubeData.MetaData(ex_dict)
        self.assertEqual(metadata.__str__(), '--Metadata--\n\n---\nDATES: test\n---\n')


class Test_time_analysis(Test_CubeData):
    def setUp(self):
        super(Test_time_analysis, self).setUp()
        self.CubeData.cube = create_cube((1,3,1,1))
    
    def test_mean(self):
        cube = self.CubeData.time_analysis()
        mean_val = numpy.mean(self.CubeData.cube.data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('time_analysis', self.CubeData.processes)
        

class Test_ensemble_mean(Test_CubeData):
    def setUp(self):
        super(Test_ensemble_mean, self).setUp()
        self.CubeData.cube = create_cube((3,1,1,1))
    
    def test_mean(self):
        cube = self.CubeData.ensemble_mean()
        mean_val = numpy.mean(self.CubeData.cube.data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('ensemble_mean', self.CubeData.processes)


class Test_area_analysis(Test_CubeData):
    def setUp(self):
        super(Test_area_analysis, self).setUp()
        self.CubeData.cube = create_cube((1,1,20,20))
    
    def test_mean(self):
        cube = self.CubeData.area_analysis()
        mean_val = numpy.mean(self.CubeData.cube.data)
        # Area weighting is applied so numpy.mean is an approximation.
        self.assertAlmostEqual(cube.data[0][0], mean_val, places=1)
        self.assertIn('area_analysis', self.CubeData.processes)


class Test_coordinate_analysis(Test_CubeData):
    def setUp(self):
        super(Test_coordinate_analysis, self).setUp()
        self.CubeData.cube = create_cube((1,3,1,1))
    
    def test_mean(self):
        cube = self.CubeData.coordinate_analysis('time')
        mean_val = numpy.mean(self.CubeData.cube.data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('time_analysis', self.CubeData.processes)


class Test_time_analysis_cubelist(Test_CubeListData):
    def setUp(self):
        super(Test_time_analysis_cubelist, self).setUp()
        self.CubeListData.cubelist = iris.cube.CubeList([create_cube((1,3,1,1))])
    
    def test_mean(self):
        cube = self.CubeListData.time_analysis()[0]
        mean_val = numpy.mean(self.CubeListData.cubelist[0].data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('time_analysis', self.CubeListData.processes)
        

class Test_ensemble_mean_cubelist(Test_CubeListData):
    def setUp(self):
        super(Test_ensemble_mean_cubelist, self).setUp()
        self.CubeListData.cubelist = iris.cube.CubeList([create_cube((3,1,1,1))])
    
    def test_mean(self):
        cube = self.CubeListData.ensemble_mean()[0]
        mean_val = numpy.mean(self.CubeListData.cubelist[0].data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('ensemble_mean', self.CubeListData.processes)


class Test_area_analysis_cubelist(Test_CubeListData):
    def setUp(self):
        super(Test_area_analysis_cubelist, self).setUp()
        self.CubeListData.cubelist = iris.cube.CubeList([create_cube((1,1,20,20))])
    
    def test_mean(self):
        cube = self.CubeListData.area_analysis()[0]
        mean_val = numpy.mean(self.CubeListData.cubelist[0].data)
        # Area weighting is applied so numpy.mean is an approximation.
        self.assertAlmostEqual(cube.data[0][0], mean_val, places=1)
        self.assertIn('area_analysis', self.CubeListData.processes)


class Test_coordinate_analysis_cubelist(Test_CubeListData):
    def setUp(self):
        super(Test_coordinate_analysis_cubelist, self).setUp()
        self.CubeListData.cubelist = iris.cube.CubeList([create_cube((1,3,1,1))])
    
    def test_mean(self):
        cube = self.CubeListData.coordinate_analysis('time')[0]
        mean_val = numpy.mean(self.CubeListData.cubelist[0].data)
        self.assertEqual(cube.data[0][0][0], mean_val)
        self.assertIn('time_analysis', self.CubeListData.processes)


class Test_period_check(unittest.TestCase):
    def setUp(self):
        self.init_dates = [datetime.datetime(2000,1,2), datetime.datetime(2000,1,3)]
    
    def test_overlap_dates(self):
        fcst_dates = [datetime.datetime(2000,1,2), datetime.datetime(2000,1,3),datetime.datetime(2000,1,4)]
        self.assertRaises(ValueError, ForecastData.period_check, self.init_dates, fcst_dates)
        
    def test_hours_not_zero(self):
        fcst_dates = [datetime.datetime(2000,1,4)]
        init_dates, fcst_dates = ForecastData.period_check(self.init_dates, fcst_dates)
        self.assertNotEqual(fcst_dates[0].hour, 0)
        

class Test_load_fcst(Test_ForecastData):
    def setUp(self):
        super(Test_load_fcst, self).setUp()
        self.Fdata = self.initialise()
        
    def test_load(self):
        cube = self.Fdata.load(files_to_search_for='ensemble')
        self.assertTrue(cube)


class Test_get_metadata_fcst(Test_ForecastData):
    def setUp(self):
        super(Test_get_metadata_fcst, self).setUp()
        self.Fdata = self.initialise()
        self.Fdata.cube = create_cube((3,3,20,20))
        self.Fdata.xy_coords = [coord.name() for coord in get_xy_coords(self.Fdata.cube)]
        
    def test_metadata(self):
        metadata = self.Fdata._get_metadata()
        self.assertEqual(type(metadata), self.Fdata.MetaData)
        self.assertEqual(metadata.meta_dict['DATA_TYPE'], 'Forecast Data')


class Test__sort_hcst_dates(Test_HindcastData):
    def setUp(self):
        super(Test__sort_hcst_dates, self).setUp()
        self.Hdata = self.initialise()
        
    def test_years_in_dates_dont_matter(self):
        init_dates = [datetime.datetime(1990, 1, 1), datetime.datetime(1990, 1, 2)]
        fcst_dates = [datetime.datetime(1990, 2, 1), datetime.datetime(1990, 2, 2)]
        years = [2000, 2001]
        init_dates, fcst_dates = self.Hdata._sort_hcst_dates(init_dates, fcst_dates, years)
        self.assertEqual(init_dates, [datetime.datetime(2000, 1, 1, 0, 0), datetime.datetime(2000, 1, 2, 0, 0), datetime.datetime(2001, 1, 1, 0, 0), datetime.datetime(2001, 1, 2, 0, 0)])
        self.assertEqual(fcst_dates, [[datetime.datetime(2000, 2, 1, 12, 0), datetime.datetime(2000, 2, 2, 12, 0)], [datetime.datetime(2001, 2, 1, 12, 0), datetime.datetime(2001, 2, 2, 12, 0)]])

    def test_overlapping_years(self):
        # The actual year follows the first forecast date.
        init_dates = [datetime.datetime(1990, 11, 1)]
        fcst_dates = [datetime.datetime(1990, 12, 30), datetime.datetime(1991, 1, 1)]
        years = [2000, 2001]
        init_dates, fcst_dates = self.Hdata._sort_hcst_dates(init_dates, fcst_dates, years)
        self.assertEqual(fcst_dates, [[datetime.datetime(2000, 12, 30, 12, 0), datetime.datetime(2001, 1, 1, 12, 0)], [datetime.datetime(2001, 12, 30, 12, 0), datetime.datetime(2002, 1, 1, 12, 0)]])

        
class Test__sort_data_hcst(Test_HindcastData):
    def setUp(self):
        super(Test__sort_data_hcst, self).setUp()
        self.init_dates = [datetime.datetime(2013,3,10,12), datetime.datetime(2013,3,12,12)]
        self.fcst_dates = [datetime.datetime(2013,3,18,12), datetime.datetime(2013,3,20,12)]
        self.years = [2012,2013,2014]
        self.Hdata = self.initialise()
        cube1 = create_cube((3,2,40,40))
        cube2 = create_cube((3,2,40,40), fcst_date_from=datetime.datetime(2013,3,18,12),
                            inti_date_from=datetime.datetime(2013,3,10,12))
        cube3 = create_cube((3,2,40,40), fcst_date_from=datetime.datetime(2012,3,18,12),
                            inti_date_from=datetime.datetime(2012,3,10,12))
        self.cubelist = iris.cube.CubeList([cube1, cube2, cube3])
        self.time_unit = iris.unit.Unit('hours since 1970-01-01 00:00:00')
        
    def test_cubelist_split_into_years(self):
        cubelist = self.Hdata._sort_data(self.cubelist)
        years = [cube_time_converter(cube.coord('time').points[0], self.time_unit).year for cube in cubelist]
        self.assertEqual(years, self.years)


class Test_load_hcst(Test_HindcastData):
    def setUp(self):
        super(Test_load_hcst, self).setUp()
        self.Hdata = self.initialise()
        
    def test_load(self):
        cubelist = self.Hdata.load(files_to_search_for='ensemble')
        self.assertTrue(cubelist)


class Test_get_metadata_hcst(Test__sort_data_hcst):
    def setUp(self):
        super(Test_get_metadata_hcst, self).setUp()
        self.Hdata = self.initialise()
        self.Hdata.cubelist = self.cubelist
        self.Hdata.xy_coords = [coord.name() for coord in get_xy_coords(self.Hdata.cubelist[0])]
        
    def test_metadata(self):
        metadata = self.Hdata._get_metadata()
        self.assertEqual(type(metadata), self.Hdata.MetaData)
        self.assertEqual(metadata.meta_dict['DATA_TYPE'], 'Hindcast Data')


class Test_create_climatology_hcst(Test_HindcastData):
    def setUp(self):
        super(Test_create_climatology_hcst, self).setUp()
        self.Hdata = self.initialise()
        self.Hdata.load('ensemble')
        
    def test_clim_cube(self):
        clim_cube = self.Hdata.create_climatology()
        self.assertEqual(len(clim_cube.coord('realization').points), 1)
        self.assertEqual(len(clim_cube.coord('time').points), 3)


class Test__sort_dates(Test_ObservationData):
    def setUp(self):
        super(Test__sort_dates, self).setUp()
        self.Obsdata = self.initialise()
        
    def test_years_in_dates_dont_matter(self):
        dates = [datetime.datetime(1990, 2, 1), datetime.datetime(1990, 2, 2)]
        years = [2000, 2001]
        dates = self.Obsdata._sort_dates(dates, years)
        self.assertEqual(dates, [[datetime.datetime(2000, 2, 1, 12, 0), datetime.datetime(2000, 2, 2, 12, 0)], [datetime.datetime(2001, 2, 1, 12, 0), datetime.datetime(2001, 2, 2, 12, 0)]])

    def test_overlapping_years(self):
        dates = [datetime.datetime(1990, 12, 30), datetime.datetime(1991, 1, 1)]
        years = [2000, 2001]
        dates = self.Obsdata._sort_dates(dates, years)
        self.assertEqual(dates, [[datetime.datetime(2000, 12, 30, 12, 0), datetime.datetime(2001, 1, 1, 12, 0)], [datetime.datetime(2001, 12, 30, 12, 0), datetime.datetime(2002, 1, 1, 12, 0)]])

        
class Test__sort_data_obs(Test_ObservationData):
    def setUp(self):
        super(Test__sort_data_obs, self).setUp()
        self.dates = [datetime.datetime(2013,3,18,12), datetime.datetime(2013,3,20,12)]
        self.years = [2012,2013,2014]
        self.Obsdata = self.initialise()
        cube1 = create_cube((1,2,40,40))
        cube2 = create_cube((1,2,40,40), fcst_date_from=datetime.datetime(2013,3,18,12))
        cube3 = create_cube((1,2,40,40), fcst_date_from=datetime.datetime(2012,3,18,12))
        self.cubelist = iris.cube.CubeList([cube1, cube2, cube3])
        for cube in self.cubelist:
            cube.remove_coord('forecast_reference_time')
            cube.remove_coord('realization')
        self.time_unit = iris.unit.Unit('hours since 1970-01-01 00:00:00')
        
    def test_cubelist_split_into_years(self):
        cubelist = self.Obsdata._sort_data(self.cubelist)
        years = [cube_time_converter(cube.coord('time').points[0], self.time_unit).year for cube in cubelist]
        self.assertEqual(years, self.years)


class Test_load_obs(Test_ObservationData):
    def setUp(self):
        super(Test_load_obs, self).setUp()
        self.Obsdata = self.initialise()
        
    def test_load(self):
        cubelist = self.Obsdata.load(files_to_search_for='ensemble_000')
        self.assertTrue(cubelist)


class Test_get_metadata_obs(Test__sort_data_obs):
    def setUp(self):
        super(Test_get_metadata_obs, self).setUp()
        self.Obsdata = self.initialise()
        self.Obsdata.cubelist = self.cubelist
        self.Obsdata.xy_coords = [coord.name() for coord in get_xy_coords(self.Obsdata.cubelist[0])]
        
    def test_metadata(self):
        metadata = self.Obsdata._get_metadata()
        self.assertEqual(type(metadata), self.Obsdata.MetaData)
        self.assertEqual(metadata.meta_dict['DATA_TYPE'], 'Observation Data')


class Test_create_climatology_obs(Test_ObservationData):
    def setUp(self):
        super(Test_create_climatology_obs, self).setUp()
        self.Obsdata = self.initialise()
        self.Obsdata.load('ensemble_000')
        
    def test_clim_cube(self):
        clim_cube = self.Obsdata.create_climatology()
        self.assertEqual(len(clim_cube.coord('time').points), 3)


class Test_calculate_day_of_year(unittest.TestCase):
    def test_day_number(self):
        day_1 = FourierClimatology.calculate_day_of_year(datetime.datetime(2014, 1, 1))
        day_300 = FourierClimatology.calculate_day_of_year(datetime.datetime(2014, 10, 27))
        self.assertEqual(day_1, 1)
        self.assertEqual(day_300, 300)
        
    def test_leap(self):
        # Leap day 29th Feb should resort to the same as 28th
        feb_29 = FourierClimatology.calculate_day_of_year(datetime.datetime(2000, 2, 29))
        feb_28 = FourierClimatology.calculate_day_of_year(datetime.datetime(2000, 2, 28))
        self.assertEqual(feb_29, feb_28)


class Test__calculate_date(unittest.TestCase):
    def test_date_str(self):
        jan_1 = FourierClimatology._calculate_date(1)
        oct_27 = FourierClimatology._calculate_date(300)
        self.assertEqual(jan_1, '01-Jan')
        self.assertEqual(oct_27, '27-Oct')
        

class Test__get_lags_dict(Test_FourierClimatology):
    def setUp(self):
        super(Test__get_lags_dict, self).setUp()
        self.FClim = self.initialise()
    
    def test_dict(self):
        lags_dict = self.FClim._get_lags_dict()
        self.assertEqual(lags_dict, {1: [202], 2: [202, 203], 3: [203]})
        

class Test__calculate_clim_data(unittest.TestCase):
    def test_clim_data(self):
        fourier_coeffs = []
        for data in range(9):
            fourier_coeffs.append(iris.cube.Cube(data))
        day_of_year = 20
        self.assertAlmostEqual(FourierClimatology._calculate_clim_data(day_of_year, fourier_coeffs),
                               -29.13891, 4)


class Test__get_metadata(Test_FourierClimatology):
    def setUp(self):
        super(Test__get_metadata, self).setUp()
        self.FClim = self.initialise()
        self.FClim.cube = create_cube((1,3,20,20))
        self.FClim.cube.remove_coord('realization')
        self.FClim.cube_init_dates = ['test']
        self.FClim.cube_dates = ['test']
        self.FClim.xy_coords = [coord.name() for coord in get_xy_coords(self.FClim.cube)]
        
    def test_metadata(self):
        metadata = self.FClim._get_metadata()
        self.assertEqual(type(metadata), self.FClim.MetaData)
        self.assertEqual(metadata.meta_dict['DATA_TYPE'], 'Fourier Climatology')


class Test_advanced_search_pattern(Test_FourierClimatology):
    def setUp(self):
        super(Test_advanced_search_pattern, self).setUp()
        self.FClim = self.initialise()
    
    def test_string_updates(self):
        self.assertEqual(self.FClim.prepend_lag_string, '')
        self.FClim.advanced_search_pattern(prepend_lag_string='test')
        self.assertEqual(self.FClim.prepend_lag_string, 'test')
    
    
class Test_load_fourier(Test_FourierClimatology):
    def setUp(self):
        super(Test_load_fourier, self).setUp()
        self.FClim = self.initialise()
        
    def test_load(self):
        cube = self.FClim.load()
        self.assertTrue(cube)
        