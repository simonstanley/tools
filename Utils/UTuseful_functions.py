import unittest
from useful_functions import *
import inspect

this_file = inspect.getfile(inspect.currentframe())

class Test_iso_time_converter(unittest.TestCase):
    def setUp(self):
        self.datetime = datetime.datetime(2000,2,3,4)
        self.iso_time = '2000-02-03T04:00:00'
        
    def test_convert_to_iso(self):
        iso_time = iso_time_converter(self.datetime)
        self.assertEqual(iso_time, self.iso_time)
        
    def test_convert_from_iso(self):
        datetime = iso_time_converter(self.iso_time)
        self.assertEqual(datetime, self.datetime)
               

class Test_get_all_dates_between_time_bounds(unittest.TestCase):
    def setUp(self):
        self.time_start = datetime.datetime(2000,2,3)
        self.time_end   = datetime.datetime(2000,2,6)
        self.format     = '%d/%m/%Y'
        
    def test_days_between_bounds(self):
        dates = get_all_dates_between_time_bounds(self.time_start, self.time_end)
        for i, date in enumerate(dates):
            day = i+3
            self.assertEqual(datetime.datetime(2000,2,day), date)
    
    def test_formatting(self):
        dates = get_all_dates_between_time_bounds(self.time_start, self.time_end,
                                                     self.format)
        self.assertEqual(dates[1], '04/02/2000')
        
    def test_dates_in_wrong_order(self):
        self.assertRaises(AssertionError, get_all_dates_between_time_bounds,
                          self.time_end, self.time_start)
                

class Test_is_leap_year(unittest.TestCase):
    def test_leaps(self):
        self.assertTrue(is_leap_year(2000))
        self.assertTrue(is_leap_year(1996))
    
    def test_non_leaps(self):
        self.assertFalse(is_leap_year(1900))
        self.assertFalse(is_leap_year(1987))


class Test_round_data_limits(unittest.TestCase):
    def test_scalar_rounding(self):
        dmin, dmax = round_data_limits(2.34, 5.66, 0.1)
        self.assertAlmostEqual(dmin, 2.3)
        self.assertAlmostEqual(dmax, 5.7)


class Test_find_integer_divisor(unittest.TestCase):
    def test_result(self):
        result = find_integer_divisor(21, 5, 1)
        self.assertEqual(result, 7)


class Test_DirectorySearcher(unittest.TestCase):
    def setUp(self):
        self.path_list = this_file.split('/')
        self.directory = ('/').join(self.path_list[:-1])
        self.searcher = DirectorySearcher(self.directory, 
                                          look_in='', 
                                          search_for=self.path_list[:-1], 
                                          file_type='both',
                                          current_dir=False, 
                                          full_name=True)
        

class Test__get_directory(Test_DirectorySearcher):
    def setUp(self):
        super(Test__get_directory, self).setUp()
        self.searcher.search_in = self.searcher.look_in[0]
        
    def test_directory_page(self):
        dir_contents = os.listdir(self.directory)
        dir_page = self.searcher._get_directory()
        self.assertEqual(dir_page, dir_contents)
        

class Test__directory_search_in(Test_DirectorySearcher):
    def setUp(self):
        super(Test__directory_search_in, self).setUp()
        
    def test_files(self):
        self.searcher.directory = False
        self.searcher.files = True
        dir_results = self.searcher._directory_search_in()
        for dir_file in dir_results:
            self.assertNotEqual(dir_file[-1], '/')
    
    def test_directories(self):
        self.searcher.directory = True
        self.searcher.files = False
        results = self.searcher._directory_search_in()
        for dir_file in results:
            self.assertEqual(dir_file[-1], '/')
        
    
class Test__search_directory(Test_DirectorySearcher):
    def setUp(self):
        super(Test__search_directory, self).setUp()
        self.searcher.search_in = self.searcher.look_in[0]
        
    def test_result(self):
        results, directories = self.searcher._search_directory()
        for dir_file in directories:
            self.assertEqual(dir_file[-1], '/')
        self.assertIn(self.path_list[-1], results)
        

class Test__OPeNDAP_get_directory(unittest.TestCase):
    def setUp(self):
        self.searcher = DirectorySearcher('http://ukmo-01.cems.rl.ac.uk:8080/opendap/', 
                                          look_in='glosea5/', 
                                          search_for='20140205', 
                                          file_type='both',
                                          current_dir=False, 
                                          full_name=True)

    def test_directory_page(self):
        self.searcher.search_in = self.searcher.look_in[0]
        dir_page = self.searcher._OPeNDAP_get_directory()
        self.assertTrue(dir_page)
        

class Test__OPeNDAP_search_in(Test__OPeNDAP_get_directory):
    def setUp(self):
        super(Test__OPeNDAP_search_in, self).setUp()
        
    def test_files(self):
        self.searcher.directory = False
        self.searcher.files = True
        dir_results = self.searcher._OPeNDAP_search_in()
        for dir_file in dir_results:
            self.assertNotEqual(dir_file[-1], '/')
    
    def test_directories(self):
        self.searcher.directory = True
        self.searcher.files = False
        results = self.searcher._OPeNDAP_search_in()
        for dir_file in results:
            self.assertEqual(dir_file[-1], '/')


class Test__OPeNDAP_search_directory(Test__OPeNDAP_get_directory):
    def setUp(self):
        super(Test__OPeNDAP_search_directory, self).setUp()
        
    def test_result(self):
        self.searcher.search_in = self.searcher.look_in[0]
        results, directories = self.searcher._OPeNDAP_search_directory()
        for dir_file in directories:
            self.assertEqual(dir_file[-1], '/')
        self.assertEqual(results, ['glosea5/20140205/'])


class Test__directory_search_for(Test_DirectorySearcher):
    def setUp(self):
        super(Test__directory_search_for, self).setUp()
        
    def test_search_results(self):
        results = self.searcher._directory_search_for()
        self.assertIn(self.path_list[-1], results)
        
        
class Test_search(Test_DirectorySearcher):
    def setUp(self):
        super(Test_search, self).setUp()
        
    def test_search_results(self):
        # full_path = True
        results = self.searcher.search()
        self.assertIn(this_file, results)


class Test_search_directory(unittest.TestCase):
    def test_search_results(self):
        path_list = this_file.split('/')
        directory = ('/').join(path_list[:-1])
        results = search_directory(directory, 
                                   look_in='', 
                                   search_for=path_list[:-1])
        self.assertIn(path_list[-1], results)
             