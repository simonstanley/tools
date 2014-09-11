import datetime
import dateutil.parser
import os
import re
import math
from xml.etree.ElementTree import ElementTree
import urllib2
from urlparse import urlparse

def iso_time_converter(time):
    """
    Convert time between datetime object and iso string format.

    Args:
    
    * time: datetime, string
        Which ever type is given, it is converted to the other. If a string is
        given, it must be in an iso format.

    Returns:
        datetime or iso string

    """
    if type(time) == datetime.datetime:
        converted_time = time.isoformat()
    elif type(time) in [str, unicode]:
        converted_time = dateutil.parser.parse(time)
    else:
        raise UserWarning('Invalid time type.')
    return converted_time

def get_all_dates_between_time_bounds(time_start, time_end, str_format=None):
    """
    Return a list of all dates (i.e. days) between two points in time.
    
    Args:
    
    * time_start: datetime
    
    * time_end: datetime
    
    Kwargs:
    
    * str_format: string
        Converts the datetime objects into string formats using datetime 
        formatting codes. Default is to leave as datetimes.

    Returns:
        List of dates
    """
    assert time_start <= time_end, 'Start time must come before end time'
    all_dates = []
    current_date = time_start
    while current_date <= time_end:
        if str_format:
            formatted_current_date = current_date.strftime(str_format)
        else:
            formatted_current_date = current_date
        all_dates.append(formatted_current_date)
        current_date += datetime.timedelta(days=1)
    
    return all_dates

def is_leap_year(year):
    """
    Determine whether year is leap.
    
    Args:
    
    * year: integer
    
    Returns:
        Boolean
    
    """
    if year % 400 == 0:
        return True
    if year % 100 == 0:
        return False
    if year % 4 == 0:
        return True
    else:
        return False

def round_data_limits(dmin, dmax, scalar=1.):
    """
    Round the data limits according to the scalar. E.g. round 100.11 to 100 
    on scalar 1, and to 100.1 on scalar 0.1.  
    
    Args:
    
    * dmin: float
    
    * dmax: float
    
    Kwargs:
    
    * scalar: float
    
    Returns:
        tuple
    
    """   
    adj_dmin = math.floor(dmin / scalar) * scalar
    adj_dmax = math.ceil(dmax / scalar) * scalar
    return adj_dmin, adj_dmax

def find_integer_divisor(dividend, requested_divisor, multiple_of=0.5, 
                           search_range=4):
    """
    Find an integer value close to the requested divisor which divides the 
    dividend leaving a requested multiple. While preference is taken on 
    searching for a divisor above the requested divisor, distance from 
    requested divisor is still first priority.
    
    Args:
    
    * dividend: float
        Value to be divided.
        
    * requested_divisor: positive integer
         Value to search around looking for an even divide.
    
    Kwargs:
    
    * multiple_of: float
        Value the result of the division should be a multiple of.
        
    * search_range: integer
        How far +- from the requested divisor to search for a value which
        fulfils the request.
    
    Returns:
        The result or None if none found.
        
    
    """
    if requested_divisor < 1:
        raise UserWarning('Requested divisor must be positive.')
    if (float(dividend) / float(requested_divisor)) % multiple_of == 0:
        return requested_divisor
    
    divisor_found = False
    num_adj = 1
    while num_adj <= search_range:
        for sign in [1., -1.]:
            test_divisor = requested_divisor + (num_adj * sign)
            # Do not divide by 0.
            if test_divisor == 0:
                continue
            if (float(dividend) / test_divisor) % multiple_of == 0:
                divisor_found = True
                result = test_divisor
                break
        if divisor_found:
            break
        num_adj += 1
    if not divisor_found:
        result = None
    
    return result
    
    
class DirectorySearcher(object):
    """
    Search through a given directory. Also excepts valid OPeNDAP URL's.

    Args:
    
    * dir_path: string
        Specify the path or OPeNDAP URL to the top of the directory to be 
        searched.

    Kwargs:

    * look_in: string, list, tuple
        Specify the directory(s) to search in. Default is the top directory.

    * search_for: string, list, tuple
        Specify the string pattern(s) to search for. The search will be 
        performed in the directory(s) given by look_in and all directories 
        beneath, unless current_dir is set to True. Default is None, this
        returns all files in the look_in directory(s) only.
        
    * file_type: 'dirs', 'files'
        If set to 'dirs', only directory results are returned, if 'files', only
        files are returned. Default is both.
        
    * current_dir: boolean
        Only return matches which occur inside the directory(s) given by 
        look_in. Default is to include all directories beneath, as well.
        
    * full_name: boolean
        If set to True, the full URL's are returned (for use in loading). 
        Default is False which returns file paths relative to the top
        directory (for use in continued searching).
        
    """
    def __init__(self, dir_path, look_in='', search_for=None, file_type='both',
                  current_dir=False, full_name=False, hidden_files=False):
        
        self.dir_path = self._append_slash(dir_path)
        urlparse_result = urlparse(dir_path)
        if not urlparse_result[0]:
            self.path_type = 'local'
        else:
            self.path_type = 'opendap'
        
        if type(look_in) == str:
            look_in = [look_in]
        updated_look_in = []
        for look_in_dir in look_in:
            if look_in_dir:
                look_in_dir = self._remove_start_slash(look_in_dir)
                look_in_dir = self._append_slash(look_in_dir)
            updated_look_in.append(look_in_dir)
        self.look_in = updated_look_in
        
        self.search_for = search_for
        if type(self.search_for) == str:
            self.search_for = [self.search_for]
            
        self.current_dir = current_dir
        
        # Specify what type of file to look for.
        self.directory   = True
        self.files       = True
        if file_type == 'dirs':
            self.files = False
        elif file_type == 'files':
            self.directory = False
        self.full_name = full_name
        self.hidden_files = hidden_files
        
        self.items_found = []

    @staticmethod
    def _append_slash(path):
        """
        Add '/' to the end of path if not there.
        
        """
        if path[-1] != '/':
            path += '/'
        return path
    
    @staticmethod
    def _remove_start_slash(path):
        """
        Remove '/' from the start of path if there.
        
        """
        if path[0] == '/':
            path = path[1:]
        return path

    def _get_directory(self):
        path = '{d}{p}'.format(d=self.dir_path, p=self.search_in)
        try:
            directory_page = os.listdir(path)
        except OSError:
            raise UserWarning('Invalid directory path: %s' % path)
        return directory_page
    
    def _directory_search_in(self):
        """
        Return all files within 'look_in' directories. Unless only files or
        directories are requested.
        
        """
        results = []
        for self.search_in in self.look_in:
            directory = self._get_directory()
            for dir_file in directory:
                if self.directory and os.path.isdir(self.dir_path + self.search_in + dir_file):
                    results.append(self.search_in + dir_file + '/')
                if self.files and not os.path.isdir(self.dir_path + self.search_in + dir_file):
                    if self.hidden_files: 
                        results.append(self.search_in + dir_file)
                    else:
                        if dir_file[0] != '.':
                            results.append(self.search_in + dir_file)
        return results
    
    def _search_directory(self):
        """
        Search for 'search_for' items in directory, return all directories 
        found so search can continue in each.
        
        """
        results = []
        directories = []
        
        directory = self._get_directory()
        
        for dir_file in directory:
            for item in self.search_for:
                match = re.findall(item, dir_file)
                if match:
                    if item not in self.items_found:
                        self.items_found.append(item)
                    if self.directory and os.path.isdir(self.dir_path + self.search_in + dir_file):
                        results.append(self.search_in + dir_file + '/')
                    if self.files and not os.path.isdir(self.dir_path + self.search_in + dir_file):
                        if self.hidden_files: 
                            results.append(self.search_in + dir_file)
                        else:
                            if dir_file[0] != '.':
                                results.append(self.search_in + dir_file)
            if os.path.isdir(self.dir_path + self.search_in + dir_file):
                directories.append(self.search_in + dir_file + '/')
                
        return results, directories
        
    def _OPeNDAP_get_directory(self):
        path = '{o}{p}'.format(o=self.dir_path, p=self.search_in)
        try:
            directory_page = urllib2.urlopen(path + 'catalog.xml')
        except urllib2.HTTPError:
            raise UserWarning('Invalid OPeNDAP URL: %s' % path)
        tree = ElementTree()    
        tree.parse(directory_page)
        root = tree.getroot()
        tag = re.findall('{.*}', root.tag)[0]
        
        return root, tag
    
    def _OPeNDAP_search_in(self):
        
        results = []
        
        for self.search_in in self.look_in:
            root, tag = self._OPeNDAP_get_directory()
            if self.directory:
                for directory in root.iter('{tag}catalogRef'.format(tag=tag)):
                    results.append(self.search_in + directory.get('name') + '/')
            
            if self.files:
                for dataset in root.iter('{tag}dataset'.format(tag=tag)):
                    # If there is a '.' in last element of the path, assume it 
                    # to be a file.
                    if '.' in dataset.get('name'):
                        results.append(self.search_in + dataset.get('name'))
                    
        return results
    
    def _OPeNDAP_search_directory(self):
        
        results = []
        directories = []
        
        root, tag = self._OPeNDAP_get_directory()
        
        
        for directory in root.iter('{tag}catalogRef'.format(tag=tag)):
            if self.directory:
                for item in self.search_for:
                    match = re.findall(item, directory.get('name'))
                    if match:
                        if item not in self.items_found:
                            self.items_found.append(item)
                        results.append(self.search_in + directory.get('name') + '/')
            directories.append(self.search_in + directory.get('name') + '/')
        
        if self.files:
            for dataset in root.iter('{tag}dataset'.format(tag=tag)):
                # If there is a '.' in last element of the path, assume it to 
                # be a file.
                if '.' in dataset.get('name'):
                    for item in self.search_for:
                        match = re.findall(item, dataset.get('name'))
                        if match:
                            if item not in self.items_found:
                                self.items_found.append(item)
                            results.append(self.search_in + dataset.get('name'))
                
        return results, directories

    def _directory_search_for(self):
        """
        Search all directories beneath 'look_in' level (unless 'current_dir' is
        True) looking for 'search_for' items.
        
        """
        results = []
        directories = self.look_in
        if self.path_type == 'local':
            search_directory = self._search_directory
        elif self.path_type == 'opendap':
            search_directory = self._OPeNDAP_search_directory
        
        while True:
            next_directories = []
            for self.search_in in directories:
                search_results, directory_results = search_directory()
                for result in search_results:
                    results.append(result)
                for directory in directory_results:
                    if not self.current_dir:
                        next_directories.append(directory)      
            if not next_directories:
                break
            else:       
                directories = next_directories
        return results

    def search(self):
        """
        Perform search.
        
        Returns:
            All search results in a list.
            
        """
        if self.search_for:
            results = self._directory_search_for()
            items_not_found = set(self.search_for) - set(self.items_found)
            if items_not_found:
                print 'The following items were not found: %s' \
                       % ', '.join(items_not_found)
        else:
            if self.path_type == 'local':
                results = self._directory_search_in()
            elif self.path_type == 'opendap':
                results = self._OPeNDAP_search_in()
                
        if self.full_name:
            results = [self.dir_path+filename for filename in results]
                
        return results

def search_directory(dir_path, look_in='', search_for=None, file_type='both', 
                       current_dir=False, full_name=False, hidden_files=False):
    """
    Search through a given directory. Also excepts valid OPeNDAP URL's.

    Args:
    
    * dir_path: string
        Specify the path or OPeNDAP URL to the top of the directory to be 
        searched.

    Kwargs:

    * look_in: string, list, tuple
        Specify the directory(s) to search in. Default is the top directory.

    * search_for: string, list, tuple
        Specify the string pattern(s) to search for. The search will be 
        performed in the directory(s) given by look_in and all directories 
        beneath, unless current_dir is set to True. Default is None, this
        returns all files in the look_in directory(s) only.
        
    * file_type: 'dirs', 'files'
        If set to 'dirs', only directory results are returned, if 'files', only
        files are returned. Default is both.
        
    * current_dir: boolean
        Only return matches which occur inside the directory(s) given by 
        look_in. Default is to include all directories beneath, as well.
        
    * full_name: boolean
        If set to True, the full URL's are returned (for use in loading). 
        Default is False which returns file paths relative to the top
        directory (for use in continued searching).

    Returns:
        All search results in a list.

    """
    searcher = DirectorySearcher(dir_path, look_in, search_for, file_type, 
                                 current_dir, full_name, hidden_files)
    results = searcher.search()    
    return results
