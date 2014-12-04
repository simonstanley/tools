"""
Module for gaining quick access to specific data kept in known files and 
formats. Also for performing basic processing.

Author: S Stanley
Date: May 2014
"""
import numpy
import datetime
import inspect
from urllib2 import urlopen
from calendar import monthrange
from GetData import ObservationData

NAO_PATH      = '/project/hadobs1/OBS/pressure/NAO/actual/nao_jones_rodwell_monthly.dat'
FLOOD_PATH    = '/home/h02/sstanley/Documents/data/flood_event_count_89_13.txt'
NCICTEXT_PATH = 'http://www01/obs_dev/od4/series/text/'
EObs_PATH     = '/data/cr1/sstanley/eobs/'
NCIC_PATH     = '/project/ncic/NCIC/projects/NCIC_Gridding_Software/Test_data/'
ERA_I_PATH    = '/project/seasonal/frgo/'
ISSFCST_PATH  = '/home/h02/frgo/TEST/jhirst_plots/new_caboff_plots/plots_N216/'
ISS_SAVEFILE  = '/net/home/h02/sstanley/Documents/data/forecast_data/'

MONS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

ncic_txt_var_dict = {'Temp'     : {'name'  : 'Tmean', 
                                   'units' : 'celsius',
                                   'type'  : 'mean'},
                     'Tmax'     : {'name'  : 'Tmax', 
                                   'units' : 'celsius',
                                   'type'  : 'mean'},
                     'Tmin'     : {'name'  : 'Tmin', 
                                   'units' : 'celsius',
                                   'type'  : 'mean'},
                     'Precip'   : {'name'  : 'Rainfall', 
                                   'units' : 'mm',
                                   'type'  : 'total'},
                     'Sunshine' : {'name'  : 'Sunshine', 
                                   'units' : 'hours',
                                   'type'  : 'total'},
                     'Airfrost' : {'name'  : 'AirFrost', 
                                   'units' : 'days',
                                   'type'  : 'total'}}

ncic_var_dict = {'Precip' : {'name'     : 'Precipitation',
                             'filename' : 'rr'},
                 'Temp'   : {'name'     : 'air_temperature',
                             'filename' : 'tmean'}}

eobs_var_dict = {'precip' : {'name'     : 'thickness_of_rainfall_amount',
                             'filename' : 'rr_0.25deg_reg_v10.0.nc'},
                 'temp'   : {'name'     : 'air_temperature',
                             'filename' : 'tg_0.25deg_reg_v10.0.nc'},
                 'mslp'   : {'name'     : 'air_pressure_at_sea_level',
                             'filename' : 'pp_0.25deg_reg_v10.0.nc'}}

era_var_dict = {'temp' : {'name'     : 'air_temperature',
                          'filename' : 't2m'},
                'mslp' : {'name'     : 'air_pressure_at_sea_level',
                          'filename' : 'mslp'}}

def sort_months(months, months_are_bounds):
    """
    Return a sorted month list accounting for overlapping years.
    
    """
    if sorted(months) == months:
        if months_are_bounds:
            months = range(min(months), max(months)+1)
        return months
    else:
        # Months overlap years.
        if months_are_bounds:
            all_months = []
            month = months[0]
            while month != months[-1] + 1:
                all_months.append(month)
                if month == 12:
                    month = 1
                else:
                    month += 1
            months = all_months
        return months

def check(arg, lst):
    if arg in lst:
        return arg
    else:
        raise UserWarning('"%s" is not valid. Valid options are %s' 
                          % (arg, lst))

class TextData(object):
    def analysis(self, method, axis=1):
        """
        Process the data with the appropriate method.
        
        Args:
        
        * method: string
            Name of numpy analysis method, e.g. mean or sum.
        
        * axis: integer
            The axis over which to perform the analysis. Default is 1 which 
            refers to analysis on each year separately.
        
        Returns:
            numpy array
        
        """
        if self.data is not None:
            return getattr(numpy, method.lower())(self.data, axis=axis)

    def climatology(self, clim_period=[1981,2010]):
        """
        Return a class instance for the same variable and season as the current
        instance but for a climatological period.
        
        Args:
        
        * clim_period: list
            All years in between the minimum and maximum year define the
            climatological period.

        Returns:
            numpy array
            
        """
        arg_dict = {}       
        for arg in inspect.getargspec(self.__init__).args:
            if arg not in ['self', 'years', 'years_are_bounds']:
                arg_dict[arg] = getattr(self, arg)   
        arg_dict['years'] = clim_period
        arg_dict['years_are_bounds'] = True
        climatology = self.__class__(**arg_dict)
        return climatology

    def anomalies(self, clim_period=[1981,2010], analysis_method=None):
        """
        Convert the data into anomaly space using its climatology.
        
        Args:
        
        * clim_period: list
            All years in between the minimum and maximum year define the
            climatological period to be removed.
        
        * analysis_method: string
            If yearly analysis has been performed on the data e.g. mean or sum
            of the months, the climatological data must also have yearly 
            analysis performed to match (it does not have to be the same 
            analysis but should be for true anomalies). The string must be the
            name of a numpy analysis method.

        Returns:
            numpy array

        """
        climatology = self.climatology(clim_period)
        climatology.data = climatology.analysis('mean', axis=0)
        if self.data.ndim == 2 and self.data.shape[1] == len(climatology.data):
            return self.data - climatology.data
        elif self.data.ndim == 1:
            if analysis_method:
                climatology.data = climatology.analysis(analysis_method, 
                                                        axis=0)
                return self.data - climatology.data
            else:
                raise UserWarning('Yearly analysis has been performed on the '\
                      'data. Use the analysis_method argument (e.g. mean or '\
                      'sum) to perform analysis on the climatology data so '\
                      'anomalies can be calculated correctly.')
        else:
            raise UserWarning('Analysis has been performed on the data '\
                      'resulting in an unrecognised dimension shape. To avoid'\
                      ' this run the anomalies function before the analysis '\
                      'or use the climatology method and manually remove.')
        

class NAOIndexObs(TextData):
    """
    Load NAO index observations from text file. The year is defined by the 
    first month in the case of months crossing years.
    
    Args:
    
    * months: list
        List of month numbers
    
    * years: list
        List of years
    
    * months_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.
    
    """
    def __init__(self, months=None, years=None, months_are_bounds=True, 
                  years_are_bounds=True):
        if not months:
            months = range(1,13)
        if not years:
            years_are_bounds = True
            years = [1866, datetime.datetime.now().year]
        if max(months) > 12 or min(months) < 1:
            raise ValueError('Months must fall within 1 and 12.')
        if years_are_bounds:
            years = range(min(years), max(years)+1)
        self.years = years
        self.months = sort_months(months, months_are_bounds)
        self.months_are_bounds = months_are_bounds
        self.years_are_bounds  = years_are_bounds
        self._load()
                    
    def _load(self):
        """
        Load data from file.
        
        """
        data = numpy.genfromtxt(NAO_PATH, missing_values=-999, usemask=True)
        ordered_data = []
        for year in self.years:
            year_data = []
            # The first column is the years.
            try:
                year_index = numpy.where(data[:,0] == year)[0][0]
            except IndexError:
                raise UserWarning('There is no data for the year %s. Years '\
                                  'range from 1866 to the current year.' 
                                  % year)
            last_month = self.months[0] - 1
            for month in self.months:
                if month < last_month:
                    year_index += 1
                year_data.append(data[year_index, month])
                last_month = month
            ordered_data.append(year_data)
        # Convert to Pa.
        self.data = numpy.array(ordered_data) * 10.
        

class FloodCountObs(TextData):
    """
    Load flood count observations from text file. The year is defined by the 
    first month in the case of months crossing years.
    
    Args:
    
    * months: list
        List of month numbers
    
    * years: list
        List of years
    
    * months_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.
    
    """
    def __init__(self, months=None, years=None, months_are_bounds=True, 
                  years_are_bounds=True):
        if not months:
            months = range(1,13)
        if not years:
            years_are_bounds = True
            years = [1989, 2012]
        if max(months) > 12 or min(months) < 1:
            raise ValueError('Months must fall within 1 and 12.')
        if years_are_bounds:
            years = range(min(years), max(years)+1)
        self.years = years
        self.months = sort_months(months, months_are_bounds)
        self.months_are_bounds = months_are_bounds
        self.years_are_bounds  = years_are_bounds
        self._load()
                    
    def _load(self):
        """
        Load data from file.
        
        """
        data = numpy.genfromtxt(FLOOD_PATH, dtype=int, missing_values=-999, 
                                usemask=True, skiprows=1)
        ordered_data = []
        for year in self.years:
            year_data = []
            # The first column is the years.
            try:
                year_index = numpy.where(data[:,0] == year)[0][0]
            except IndexError:
                raise UserWarning('There is no data for the year %s. Years '\
                                  'range from 1989 to 2012.' 
                                  % year)
            last_month = self.months[0] - 1
            for month in self.months:
                if month < last_month:
                    year_index += 1
                year_data.append(data[year_index, month])
                last_month = month
            ordered_data.append(year_data)
        self.data = numpy.array(ordered_data)


class NCICTextData(TextData):
    """
    Load data from NCIC text files.

    Args:
    
    * variable: string
        Run NCICTextData.print_variables() available variables.
    
    * months: list
        List of month numbers
    
    * years: list
        List of years
    
    * region: string
        Run NCICTextData.print_regions() for available regions. Default is UK.
    
    * months_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.

    """
    def __init__(self, variable=None, months=None, years=None, region='UK',
                  months_are_bounds=True, years_are_bounds=True, 
                  missing_val=-99999):
        if variable:        
            self.variable = variable.title()
            if not months:
                months = range(1,13)
            if max(months) > 12 or min(months) < 1:
                raise ValueError('Months must fall within 1 and 12.')
            if not years:
                years_are_bounds = True
                years = [1961, datetime.datetime.now().year]
            if years_are_bounds:
                years = range(min(years), max(years)+1)
            self.months = sort_months(months, months_are_bounds)
            self.years  = years
            self.region = region
            self.unit   = ncic_txt_var_dict[self.variable]['units']
            self.months_are_bounds = months_are_bounds
            self.years_are_bounds  = years_are_bounds
            self.missing_val = missing_val
                 
            # Load data.
            var_load_name = ncic_txt_var_dict[self.variable]['name']
            filedata = urlopen('{p}{v}/date/{r}'.format(p=NCICTEXT_PATH, 
                                                        v=var_load_name,
                                                        r=self.region))
            self._load(filedata)
        else:
            self.data = None
        self.variable_list = ncic_txt_var_dict.keys()
        self.days_in_months = None

    def _load(self, filedata):
        """
        Loads data from file taking into account overlapping years.
        
        """
        # First we create an array with all data in it, including +1 year
        # in case of overlaps.
        full_years = range(self.years[0], self.years[-1]+2)
        
        data = []
        for line in filedata:
            line_start = line[0:4]
            if line_start.isdigit() and int(line_start) in full_years:
                line_data = []
                # Index up to 91 as this is where the unwanted season data 
                # starts.
                for val in line[:91].split():
                    try:
                        val = float(val)
                        line_data.append(val)
                    except:
                        line_data.append(self.missing_val)
                # Append the year and the 12 months of data.
                year_data = line_data[:13]
                if len(year_data) != 13:
                    diff = 13 - len(year_data)
                    year_data = year_data + [self.missing_val]*diff
                data.append(year_data)
        data = numpy.ma.masked_values(data, value=self.missing_val)
        
        ordered_data = []
        for year in self.years:
            year_data = []
            # The first column is the years.
            year_index = numpy.where(data[:,0] == year)[0][0]
            last_month = self.months[0] - 1
            for month in self.months:
                if month < last_month:
                    year_index += 1
                year_data.append(data[year_index, month])
                last_month = month
            ordered_data.append(year_data)
        self.data = numpy.array(ordered_data)
    
    @staticmethod
    def print_regions():
        """
        Print the available regions by reading the html. (Not a robust method).
        
        """
        filedata = urlopen('{p}Tmean/date/'.format(p=NCICTEXT_PATH))
        for line in filedata:
            if 'alt="TXT">' in line.split():
                tag = line.split()[4]
                print tag.split('.')[1][1:]

    @staticmethod
    def print_variables():
        """
        Print available variables.
        
        """
        for key in ncic_txt_var_dict.keys():
            print key

    def monthly_totals_to_daily_means(self):
        """
        Convert monthly total data into daily means.
        
        """
        if ncic_txt_var_dict[self.variable]['type'] == 'mean':
            print '%s data is in monthly means, therefore daily means are the'\
                  ' same.' % self.variable
            return self.data
        if self.data.shape != (len(self.years), len(self.months)):
            raise UserWarning('Conversion cannot be made on data which has '\
                              'been processed.')
                
        if self.days_in_months is None:
            all_days_in_months = []
            for year in self.years:
                last_month = None
                days_in_months = []
                for month in self.months:
                    if last_month:
                        if month < last_month:
                            year += 1
                    days_in_months.append(monthrange(year, month)[1])
                    last_month = month
                all_days_in_months.append(days_in_months)
            self.days_in_months = numpy.array(all_days_in_months)

        return self.data / self.days_in_months


class EObsData(ObservationData):
    """
    Sub-class of GetData.ObservationData. Loads E-Obs data. Daily data is 
    available from 1950 up to and including 2012. 

    Args:
    
    * variable: string
        Run EObsData.print_variables() for available variables.
    
    * dates: datetime or list of datetimes
        Specify the dates to be loaded. If dates_are_bounds (see below) is 
        True, all dates in between the min and max dates are loaded.
        Note, If years are given (see below), the year of the dates is ignored. 
        In this case only specify a single period e.g. 5/3/2014 - 10/4/2014, 
        and this will be replicated for all given years (meaning 2014 is not
        necessarily loaded unless specified in years argument). If years are 
        not given then the exact dates are loaded.
        
    * years: list of integers
        Specify the years to be loaded. If years_bounds (see below) is True, 
        all years in between the max and min years are loaded. Default is None
        meaning only the dates given (see above) are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is None.
    
    * dates_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.

    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    """
    def __init__(self, variable, dates, years=None, area_bounds=None,
                  dates_are_bounds=True, years_are_bounds=True,
                  area_bounds_format=['x_min','x_max','y_min','y_max']):
        
        if type(dates) == datetime.datetime:
            dates = [dates, dates]
        variable = variable.lower()
        super(EObsData, self).__init__(
                              directory=EObs_PATH, 
                              variable=eobs_var_dict[variable]['name'],
                              dates=dates, 
                              years=years,
                              area_bounds=area_bounds,
                              dates_are_bounds=dates_are_bounds,
                              years_are_bounds=years_are_bounds,
                              area_bounds_format=area_bounds_format)
        self.load(eobs_var_dict[variable]['filename'], current_dir=True)
        self.variable = variable

    @staticmethod
    def print_variables():
        """
        Print available variables.
        
        """
        for key in eobs_var_dict.keys():
            print key


class NCICData(ObservationData):
    """
    Sub-class of GetData.ObservationData. Loads NCIC spatial data. Monthly data 
    is available from 1910 up to and including 2012.

    Args:
    
    * variable: string
        Run NCICData.print_variables() for available variables.
    
    * months: list
        List of month numbers
        
    * years: list of integers
        Specify the years to be loaded. If years_bounds (see below) is True, 
        all years in between the max and min years are loaded. Default is None
        meaning all years are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is None.
    
    * months_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.

    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    """
    def __init__(self, variable=None, months=None, years=None, 
                  area_bounds=None, months_are_bounds=True, 
                  years_are_bounds=True,
                  area_bounds_format=['x_min','x_max','y_min','y_max']):

        if variable:
            variable = variable.title()
            if not months:
                months = range(1,13)
            if max(months) > 12 or min(months) < 1:
                raise ValueError('Months must fall within 1 and 12.')
            self.months = sort_months(months, months_are_bounds)
            if not years:
                years_are_bounds = True
                years = [1910, 2012]
            self.years = years
            
            dates = []
            this_month = self.months[0]
            year = self.years[0]
            for month in self.months:
                if month < this_month:
                    year += 1
                dates.append(datetime.datetime(year, month, 15))
                this_month = month
                
            super(NCICData, self).__init__(
                                  directory=NCIC_PATH, 
                                  variable=ncic_var_dict[variable]['name'],
                                  dates=dates, 
                                  years=years,
                                  area_bounds=area_bounds,
                                  dates_are_bounds=months_are_bounds,
                                  years_are_bounds=years_are_bounds,
                                  area_bounds_format=area_bounds_format,
                                  multiprocess_workers=None)
            self.load(ncic_var_dict[variable]['filename'], current_dir=True)
            self.variable = variable

    @staticmethod
    def print_variables():
        """
        Print available variables.
        
        """
        for key in ncic_var_dict.keys():
            print key

class ERAIData(ObservationData):
    """
    Sub-class of GetData.ObservationData. Loads ERA-Interim data. Monthly data 
    is available from 1979 up to 2012.

    Args:
    
    * variable: string
        Run ERAIData.print_variables() for available variables.
    
    * dates: datetime or list of datetimes
        Specify the dates to be loaded. If dates_are_bounds (see below) is 
        True, all dates in between the min and max dates are loaded.
        Note, If years are given (see below), the year of the dates is ignored. 
        In this case only specify a single period e.g. 5/3/2014 - 10/4/2014, 
        and this will be replicated for all given years (meaning 2014 is not
        necessarily loaded unless specified in years argument). If years are 
        not given then the exact dates are loaded.
        
    * years: list of integers
        Specify the years to be loaded. If years_bounds (see below) is True, 
        all years in between the max and min years are loaded. Default is None
        meaning only the dates given (see above) are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is None.
    
    * dates_are_bounds: boolean
        If True, all month numbers between the highest and lowest values given
        in months are taken. Default True.
    
    * years_are_bounds: boolean
        If True, all years between the highest and lowest values given in years
        are taken. Default True.

    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    """
    def __init__(self, variable, dates, years=None, area_bounds=None,
                  dates_are_bounds=True, years_are_bounds=True,
                  area_bounds_format=['x_min','x_max','y_min','y_max'],
                  area_bounds_range=[-180,180,-90,90]):
        
        if type(dates) == datetime.datetime:
            dates = [dates, dates]
        variable = variable.lower()
        super(ERAIData, self).__init__(
                              directory=ERA_I_PATH, 
                              variable=era_var_dict[variable]['name'],
                              dates=dates, 
                              years=years,
                              area_bounds=area_bounds,
                              dates_are_bounds=dates_are_bounds,
                              years_are_bounds=years_are_bounds,
                              area_bounds_format=area_bounds_format,
                              area_bounds_range=area_bounds_range)
        filename = 'erai_{v}_mnthmn_Jan79Mar14.grb'.format(v=era_var_dict[variable]['filename'])
        self.load(filename, current_dir=True)
        self.variable = variable

    @staticmethod
    def print_variables():
        """
        Print available variables.
        
        """
        for key in era_var_dict.keys():
            print key
            

class IssuedForecastData(object):
    """
    Class for retrieving forecast and observation data for the coming forecast 
    period.
    
    Args:
    
    * variable: string 
        Specify the variable, 't2m' for air temperature at 2 metres or 'precip'
        for precipitation.
     
    * period: string:
        Specify the period type of the forecast, either 'mon' for the monthly 
        or 'seas' for the seasonal (3 months).
    
    * iss_month: string
        Specify the month the forecast was issued (not the month of the 
        forecast), this is the month before the first forecasted month. Use the
        first 3 letters of the month e.g. 'Jan' or 'Feb'.
    
    * iss_year: integer
        Specify the year the forecast was issued, this is the year of the 
        iss_month, not the year of the forecast period (if different).
    
    Kwargs:
    
    * missing_val: float or integer
         Specify the value used if data is missing.
         
    """
    def __init__(self, variable, period, iss_month, iss_year, 
                  data_dir=ISSFCST_PATH, missing_val=99999.):
        self.variable  = check(variable.lower(), ['t2m', 'precip'])
        self.period    = check(period.lower(), ['seas', 'mon'])
        self.iss_month = check(iss_month.title(), MONS)
        self.iss_year  = iss_year
        self.data_dir  = data_dir
        self.missing_val = missing_val
            
    def _create_filename(self, modified):
        if modified:
            adjusted = '_adj%s' % self.period
        else:
            adjusted = ''
        return '{P}{M}{Y}_{V}{X}.dat'.format(P=self.data_dir,
                                             M=self.iss_month,
                                             Y=str(self.iss_year),
                                             V=self.variable,
                                             X=adjusted)
    
    def _get_index_of_split(self, filename):
        """
        Read through the file to establish where monthly data ends and seasonal
        data begins. This break is where the line starts with a letter instead 
        of a number. Return the index of the first data value after the split.
        Used for raw data only.
        
        """
        text_line_numbers = []
        try:
            with open(filename) as open_file:
                for i, line in enumerate(open_file):
                    if line[0].isalpha():
                        text_line_numbers.append(i)
        except IOError:
            raise IOError("Can't find raw forecast data issued %s %s" 
                          % (self.iss_month, self.iss_year))
        # The index required is at the text line furthest down the file. The 
        # text lines are later removed therefore the required index is the last
        # text line number minus the number of other text lines in the file.
        return text_line_numbers[-1] - len(text_line_numbers[:-1])

    def _seperate_data(self, all_data, split_index):
        """
        Separate the raw data into monthly and seasonal data using the given
        index.
        
        """
        monthly  = all_data[:split_index]
        seasonal = all_data[split_index:]
        monthly  = filter(lambda val: val != self.missing_val, monthly)
        seasonal = filter(lambda val: val != self.missing_val, seasonal)
        return monthly, seasonal
    
    def _raw_data_load(self, col=4):
        """
        Both monthly and seasonal data is contained in the same file. Load the
        the whole file and extract the data for the required period.
        The col is the column within the file, default is the column containing
        forecast data. The file also contains climatology data.
        
        """
        filename = self._create_filename(modified=False)
        # Get the index (line number) which splits the monthly and seasonal 
        # data.
        split_index = self._get_index_of_split(filename)
        all_data = numpy.genfromtxt(filename, delimiter='\t', usecols=col, 
                                    invalid_raise=False, filling_values=self.missing_val)
        monthly, seasonal = self._seperate_data(all_data, split_index)
        if self.period == 'mon':
            return monthly
        else:
            return seasonal
    
    def _mod_data_load(self, col=3):
        """
        Load the modified data.
        
        """
        filename = self._create_filename(modified=True)
        try:
            data = numpy.genfromtxt(filename, delimiter='\t', usecols=col, 
                                    skiprows=2, 
                                    filling_values=self.missing_val)
        except IOError:
            raise IOError("Can't find modified forecast data issued %s %s" 
                           % (self.iss_month, self.iss_year))
        return filter(lambda val: val != self.missing_val, data)

    def _obs_from_file_load(self):
        """
        Load the observation data from the raw data file.
        
        """
        return self._raw_data_load(col=2)

    def _member_numbers_load(self):
        """
        Load the observation data from the raw data file.
        
        """
        return self._raw_data_load(col=3)

    def _obs_from_ncic_load(self, years, region):
        """
        Load observation data directly from NCIC pages.
        
        """
        # Set up the variables to load NCIC data.
        if self.variable == 't2m':
            variable = 'temp'
            method   = 'MEAN'
        elif self.variable == 'precip':
            variable = self.variable
            method   = 'SUM'
            
        iss_month_index = MONS.index(self.iss_month)
        if self.period == 'mon':
            # Add 2, 1 to adjust the index to a month number and 1 because we
            # look at month ahead.
            fcst_months = [iss_month_index + 2]
        elif self.period == 'seas':
            fcst_months = range(iss_month_index + 2, iss_month_index + 5)
        # fcst_months can have values over 12, these must be sorted.
        months = []
        for month_num in fcst_months:
            if month_num > 12:
                month_num -= 12
            months.append(month_num)
            
        ncic = NCICTextData(variable, months, years, region=region,
                            months_are_bounds=False)
        obs = ncic.analysis(method)
        for ob, year in zip(obs, years):
            if numpy.isnan(ob):
                raise UserWarning('Observation for month(s) %s in the year '\
                                  '%s, does not exist.' % (months, year))
        return ncic.analysis(method)

    def _get_savename(self, dtype, modified=False):
        """
        Create a filename to save data.
        
        """
        if self.variable == 't2m':
            var = 'temp'
        else:
            var = self.variable
        if dtype == 'mod':
            if modified:
                typ = 'mod'
            else:
                typ = 'raw' 
            return '{P}_{V}_{T}.txt'.format(P=self.period, V=var, T=typ)
        elif dtype == 'obs':
            return 'obs_{P}_{V}.txt'.format(P=self.period, V=var)

    def model_load(self, modified=False, get_member_numbers=False):
        """
        Load the model forecast members.
        
        Kwargs:
        
        * modified: boolean
            Whether to load the modified forecast of original members. Note, 
            modified members only exist after a seasonal meeting.
        
        * get_member_numbers: boolean
            Whether the return the forecast member numbers as well as the 
            members themselves. The members are returned as a seperate list 
            in which each element corresponds to the equivilant elememt in the
            returned data. Note, this is only available for unmodified data.
        
        Returns:
            numpy array (a seperate array is also returned if 
            get_member_numbers is set to True)
        
        """
        if modified:
            data = self._mod_data_load()
        else:
            data = self._raw_data_load()
            if get_member_numbers:
                mems = self._member_numbers_load()
                return data, mems
        return data
    
    def current_obs_load(self, region='UK'):
        """
        Return the observation for the forecast period.
        
        Kwargs:
        
        * region: string
            Region with the UK defined by NCIC. Default is UK.
        
        Returns:
            float
            
        """
        # Make sure the correct year is used when loading the current 
        # observation.
        if self.iss_month == "Dec":
            year = self.iss_year + 1
        else:
            year = self.iss_year
        return self._obs_from_ncic_load([year], region)[0]

    def climatology_obs_load(self, source='ncic', clim_period=[1981, 2010], 
                               region='UK'):
        """
        Return the climatology for the forecast period.
        
        Kwargs:
        
        * source: 'file' or 'ncic'
            Data from 'file' does come from NCIC indirectly but is guaranteed 
            to be consistent (where any small adjustments are made). 'ncic' 
            takes data directly from the NCIC pages. Note, the kwargs 
            clim_period and region are ignored if source is set to 'file'.
        
        * clim_period: list of 2 years
            Specify the climatological period.
        
        * region: string
            Region with the UK defined by NCIC. Default is UK.
        
        
        Returns:
            numpy array
        
        """
        if source == 'ncic':
            return self._obs_from_ncic_load(clim_period, region)
        else:
            return self._obs_from_file_load()
        
    def print_data(self, data, dtype, modified=False):
        """
        Used to print current data and all previous saved data in the save 
        file. Note, this is useful specifically for the operational 
        verification work.
        
        Args:
        
        * data: array like
            Normally an array produced by this class
        
        * dtype: string 'model' or 'obs'
            Specifies the type of data, i.e. which file to look in.
        
        Kwargs:
        
        * modified: boolean
            If looking a past forecast data, specify whether to look at data
            before or after modification.
        
        """
        print 'Data in saved file:\n'
        data_fin = self.save_data(data, dtype, modified, save=False, 
                                  check_exists=False, print_data=True)            
        print '\n\nCurrent loaded data:\n'
        print data_fin
    
    def save_data(self, data, dtype, modified=False, save=False, 
                   check_exists=True, print_data=False):
        """
        Fill any missing members with the set missing_val, check the file has 
        not been updated already and save data to file. Note, this useful 
        specifically for the operational verification work.
        
        Args:
        
        * data: array like
            Normally an array produced by this class
        
        * dtype: string 'model' or 'obs'
            Specifies the type of data, i.e. which file to look in.
        
        Kwargs:
        
        * modified: boolean
            If looking a past forecast data, specify whether to look at data
            before or after modification.
        
        * save: boolean
            Whether to actually make the save.
        
        * check_exists: boolean
            Read through the save file looking for data which matches the 
            current data. This just checks you are not saving the same data
            twice as this causes problems in verification work.
        
        * print_data: boolean
            Prints the data in the file. Use print_data method to print this
            and current data.
        
        """
        dtype = check(dtype, ['model', 'obs'])
        if dtype == 'model':
            filename = self._get_savename('mod', modified)
            while len(data) < 42:
                data.append(self.missing_val)
            data_fin = ' '.join([str(x) for x in data])
        
        elif dtype == 'obs':
            filename = self._get_savename('obs')
            data_fin = str(data)
        
        with open(ISS_SAVEFILE+filename, 'r+') as outfile:
            curr_data = outfile.readlines()
            while curr_data[-1] == '\n':
                del curr_data[-1]
            
            if print_data:
                for prev_data in curr_data:
                    print prev_data,
                return data_fin
            
            if check_exists:
                unique = True
                for prev_data in curr_data:
                    if dtype == 'model':
                        if prev_data.split()[0:42] == data_fin.split()[0:42]:
                            print 'This data, from forecast issued %s %s, has been saved.'\
                                  % (self.iss_month, self.iss_year)
                            unique = False
                    if dtype == 'obs':
                        if prev_data.strip() == data_fin:
                            print 'This observation value has been saved, however, it may be another period with the same value.'
                            unique = False
                if unique:
                    print 'Unsaved %s data' % dtype
            if save:
                outfile.write('\n'+data_fin)
                