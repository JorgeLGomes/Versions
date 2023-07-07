import cdsapi
import sys
yeari = sys.argv[1]
monthi = sys.argv[2]
dayi = sys.argv[3]
timei = sys.argv[4]
c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': [
            'geopotential', 'relative_humidity', 'specific_humidity',
            'temperature', 'u_component_of_wind', 'v_component_of_wind',
        ],
        'pressure_level': [
            '20', '30', '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '400', '450',
            '500', '550', '600',
            '650', '700', '750',
            '775', '800', '825',
            '850', '875', '900',
            '925', '950', '975',
            '1000',
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
        'area': [
            30, -150, -70,
            30,
        ],
    },
    'ERA5_Pressure_'+yeari+monthi+dayi+timei+'.grib')


c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': [
            'land_sea_mask', 'mean_sea_level_pressure','surface_pressure', 'geopotential',
            'soil_temperature_level_1','volumetric_soil_water_layer_1',
            'soil_temperature_level_2','volumetric_soil_water_layer_2',
            'soil_temperature_level_3','volumetric_soil_water_layer_3',
            'soil_temperature_level_4','volumetric_soil_water_layer_4',
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
        'area': [
            30, -150, -70,
            30,
        ],
    },
    'ERA5_Surface_'+yeari+monthi+dayi+timei+'.grib')

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': [
            'sea_surface_temperature', 
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
        'area': [
            30, -150, -70,
            30,
        ],
    },
    'ERA5_SST_'+yeari+monthi+dayi+timei+'.grib')
