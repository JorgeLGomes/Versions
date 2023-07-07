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

        'variable': [
            'geopotential', 'relative_humidity', 'specific_humidity',
            'temperature', 'u_component_of_wind', 'v_component_of_wind',
        ],
        'pressure_level': [
            '1000','975','950','925','900','875','850','825','800','775',
            '750','700','650','600','550','500','450','400','350','300',
            '250','225','200','175','150','125','100','70','50','30', '20',
            '10','7','5','3','2','1',
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
         'format': 'netcdf',
    },
    'ERA5_Pressure_'+yeari+monthi+dayi+timei+'.nc')

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',

        'variable': [
            'ozone_mass_mixing_ratio',
        ],
        'pressure_level': [
            '1000','975','950','925','900','875','850','825','800','775',
            '750','700','650','600','550','500','450','400','350','300',
            '250','225','200','175','150','125','100','70','50','30', '20',
            '10','7','5','3','2','1',
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
         'format': 'netcdf',
    },
    'ERA5_Ozone_'+yeari+monthi+dayi+timei+'.nc')


c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',

        'variable': [
            'land_sea_mask','mean_sea_level_pressure','surface_pressure',
            'geopotential', 
            'soil_temperature_level_1','volumetric_soil_water_layer_1',
            'soil_temperature_level_2','volumetric_soil_water_layer_2',
            'soil_temperature_level_3','volumetric_soil_water_layer_3',
            'soil_temperature_level_4','volumetric_soil_water_layer_4',
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
         'format': 'netcdf',
    },
    'ERA5_Surface_'+yeari+monthi+dayi+timei+'.nc')

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            'sea_surface_temperature', 
        ],
        'year': yeari,
        'month': monthi,
        'day': dayi,
        'time': timei+':00',
        'format': 'netcdf',
    },
    'ERA5_SST_'+yeari+monthi+dayi+timei+'.nc')
