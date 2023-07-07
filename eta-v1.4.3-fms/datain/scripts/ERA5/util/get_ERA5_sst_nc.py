import cdsapi
import sys
yeari = sys.argv[1]
monthi = sys.argv[2]
dayi = sys.argv[3]
timei = sys.argv[4]
c = cdsapi.Client()

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
        'area': [
            35, -150, -75,
            35,
        ],
        'format': 'netcdf',
    },
    'ERA5_SST_'+yeari+monthi+dayi+timei+'.nc')
