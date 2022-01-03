import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            'divergence', 'relative_humidity', 'specific_humidity',
            'temperature', 'u_component_of_wind', 'v_component_of_wind',
            'vertical_velocity', 'vorticity',
        ],
        'pressure_level': [
            '200', '225', '250',
            '300', '350', '400',
            '450', '500', '550',
            '600', '650', '700',
        ],
        'year': '2021',
        'month': '11',
        'time': '00:00',
        'area': [
            30, -120, -60,
            -30,
        ],
    },
    'nov21.nc')
