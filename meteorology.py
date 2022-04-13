import ee
import math

GSC = 0.0820 # Solar constant MJ m2
DEG2RAD = math.pi / 180.0 # Degree to rad

def meteorology(time_start, meteo_source):

    meteo_data = ee.ImageCollection(meteo_source)\
                .filterDate(ee.Date(time_start),ee.Date(time_start).advance(1,'day'))

    if 'gldas' in meteo_source.lower():

        bands = {
            'shortwave_radiation':'SWdown_f_tavg',
            'specific_humidity': 'Specific humidity',
            'temperature': 'Tair_f_inst',
            'ux': 'Wind_f_inst',
            'pressure': 'Pressure'
        }

    elif 'era5' in meteo_source.lower():

        bands = {
            'shortwave_radiation':'surface_solar_radiation_downwards_hourly',
            'dew_point_temp': 'dewpoint_temperature_2m',
            'temperature': 'temperature_2m',
            'u_wind': 'u_component_of_wind_10m',
            'v_wind': 'v_component_of_wind_10m',
            'pressure': 'surface_pressure'
        }

    else:
        raise ValueError('Meteorology dataset not supported.')

    #wind speed
    if 'u_wind' in bands:

        wind_u = meteo_data.select(bands['u_wind']).mean()
        wind_v = meteo_data.select(bands['v_wind']).mean()

        wind_med = wind_u.expression(
            'sqrt(ux_u ** 2 + ux_v ** 2)', {'ux_u': wind_u, 'ux_v': wind_v},
        ).rename('ux')
    
    else:

        wind_med = meteo_data.select(bands['ux']).mean().rename('ux')

    # tmin and tmax
    tmin = meteo_data.select(bands['temperature']).min().rename('tmin')
    tmax = meteo_data.select(bands['temperature']).max().rename('tmax')

    # daily shortwave
    if 'gldas' in meteo_source.lower():

        swdown =  ee.ImageCollection(meteo_source)\
                .filterDate(ee.Date(time_start).advance(3,'hour'),
                            ee.Date(time_start).advance(27,'hour'))\
                .select(bands['shortwave_radiation'])\
                .mean()\
                .rename('shortwave_radiation')

    else: 

        swdown =  ee.ImageCollection(meteo_source)\
                .filterDate(ee.Date(time_start),
                            ee.Date(time_start).advance(1,'day'))\
                .select("surface_solar_radiation_downwards_hourly")\
                .sum()\
                .divide(86400)\
                .rename('shortwave_radiation')

    # temperature
    tair = meteo_data.select(bands['temperature']).mean()


    return [tmin, tmax, tair, swdown, wind_med]



def daily_rad(time_start, tmax, tmin, elev_source, rso24h):

    rs = rso24h.multiply(0.0864).rename('Rs')

    doy = ee.Date(time_start).getRelative('day', 'year').add(1)

    dr = rso24h.expression(
        '1 + (0.033 * cos((2 * pi / 365) * doy))', {'doy': doy, 'pi': math.pi})

    sd = rso24h.expression(
        '0.40928 * sin(((2 * pi / 365) * doy) - 1.39)', {'doy': doy, 'pi': math.pi})

    lat = rso24h.pixelLonLat().select(['latitude']).multiply(DEG2RAD)\
        .rename('latitude')

    ws = rso24h.expression('acos(-tan(Lat) * tan(Sd))', {'Lat': lat, 'Sd': sd})

    rad_a = rso24h.expression(
        'Ws * sin(Lat) * sin(Sd) + cos(Lat) * cos(Sd) * sin(Ws)',
        {'Ws': ws, 'Lat': lat, 'Sd': sd},
    )

    ra = rso24h.expression(
        '((24 * 60) / pi) * Gsc * Dr * rad_a',
        {'pi': math.pi, 'Gsc': GSC, 'Dr': dr, 'rad_a': rad_a},
    )

    rso = rso24h.expression(
        '(0.75 + 2E-5 * z) * Ra', {'z': elev_source, 'Ra': ra})

    rns = rso24h.expression(
        '(1 - albedo) * Rs', {'Rs': rs, 'albedo': 0.23})

    ea = rso24h.expression(
        '0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))', {'T_air': tmin.subtract(273.15)})

    rnl = rso24h.expression(
        '4.901E-9 * ((Tmax ** 4 + Tmin ** 4) / 2) * (0.34 - 0.14 * sqrt(ea)) * '
        '(1.35 * (Rs / Rso) - 0.35)',
        {'Tmax': tmax, 'Tmin': tmin, 'ea': ea,
         'Rs': rs, 'Rso': rso},
    )

    rn = rso24h.expression('Rns - Rnl', {'Rns': rns, 'Rnl': rnl})

    return rn.rename('rad_24h')


def eto_grass(time_start, tmin, tmax, tair, ws, rso24h, elev_source):

    rad24h = daily_rad(time_start,tmax, tmin, elev_source, rso24h)

    press = tmin.expression(
    '101.3 * pow(((293 - (0.0065 * Z))/ 293),5.26) ', {
       'Z' : elev_source,
       }).rename('p_atm')
    
    const_psci=tmin.expression(
            '0.000665* P ',{
                    'P':press}).rename('cte_psi')
    
    ea_max = tmax.expression(
     ' 0.6108 *(exp( (17.27 * tmax) / (tmax + 237.3)))', {
       'tmax': tmax.subtract(273.15),  
       }).rename('ea_max')
    
    ea_min = tmin.expression(
     ' 0.6108 *(exp( (17.27 * tmin) / (tmin + 237.3)))', {
       'tmin': tmin.subtract(273.15),  
       }).rename('ea_min')

    es_mean = tmin.expression(
     '(ea_min +ea_max)/2', {
       'ea_min': ea_min,  
       'ea_max':ea_max
     }).rename('es_mean')

    delta= tmin.expression(
            '(4098*es)/(tair)**2',{
            'es': es_mean,
            'tair':tair.subtract(273.15)}).rename('delta')

    eto24h = tmin.expression(
            '(0.408*delta*(Rn)+ cte_psi*(Cn/(Tair + 273))*ws*(es -ea))/'+
                                          '(delta + cte_psi*(1 + Cd*ws))',{
            'delta':delta,
            'Rn':rad24h,
            'cte_psi':const_psci,
            'Cn':ee.Number(900),
            'Cd':ee.Number(0.34),
            'Tair':tair.subtract(273.15),
            'ws':ws,
            'es':es_mean,
            'ea':ea_min}).rename('eto24h')

    return eto24h

