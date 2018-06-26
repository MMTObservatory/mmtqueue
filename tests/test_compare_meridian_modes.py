    
    
def test_compare_meridian_modes():
    
    vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg),
                   name="Vega")
    rigel = FixedTarget(coord=SkyCoord(ra=78.63446707*u.deg, dec=8.20163837*u.deg),
                    name="Rigel")
    polaris = FixedTarget(coord=SkyCoord(ra=37.95456067*u.deg,
                                     dec=89.26410897*u.deg), name="Polaris")
    time = Time('2001-02-03 04:05:06')
    time_ranges = [Time([time, time+1*u.hour]) + offset
                   for offset in np.arange(0, 400, 100)*u.day]
    
    subaru = Observer.at_site("Subaru")
    targets = [vega, rigel, polaris]
    max_hrs = 3.5*u.hour
    for time_range in time_ranges:
         always_with_python = is_always_observable(MeridianConstraint(max=max_hrs,mode='python'),
                                                      subaru, targets,
                                                      time_range=time_range)
         always_with_numpy = is_always_observable(MeridianConstraint(max=max_hrs,mode='numpy'),
                                                      subaru, targets,
                                                      time_range=time_range)
         assert all(always_with_python == always_with_numpy)

