
def test_compare_meridian_constraint_and_observer():
    
    vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg),
                   name="Vega")
    rigel = FixedTarget(coord=SkyCoord(ra=78.63446707*u.deg, dec=8.20163837*u.deg),
                    name="Rigel")
    polaris = FixedTarget(coord=SkyCoord(ra=37.95456067*u.deg,
                                     dec=89.26410897*u.deg), name="Polaris")
    time = Time('2001-02-03 04:05:06')
    time_ranges = [Time([time, time+1*u.hour]) + offset
                   for offset in np.arange(0, 400, 100)*u.day]
    for time_range in time_ranges:
        subaru = Observer.at_site("Subaru")
        targets = [vega, rigel, polaris]
        hours2seconds = 60 * 60
        max_hrs = 3.5*u.hour
        # Check if each target meets meridian constraints using Observer
        always_from_observer = [all([abs(time - subaru.target_meridian_transit_time(time, target)) < max_hrs * hours2seconds
                                     for time in time_grid_from_range(time_range)])
                                for target in targets]
    
        # Check if each target meets meridian constraints using
        # is_always_observable and MeridianConstraint
        always_from_constraint = is_always_observable(MeridianConstraint(max=max_hrs),
                                                      subaru, targets,
                                                      time_range=time_range)
        assert all(always_from_observer == always_from_constraint)
