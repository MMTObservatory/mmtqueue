# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specify and constraints to determine which targets are observable for
an observer.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
from abc import ABCMeta, abstractmethod
import datetime
import warnings
import traceback

# Third-party
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_body, get_sun, get_moon, SkyCoord
from astropy import table

import sqlite3
import numpy as np
from numpy.lib.stride_tricks import as_strided

# Package
from mmtqueue.moon import moon_illumination
from mmtqueue.utils import time_grid_from_range
from mmtqueue.target import get_skycoord

# from mmtqueue.utilities import get_constraint_name
from mmtqueue.utilities import get_key
from mmtqueue.utilities import get_score
from mmtqueue.utilities import set_score
from mmtqueue.utilities import localize_time
from mmtqueue.utilities import dict_factory
from mmtqueue.utilities import get_target_date_key
from mmtqueue.utilities import get_sun_moon_date_key

# Required by new constraints
from astropy.coordinates import Angle
from astropy.time import TimeDelta

# This __all__ contains new and modified constraints from astroplan 0.4.
__all__ = ["AltitudeConstraint", "AirmassConstraint", "AtNightConstraint",
           "is_observable", "is_always_observable", "time_grid_from_range",
           "SunSeparationConstraint", "MoonSeparationConstraint",
           "MoonIlluminationConstraint", "LocalTimeConstraint",
           "PrimaryEclipseConstraint", "SecondaryEclipseConstraint",
           "Constraint", "TimeConstraint", "observability_table",
           "months_observable", "max_best_rescale", "min_best_rescale",
           "PhaseConstraint", "is_event_observable",
           "MaskAngleConstraint", "MaskNumberConstraint",
           "MeridianConstraint", "PIPriorityConstraint",
           "RotatorConstraint", "TimeAllocationConstraint",
           "TransitConstraint", "MaskReadyConstraint",
           "ObjectTypeConstraint", "TargetOfOpportunityConstraint"]

# =============================================================================
# __all__ = ["AltitudeConstraint", "AirmassConstraint", "AtNightConstraint",
#            "is_observable", "is_always_observable", "time_grid_from_range",
#            "SunSeparationConstraint", "MoonSeparationConstraint",
#            "MoonIlluminationConstraint", "LocalTimeConstraint",
#            "PrimaryEclipseConstraint", "SecondaryEclipseConstraint",
#            "Constraint", "TimeConstraint", "observability_table",
#            "months_observable", "max_best_rescale", "min_best_rescale",
#            "PhaseConstraint", "is_event_observable"]
# =============================================================================

use_daily_events = True
sqlite_file = 'mmtqueue.sqlite'
doing_line_profiling = False


def _make_cache_key(times, targets):
    """
    Make a unique key to reference this combination of ``times`` and ``targets``.

    Often, we wish to store expensive calculations for a combination of
    ``targets`` and ``times`` in a cache on an ``observer``` object. This
    routine will provide an appropriate, hashable, key to store these
    calculations in a dictionary.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint.
    targets : `~astropy.coordinates.SkyCoord`
        Target or list of targets.

    Returns
    -------
    cache_key : tuple
        A hashable tuple for use as a cache key
    """
    # make a tuple from times
    try:
        timekey = tuple(times.jd) + times.shape
    except BaseException:        # must be scalar
        timekey = (times.jd,)
    # make hashable thing from targets coords
    try:
        if hasattr(targets, 'frame'):
            # treat as a SkyCoord object. Accessing the longitude
            # attribute of the frame data should be unique and is
            # quicker than accessing the ra attribute.
            targkey = tuple(targets.frame.data.lon.value.ravel()) + targets.shape
        else:
            # assume targets is a string.
            targkey = (targets,)
    except BaseException:
        targkey = (targets.frame.data.lon,)
    return timekey + targkey


def _get_altaz(times, observer, targets, force_zero_pressure=False):
    """
    Calculate alt/az for ``target`` at times linearly spaced between
    the two times in ``time_range`` with grid spacing ``time_resolution``
    for ``observer``.

    Cache the result on the ``observer`` object.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint.
    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets.
    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``.
    force_zero_pressure : bool
        Forcefully use 0 pressure.

    Returns
    -------
    altaz_dict : dict
        Dictionary containing two key-value pairs. (1) 'times' contains the
        times for the alt/az computations, (2) 'altaz' contains the
        corresponding alt/az coordinates at those times.
    """
    if not hasattr(observer, '_altaz_cache'):
        observer._altaz_cache = {}

    # convert times, targets to tuple for hashing
    aakey = _make_cache_key(times, targets)

    if aakey not in observer._altaz_cache:
        try:
            if force_zero_pressure:
                observer_old_pressure = observer.pressure
                observer.pressure = 0

            altaz = observer.altaz(times, targets, grid_times_targets=False)
            observer._altaz_cache[aakey] = dict(times=times,
                                                altaz=altaz)
        finally:
            if force_zero_pressure:
                observer.pressure = observer_old_pressure

    return observer._altaz_cache[aakey]


def _get_moon_data(times, observer, force_zero_pressure=False):
    """
    Calculate moon altitude az and illumination for an array of times for
    ``observer``.

    Cache the result on the ``observer`` object.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint.
    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``.
    force_zero_pressure : bool
        Forcefully use 0 pressure.

    Returns
    -------
    moon_dict : dict
        Dictionary containing three key-value pairs. (1) 'times' contains the
        times for the computations, (2) 'altaz' contains the
        corresponding alt/az coordinates at those times and (3) contains
        the moon illumination for those times.
    """
    if not hasattr(observer, '_moon_cache'):
        observer._moon_cache = {}

    # convert times to tuple for hashing
    aakey = _make_cache_key(times, 'moon')

    if aakey not in observer._moon_cache:
        try:
            if force_zero_pressure:
                observer_old_pressure = observer.pressure
                observer.pressure = 0

            altaz = observer.moon_altaz(times)
            illumination = np.array(moon_illumination(times))
            observer._moon_cache[aakey] = dict(times=times,
                                               illum=illumination,
                                               altaz=altaz)
        finally:
            if force_zero_pressure:
                observer.pressure = observer_old_pressure

    return observer._moon_cache[aakey]


def _get_meridian_transit_times(times, observer, targets):
    """
    Calculate next meridian transit for an array of times for ``targets`` and
    ``observer``.

    Cache the result on the ``observer`` object.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint
    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``
    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    Returns
    -------
    time_dict : dict
        Dictionary containing a key-value pair. 'times' contains the
        meridian_transit times.
    """
    if not hasattr(observer, '_meridian_transit_cache'):
        observer._meridian_transit_cache = {}

    # convert times to tuple for hashing
    aakey = _make_cache_key(times, targets)

    if aakey not in observer._meridian_transit_cache:
        meridian_transit_times = observer.target_meridian_transit_time(times, targets)
        observer._meridian_transit_cache[aakey] = dict(times=meridian_transit_times)

    return observer._meridian_transit_cache[aakey]


@abstractmethod
class Constraint(object):
    """
    Abstract class for objects defining observational constraints.
    """
    __metaclass__ = ABCMeta

    def __call__(self, observer, targets, times=None,
                 time_range=None, time_grid_resolution=0.5 * u.hour,
                 grid_times_targets=False):
        """
        Compute the constraint for this class

        Parameters
        ----------
        observer : `~astroplan.Observer`
            the observation location from which to apply the constraints
        targets : sequence of `~astroplan.Target`
            The targets on which to apply the constraints.
        times : `~astropy.time.Time`
            The times to compute the constraint.
            WHAT HAPPENS WHEN BOTH TIMES AND TIME_RANGE ARE SET?
        time_range : `~astropy.time.Time` (length = 2)
            Lower and upper bounds on time sequence.
        time_grid_resolution : `~astropy.units.quantity`
            Time-grid spacing
        grid_times_targets : bool
            if True, grids the constraint result with targets along the first
            index and times along the second. Otherwise, we rely on broadcasting
            the shapes together using standard numpy rules.
        Returns
        -------
        constraint_result : 1D or 2D array of float or bool
            The constraints. If 2D with targets along the first index and times along
            the second.
        """

        if times is None and time_range is not None:
            times = time_grid_from_range(time_range,
                                         time_resolution=time_grid_resolution)

        if grid_times_targets:
            targets = get_skycoord(targets)
            # TODO: these broadcasting operations are relatively slow
            # but there is potential for huge speedup if the end user
            # disables gridding and re-shapes the coords themselves
            # prior to evaluating multiple constraints.
            if targets.isscalar:
                # ensure we have a (1, 1) shape coord
                targets = SkyCoord(np.tile(targets, 1))[:, np.newaxis]
            else:
                targets = targets[..., np.newaxis]
        times, targets = observer._preprocess_inputs(times, targets, grid_times_targets=False)
        result = self.compute_constraint(times, observer, targets)

        # make sure the output has the same shape as would result from
        # broadcasting times and targets against each other
        if targets is not None:
            # broadcasting times v targets is slow due to
            # complex nature of these objects. We make
            # to simple numpy arrays of the same shape and
            # broadcast these to find the correct shape
            shp1, shp2 = times.shape, targets.shape
            x = np.array([1])
            a = as_strided(x, shape=shp1, strides=[0] * len(shp1))
            b = as_strided(x, shape=shp2, strides=[0] * len(shp2))
            output_shape = np.broadcast(a, b).shape
            if output_shape != np.array(result).shape:
                result = np.broadcast_to(result, output_shape)

        return result

    @abstractmethod
    def compute_constraint(self, times, observer, targets):
        """
        Actually do the real work of computing the constraint.  Subclasses
        override this.

        Parameters
        ----------
        times : `~astropy.time.Time`
            The times to compute the constraint
        observer : `~astroplan.Observer`
            the observaton location from which to apply the constraints
        targets : sequence of `~astroplan.Target`
            The targets on which to apply the constraints.

        Returns
        -------
        constraint_result : 2D array of float or bool
            The constraints, with targets along the first index and times along
            the second.
        """
        # Should be implemented on each subclass of Constraint
        raise NotImplementedError

    def process_comment(self, targets, times, score, comment):
        if len(comment) > 0:
            try:
                if targets.isscalar:
                    target = targets
                else:
                    target = targets[0]

                if target.isscalar:
                    pass
                else:
                    target = target[0]

                if times.isscalar:
                    time = times
                else:
                    time = times[0]

                if time.isscalar:
                    pass
                else:
                    time = time[0]

                key = get_key(self, target, time)
                # create a database connection
                conn = sqlite3.connect(sqlite_file)
                schedule_time = localize_time(time)
                set_score(conn, key, score, time=schedule_time, comment=comment)
            except:
                print("Traceback from process_comment: ", traceback.format_exc())
        return score
#
# This constraint is unchanged from astroplan 0.4 code
#


class AltitudeConstraint(Constraint):
    """
    Constrain the altitude of the target.
    .. note::
        This can misbehave if you try to constrain negative altitudes, as
        the `~astropy.coordinates.AltAz` frame tends to mishandle negative
    Parameters
    ----------
    min : `~astropy.units.Quantity` or `None`
        Minimum altitude of the target (inclusive). `None` indicates no limit.
    max : `~astropy.units.Quantity` or `None`
        Maximum altitude of the target (inclusive). `None` indicates no limit.
    boolean_constraint : bool
        If True, the constraint is treated as a boolean (True for within the
        limits and False for outside).  If False, the constraint returns a
        float on [0, 1], where 0 is the min altitude and 1 is the max.
    """

    def __init__(self, min=None, max=None, boolean_constraint=True, debug=False):
        if min is None:
            self.min = -90 * u.deg
        else:
            self.min = min
        if max is None:
            self.max = 90 * u.deg
        else:
            self.max = max

        self.boolean_constraint = boolean_constraint
        self.debug = debug
        self.comments = ''

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.
         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~astroplan.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
          The targets on which to apply the constraints.
         Returns
         -------
         (scores) : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """
        cached_altaz = _get_altaz(times, observer, targets)
        alt = cached_altaz['altaz'].alt
        if self.boolean_constraint:
            lowermask = self.min <= alt
            uppermask = alt <= self.max

            # New code for logging comments for special cases.
            try:
                lower_value = np.prod(lowermask)
                upper_value = np.prod(uppermask)
                comment = ''
                if lower_value == 0.0:
                    # comment = 'Target below minimum altitude: ', repr(alt)
                    comment = 'Target below minimum altitude'
                    score = lower_value
                elif upper_value == 0:
                    # comment = 'Target above maximum altitude: ', repr(alt)
                    comment = 'Target above maximum altitude'
                    score = upper_value
                if len(comment) > 0:
                    self.process_comment(targets, times, score, comment)

            except:
                print("traceback: ", traceback.format_exc())

            return lowermask & uppermask
        else:
            return max_best_rescale(alt, self.min, self.max)
#
# This constraint is unchanged from astroplan 0.4 code
#


class AirmassConstraint(AltitudeConstraint):
    """
    Constrain the airmass of a target.

    In the current implementation the airmass is approximated by the secant of
    the zenith angle.

    .. note::
        The ``max`` and ``min`` arguments appear in the order (max, min)
        in this initializer to support the common case for users who care
        about the upper limit on the airmass (``max``) and not the lower
        limit.

    Parameters
    ----------
    max : float or `None`
        Maximum airmass of the target. `None` indicates no limit.
    min : float or `None`
        Minimum airmass of the target. `None` indicates no limit.
    boolean_contstraint : bool

    Examples
    --------
    To create a constraint that requires the airmass be "better than 2",
    i.e. at a higher altitude than airmass=2::

        AirmassConstraint(2)
    """

    def __init__(self, max=None, min=1, boolean_constraint=True):
        self.min = min
        self.max = max
        self.boolean_constraint = boolean_constraint

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~astroplan.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~mmtqueue.Target`
          The targets on which to apply the constraints.

         Returns
         -------
         (scores) : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """

        cached_altaz = _get_altaz(times, observer, targets)
        secz = cached_altaz['altaz'].secz.value
        if self.boolean_constraint:
            if self.min is None and self.max is not None:
                mask = secz <= self.max
                upper_mask = mask
                lower_mask = False
            elif self.max is None and self.min is not None:
                mask = self.min <= secz
                lower_mask = mask
                upper_mask = False
            elif self.min is not None and self.max is not None:
                # mask = (self.min <= secz) & (secz <= self.max)
                lower_mask = self.min <= secz
                upper_mask = secz <= self.max
                mask = lower_mask & upper_mask
            else:
                raise ValueError("No max and/or min specified in "
                                 "AirmassConstraint.")

            # New code for logging comments for special cases.
            try:
                lower_value = np.prod(lower_mask)
                upper_value = np.prod(upper_mask)
                comment = ''
                if lower_value == 0.0:
                    # comment = 'Target below lower airmass limit, airmass: ', str(secz)
                    comment = 'Target below lower airmass limit'
                    score = lower_value
                elif upper_value == 0:
                    # comment = 'Target above upper airmass limit, airmass: ', str(secz)
                    comment = 'Target above upper airmass limit'
                    score = upper_value
                if len(comment) > 0:
                    self.process_comment(targets, times, score, comment)
            except:
                print("traceback: ", traceback.format_exc())

            return mask
        else:
            if self.max is None:
                raise ValueError("Cannot have a float AirmassConstraint if max is None.")
            else:
                mx = self.max

            mi = 1 if self.min is None else self.min
            # values below 1 should be disregarded
            return min_best_rescale(secz, mi, mx, less_than_min=0)

#
# This constraint is unchanged from astroplan 0.4 code
#


class AtNightConstraint(Constraint):
    """
    Constrain the Sun to be below ``horizon``.
    """
    @u.quantity_input(horizon=u.deg)
    def __init__(self, max_solar_altitude=0 * u.deg, force_pressure_zero=True):
        """
        Parameters
        ----------
        max_solar_altitude : `~astropy.units.Quantity`
            The altitude of the sun below which it is considered to be "night"
            (inclusive).
        force_pressure_zero : bool (optional)
            Force the pressure to zero for solar altitude calculations. This
            avoids errors in the altitude of the Sun that can occur when the
            Sun is below the horizon and the corrections for atmospheric
            refraction return nonsense values.
        """
        self.max_solar_altitude = max_solar_altitude
        self.force_pressure_zero = force_pressure_zero

    @classmethod
    def twilight_civil(cls, **kwargs):
        """
        Consider nighttime as time between civil twilights (-6 degrees).
        """
        return cls(max_solar_altitude=-6 * u.deg, **kwargs)

    @classmethod
    def twilight_nautical(cls, **kwargs):
        """
        Consider nighttime as time between nautical twilights (-12 degrees).
        """
        return cls(max_solar_altitude=-12 * u.deg, **kwargs)

    @classmethod
    def twilight_astronomical(cls, **kwargs):
        """
        Consider nighttime as time between astronomical twilights (-18 degrees).
        """
        return cls(max_solar_altitude=-18 * u.deg, **kwargs)

    def _get_solar_altitudes(self, times, observer, targets):
        if not hasattr(observer, '_altaz_cache'):
            observer._altaz_cache = {}

        aakey = _make_cache_key(times, 'sun')

        if aakey not in observer._altaz_cache:
            try:
                if self.force_pressure_zero:
                    observer_old_pressure = observer.pressure
                    observer.pressure = 0

                # find solar altitude at these times
                altaz = observer.altaz(times, get_sun(times))
                altitude = altaz.alt
                # cache the altitude
                observer._altaz_cache[aakey] = dict(times=times,
                                                    altitude=altitude)
            finally:
                if self.force_pressure_zero:
                    observer.pressure = observer_old_pressure
        else:
            altitude = observer._altaz_cache[aakey]['altitude']

        return altitude

    # @profile
    def compute_constraint(self, times, observer, targets):
        solar_altitude = self._get_solar_altitudes(times, observer, targets)
        mask = solar_altitude <= self.max_solar_altitude

        try:
            score = np.prod(mask)
            if score > 0.0:
                # AtNightConstraint passes.  The normal case.
                # comment = 'AtNightConstraint passes'
                comment = ''
            else:
                comment = 'Constraint is zero'

            if len(comment) > 0:
                self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return mask

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#


class SunSeparationConstraint(Constraint):
    """
    Constrain the distance between the Sun and some targets.
    """

    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between Sun and target (inclusive).
            `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between Sun and target (inclusive).
            `None` indicates no limit.
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        # use get_body rather than get sun here, since
        # it returns the Sun's coordinates in an observer
        # centred frame, so the separation is as-seen
        # by the observer.
        # 'get_sun' returns ICRS coords.
        sun = get_body('sun', times, location=observer.location)
        solar_separation = sun.separation(targets)

        if self.min is None and self.max is not None:
            mask = self.max >= solar_separation
        elif self.max is None and self.min is not None:
            mask = self.min <= solar_separation
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= solar_separation) &
                    (solar_separation <= self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "SunSeparationConstraint.")
        return mask

#
# This constraint is unchanged from astroplan 0.4 code
#


class MoonSeparationConstraint(Constraint):
    """
    Constrain the distance between the Earth's moon and some targets.
    """

    def __init__(self, min=None, max=None, ephemeris=None):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between moon and target (inclusive).
            `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable separation between moon and target (inclusive).
            `None` indicates no limit.
        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            ``astropy.coordinates.solar_system_ephemeris.set`` (which is
            set to 'builtin' by default).
        """
        self.min = min
        self.max = max
        self.ephemeris = ephemeris

    # @profile
    def compute_constraint(self, times, observer, targets):
        # removed the location argument here, which causes small <1 deg
        # innacuracies, but it is needed until astropy PR #5897 is released
        # which should be astropy 1.3.2
        moon = get_moon(times,
                        ephemeris=self.ephemeris)
        # note to future editors - the order matters here
        # moon.separation(targets) is NOT the same as targets.separation(moon)
        # the former calculates the separation in the frame of the moon coord
        # which is GCRS, and that is what we want.
        moon_separation = moon.separation(targets)

        if self.min is None and self.max is not None:
            mask = self.max >= moon_separation
        elif self.max is None and self.min is not None:
            mask = self.min <= moon_separation
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= moon_separation) &
                    (moon_separation <= self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "MoonSeparationConstraint.")

        if False:
            print("MoonSeparation: ", repr(mask))

        # New code for special case logging.
        try:
            # Testing if all are true
            score = np.all(mask)
            if False:
                print("MoonSeparation score", repr(score))
            if score:
                # MoonSeparationConstraint passes.  The normal case.
                # comment = 'MoonSeparationConstraint passes'
                comment = ''
            else:
                comment = 'MoonSeparationConstraint fails in at least some cases.'

            if False and len(comment) > 0:
                if True:
                    print("MoonSeparation comment: ", comment)
                self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return mask

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#


class MoonIlluminationConstraint(Constraint):
    """
    Constrain the fractional illumination of the Earth's moon.

    Constraint is also satisfied if the Moon has set.
    """

    def __init__(self, min=None, max=None, ephemeris=None):
        """
        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            `~astropy.coordinates.solar_system_ephemeris` (which is
            set to 'builtin' by default).
        """
        self.min = min
        self.max = max
        self.ephemeris = ephemeris

    @classmethod
    def dark(cls, min=None, max=0.25, **kwargs):
        """
        initialize a `~astroplan.constraints.MoonIlluminationConstraint`
        with defaults of no minimum and a maximum of 0.25

        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        """
        return cls(min, max, **kwargs)

    @classmethod
    def grey(cls, min=0.25, max=0.65, **kwargs):
        """
        initialize a `~astroplan.constraints.MoonIlluminationConstraint`
        with defaults of a minimum of 0.25 and a maximum of 0.65

        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        """
        return cls(min, max, **kwargs)

    @classmethod
    def bright(cls, min=0.65, max=None, **kwargs):
        """
        initialize a `~astroplan.constraints.MoonIlluminationConstraint`
        with defaults of a minimum of 0.65 and no maximum

        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        """
        return cls(min, max, **kwargs)

    def compute_constraint(self, times, observer, targets):
        # first is the moon up?
        cached_moon = _get_moon_data(times, observer)
        moon_alt = cached_moon['altaz'].alt
        moon_down_mask = moon_alt < 0
        moon_up_mask = moon_alt >= 0

        illumination = cached_moon['illum']
        if self.min is None and self.max is not None:
            mask = (self.max >= illumination) | moon_down_mask
        elif self.max is None and self.min is not None:
            mask = (self.min <= illumination) & moon_up_mask
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= illumination) &
                    (illumination <= self.max)) & moon_up_mask
        else:
            raise ValueError("No max and/or min specified in "
                             "MoonSeparationConstraint.")

        return mask

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#


class LocalTimeConstraint(Constraint):
    """
    Constrain the observable hours.
    """

    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~datetime.time`
            Earliest local time (inclusive). `None` indicates no limit.

        max : `~datetime.time`
            Latest local time (inclusive). `None` indicates no limit.

        Examples
        --------
        Constrain the observations to targets that are observable between
        23:50 and 04:08 local time:

        >>> from astroplan import Observer
        >>> from astroplan.constraints import LocalTimeConstraint
        >>> import datetime as dt
        >>> subaru = Observer.at_site("Subaru", timezone="US/Hawaii")
        >>> # bound times between 23:50 and 04:08 local Hawaiian time
        >>> constraint = LocalTimeConstraint(min=dt.time(23,50), max=dt.time(4,8))
        """

        self.min = min
        self.max = max

        if self.min is None and self.max is None:
            raise ValueError("You must at least supply either a minimum or a maximum time.")

        if self.min is not None:
            if not isinstance(self.min, datetime.time):
                raise TypeError("Time limits must be specified as datetime.time objects.")

        if self.max is not None:
            if not isinstance(self.max, datetime.time):
                raise TypeError("Time limits must be specified as datetime.time objects.")

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~astroplan.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~mmtqueue.Target`
          The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """
        timezone = None

        # get timezone from time objects, or from observer
        if self.min is not None:
            timezone = self.min.tzinfo

        elif self.max is not None:
            timezone = self.max.tzinfo

        if timezone is None:
            timezone = observer.timezone

        if self.min is not None:
            min_time = self.min
        else:
            min_time = self.min = datetime.time(0, 0, 0)

        if self.max is not None:
            max_time = self.max
        else:
            max_time = datetime.time(23, 59, 59)

        # If time limits occur on same day:
        if self.min < self.max:
            try:
                mask = np.array([min_time <= t.time() <= max_time for t in times.datetime])
            except BaseException:                # use np.bool so shape queries don't cause problems
                mask = np.bool_(min_time <= times.datetime.time() <= max_time)

        # If time boundaries straddle midnight:
        else:
            try:
                mask = np.array([(t.time() >= min_time) or
                                 (t.time() <= max_time) for t in times.datetime])
            except BaseException:
                mask = np.bool_((times.datetime.time() >= min_time) or
                                (times.datetime.time() <= max_time))
        return mask


#
#  Code reviewed: 2018-01-01.  Should be correct.  Need to write a test cases.
#
# New version of TimeConstraint that handles multiple time intervals
# instead of a single tiem interval.
# The internal logic for calculating the constraint is the same as for
# the astroplan 0.4 code.


class TimeConstraint(Constraint):
    """Constrain the observing time to be within certain time limits.
    An example use case for this class would be to associate an acceptable
    time range with a specific observing block. This can be useful if not
    all observing blocks are valid over the time limits used in calls
    to `is_observable` or `is_always_observable`.
    """

    def __init__(self, min=None, max=None, constraint_times=None):
        """
        Parameters
        ----------
        min : `~astropy.time.Time`
            Earliest time (inclusive). `None` indicates no limit.
        max : `~astropy.time.Time`
            Latest time (inclusive). `None` indicates no limit.
        Examples
        --------
        Constrain the observations to targets that are observable between
        2016-03-28 and 2016-03-30:
        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> subaru = Observer.at_site("Subaru")
        >>> t1 = Time("2016-03-28T12:00:00")
        >>> t2 = Time("2016-03-30T12:00:00")
        >>> constraint = TimeConstraint(t1,t2)
        """
        self.min = min
        self.max = max
        self.constraint_times = constraint_times

        if self.constraint_times is None:
            # Original case where there is only a single min and max time
            if self.min is None and self.max is None:
                raise ValueError("You must at least supply either a minimum or a "
                                 "maximum time.")

            if self.min is not None:
                if not isinstance(self.min, Time):
                    raise TypeError("Time limits must be specified as "
                                    "astropy.time.Time objects.")

            if self.max is not None:
                if not isinstance(self.max, Time):
                    raise TypeError("Time limits must be specified as "
                                    "astropy.time.Time objects.")

        else:
            # New case where constraint_times is an array of Time duples,
            # e.g., [[TimeStart1, TimeEnd1],[TimeStart2, TimeEnd2],[ TimeStart3, TimeEnd3]]
            #
            # Check each start Time and end Time object.
            for [min_t, max_t] in self.constraint_times:
                if min_t is None and max_t is None:
                    raise ValueError("You must at least supply either a minimum or a "
                                     "maximum time.")

                if min_t is not None:
                    if not isinstance(min_t, Time):
                        raise TypeError("Time limits must be specified as "
                                        "astropy.time.Time objects.")

                if max_t is not None:
                    if not isinstance(max_t, Time):
                        raise TypeError("Time limits must be specified as "
                                        "astropy.time.Time objects.")

    # @profile
    def get_score(self, times, min_time, max_time):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """
        return np.logical_and(times > min_time, times < max_time)

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.
         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~mmtqueue.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~mmtqueue.Target`
          The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """
        if self.constraint_times is None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                min_time = Time("1950-01-01T00:00:00") if self.min is None else self.min
                max_time = Time("2120-01-01T00:00:00") if self.max is None else self.max
                # mask = np.logical_and(times > min_time, times < max_time)
                mask = self.get_score(times, min_time, max_time)

        else:
            first = True
            for [min_t, max_t] in self.constraint_times:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    min_time = Time("1950-01-01T00:00:00") if min_t is None else min_t
                    max_time = Time("2120-01-01T00:00:00") if max_t is None else max_t
                    m2 = self.get_score(times, min_time, max_time)
                    # The first time mask will be None.
                    if first:
                        first = False
                        mask = m2
                    else:
                        # If it's True for any of the possible time periods, it's True
                        # So use "logical_or".
                        mask = np.logical_or(m2, mask)

        # New code for logging comments for special cases.
        try:
            score = np.prod(mask)
            if score > 0.0:
                # TimeConstraint passes.  The normal case.
                comment = ''
            else:
                comment = 'TimeConstraint fails'

            if len(comment) > 0:
                self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return mask


# Original TimeConstraint from astroplan 0.4
# =============================================================================
#
# class TimeConstraint(Constraint):
#     """Constrain the observing time to be within certain time limits.
#
#     An example use case for this class would be to associate an acceptable
#     time range with a specific observing block. This can be useful if not
#     all observing blocks are valid over the time limits used in calls
#     to `is_observable` or `is_always_observable`.
#     """
#
#     def __init__(self, min=None, max=None):
#         """
#         Parameters
#         ----------
#         min : `~astropy.time.Time`
#             Earliest time (inclusive). `None` indicates no limit.
#
#         max : `~astropy.time.Time`
#             Latest time (inclusive). `None` indicates no limit.
#
#         Examples
#         --------
#         Constrain the observations to targets that are observable between
#         2016-03-28 and 2016-03-30:
#
#         >>> from astroplan import Observer
#         >>> from astropy.time import Time
#         >>> subaru = Observer.at_site("Subaru")
#         >>> t1 = Time("2016-03-28T12:00:00")
#         >>> t2 = Time("2016-03-30T12:00:00")
#         >>> constraint = TimeConstraint(t1,t2)
#         """
#         self.min = min
#         self.max = max
#
#         if self.min is None and self.max is None:
#             raise ValueError("You must at least supply either a minimum or a "
#                              "maximum time.")
#
#         if self.min is not None:
#             if not isinstance(self.min, Time):
#                 raise TypeError("Time limits must be specified as "
#                                 "astropy.time.Time objects.")
#
#         if self.max is not None:
#             if not isinstance(self.max, Time):
#                 raise TypeError("Time limits must be specified as "
#                                 "astropy.time.Time objects.")
#
#     def compute_constraint(self, times, observer, targets):
#         with warnings.catch_warnings():
#             warnings.simplefilter('ignore')
#             min_time = Time("1950-01-01T00:00:00") if self.min is None else self.min
#             max_time = Time("2120-01-01T00:00:00") if self.max is None else self.max
#         mask = np.logical_and(times > min_time, times < max_time)
#         return mask
#
# =============================================================================

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#
class PrimaryEclipseConstraint(Constraint):
    """
    Constrain observations to times during primary eclipse.
    """

    def __init__(self, eclipsing_system):
        """
        Parameters
        ----------
        eclipsing_system : `~astroplan.periodic.EclipsingSystem`
            System which must be in primary eclipse.
        """
        self.eclipsing_system = eclipsing_system

    def compute_constraint(self, times, observer=None, targets=None):
        mask = self.eclipsing_system.in_primary_eclipse(times)
        return mask

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#


class SecondaryEclipseConstraint(Constraint):
    """
    Constrain observations to times during secondary eclipse.
    """

    def __init__(self, eclipsing_system):
        """
        Parameters
        ----------
        eclipsing_system : `~astroplan.periodic.EclipsingSystem`
            System which must be in secondary eclipse.
        """
        self.eclipsing_system = eclipsing_system

    def compute_constraint(self, times, observer=None, targets=None):
        mask = self.eclipsing_system.in_secondary_eclipse(times)
        return mask

#
# Not used.  This constraint is unchanged from astroplan 0.4 code
#


class PhaseConstraint(Constraint):
    """
    Constrain observations to times in some range of phases for a periodic event
    (e.g.~transiting exoplanets, eclipsing binaries).
    """

    def __init__(self, periodic_event, min=None, max=None):
        """
        Parameters
        ----------
        periodic_event : `~astroplan.periodic.PeriodicEvent` or subclass
            System on which to compute the phase. For example, the system
            could be an eclipsing or non-eclipsing binary, or exoplanet system.
        min : float (optional)
            Minimum phase (inclusive) on interval [0, 1). Default is zero.
        max : float (optional)
            Maximum phase (inclusive) on interval [0, 1). Default is one.

        Examples
        --------
        To constrain observations on orbital phases between 0.4 and 0.6,
        >>> from astroplan import PeriodicEvent
        >>> from astropy.time import Time
        >>> import astropy.units as u
        >>> binary = PeriodicEvent(epoch=Time('2017-01-01 02:00'), period=1*u.day)
        >>> constraint = PhaseConstraint(binary, min=0.4, max=0.6)

        The minimum and maximum phase must be described on the interval [0, 1).
        To constrain observations on orbital phases between 0.6 and 1.2, for
        example, you should subtract one from the second number:
        >>> constraint = PhaseConstraint(binary, min=0.6, max=0.2)
        """
        self.periodic_event = periodic_event
        if (min < 0) or (min > 1) or (max < 0) or (max > 1):
            raise ValueError('The minimum of the PhaseConstraint must be within'
                             ' the interval [0, 1).')
        self.min = min if min is not None else 0.0
        self.max = max if max is not None else 1.0

    def compute_constraint(self, times, observer=None, targets=None):
        phase = self.periodic_event.phase(times)

        mask = np.where(self.max > self.min,
                        (phase >= self.min) & (phase <= self.max),
                        (phase >= self.min) | (phase <= self.max))
        return mask

#
# Unchanged from astroplan 0.4 code
#


def is_always_observable(constraints, observer, targets, times=None,
                         time_range=None, time_grid_resolution=0.5 * u.hour):
    """
    A function to determine whether ``targets`` are always observable throughout
    ``time_range`` given constraints in the ``constraints_list`` for a
    particular ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_grid_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution,
                                      grid_times_targets=True)
                           for constraint in constraints]
    constraint_arr = np.logical_and.reduce(applied_constraints)
    return np.all(constraint_arr, axis=1)

#
# Unchanged from astroplan 0.4 code
#


def is_observable(constraints, observer, targets, times=None,
                  time_range=None, time_grid_resolution=0.5 * u.hour):
    """
    Determines if the ``targets`` are observable during ``time_range`` given
    constraints in ``constraints_list`` for a particular ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_grid_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is ever observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution,
                                      grid_times_targets=True)
                           for constraint in constraints]
    constraint_arr = np.logical_and.reduce(applied_constraints)
    return np.any(constraint_arr, axis=1)

#
# Unchanged from astroplan 0.4 code
#


def is_event_observable(constraints, observer, target, times=None,
                        times_ingress_egress=None):
    """
    Determines if the ``target`` is observable at each time in ``times``, given
    constraints in ``constraints`` for a particular ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    target : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target

    times : `~astropy.time.Time` (optional)
        Array of mid-event times on which to test the constraints

    times_ingress_egress : `~astropy.time.Time` (optional)
        Array of ingress and egress times for ``N`` events, with shape
        (``N``, 2).

    Returns
    -------
    event_observable : `~numpy.ndarray`
        Array of booleans of same length as ``times`` for whether or not the
        target is ever observable at each time, given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    if times is not None:
        applied_constraints = [constraint(observer, target, times=times,
                                          grid_times_targets=True)
                               for constraint in constraints]
        constraint_arr = np.logical_and.reduce(applied_constraints)

    else:
        times_ing = times_ingress_egress[:, 0]
        times_egr = times_ingress_egress[:, 1]
        applied_constraints_ing = [constraint(observer, target, times=times_ing,
                                              grid_times_targets=True)
                                   for constraint in constraints]
        applied_constraints_egr = [constraint(observer, target, times=times_egr,
                                              grid_times_targets=True)
                                   for constraint in constraints]

        constraint_arr = np.logical_and(np.logical_and.reduce(applied_constraints_ing),
                                        np.logical_and.reduce(applied_constraints_egr))
    return constraint_arr

#
# Unchanged from astroplan 0.4 code
#


def months_observable(constraints, observer, targets,
                      time_grid_resolution=0.5 * u.hour):
    """
    Determines which month the specified ``targets`` are observable for a
    specific ``observer``, given the supplied ``constriants``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    time_grid_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    observable_months : list
        List of sets of unique integers representing each month that a target is
        observable, one set per target. These integers are 1-based so that
        January maps to 1, February maps to 2, etc.

    """
    # TODO: This method could be sped up a lot by dropping to the trigonometric
    # altitude calculations.
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    # Calculate throughout the year of 2014 so as not to require forward
    # extrapolation off of the IERS tables
    time_range = Time(['2014-01-01', '2014-12-31'])
    times = time_grid_from_range(time_range, time_grid_resolution)

    # TODO: This method could be sped up a lot by dropping to the trigonometric
    # altitude calculations.

    applied_constraints = [constraint(observer, targets,
                                      times=times,
                                      grid_times_targets=True)
                           for constraint in constraints]
    constraint_arr = np.logical_and.reduce(applied_constraints)

    months_observable = []
    for target, observable in zip(targets, constraint_arr):
        s = set([t.datetime.month for t in times[observable]])
        months_observable.append(s)

    return months_observable

#
# Unchanged from astroplan 0.4 code
#


def observability_table(constraints, observer, targets, times=None,
                        time_range=None, time_grid_resolution=0.5 * u.hour):
    """
    Creates a table with information about observability for all  the ``targets``
    over the requested ``time_range``, given the constraints in
    ``constraints_list`` for ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_grid_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    observability_table : `~astropy.table.Table`
        A Table containing the observability information for each of the
        ``targets``. The table contains four columns with information about the
        target and it's observability: ``'target name'``, ``'ever observable'``,
        ``'always observable'``, and ``'fraction of time observable'``.  It also
        contains metadata entries ``'times'`` (with an array of all the times),
        ``'observer'`` (the `~astroplan.Observer` object), and ``'constraints'``
        (containing the supplied ``constraints``).
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution,
                                      grid_times_targets=True)
                           for constraint in constraints]
    constraint_arr = np.logical_and.reduce(applied_constraints)

    colnames = ['target name', 'ever observable', 'always observable',
                'fraction of time observable']

    target_names = [target.name for target in targets]
    ever_obs = np.any(constraint_arr, axis=1)
    always_obs = np.all(constraint_arr, axis=1)
    frac_obs = np.sum(constraint_arr, axis=1) / constraint_arr.shape[1]

    tab = table.Table(names=colnames, data=[target_names, ever_obs, always_obs,
                                            frac_obs])

    if times is None and time_range is not None:
        times = time_grid_from_range(time_range,
                                     time_resolution=time_grid_resolution)

    tab.meta['times'] = times.datetime
    tab.meta['observer'] = observer
    tab.meta['constraints'] = constraints

    return tab

#
# Unchanged from astroplan 0.4 code
#


def min_best_rescale(vals, min_val, max_val, less_than_min=1):
    """
    rescales an input array ``vals`` to be a score (between zero and one),
    where the ``min_val`` goes to one, and the ``max_val`` goes to zero.

    Parameters
    ----------
    vals : array-like
        the values that need to be rescaled to be between 0 and 1
    min_val : float
        worst acceptable value (rescales to 0)
    max_val : float
        best value cared about (rescales to 1)
    less_than_min : 0 or 1
        what is returned for ``vals`` below ``min_val``. (in some cases
        anything less than ``min_val`` should also return one,
        in some cases it should return zero)

    Returns
    -------
    array of floats between 0 and 1 inclusive rescaled so that
    ``vals`` equal to ``max_val`` equal 0 and those equal to
    ``min_val`` equal 1

    Examples
    --------
    rescale airmasses to between 0 and 1, with the best (1)
    and worst (2.25). All values outside the range should
    return 0.
    >>> from astroplan.constraints import min_best_rescale
    >>> import numpy as np
    >>> airmasses = np.array([1, 1.5, 2, 3, 0])
    >>> min_best_rescale(airmasses, 1, 2.25, less_than_min = 0)
    array([ 1. ,  0.6,  0.2,  0. , 0. ])
    """
    rescaled = (vals - max_val) / (min_val - max_val)
    below = vals < min_val
    above = vals > max_val
    rescaled[below] = less_than_min
    rescaled[above] = 0

    return rescaled

#
# Unchanged from astroplan 0.4 code
#


def max_best_rescale(vals, min_val, max_val, greater_than_max=1):
    """
    rescales an input array ``vals`` to be a score (between zero and one),
    where the ``max_val`` goes to one, and the ``min_val`` goes to zero.

    Parameters
    ----------
    vals : array-like
        the values that need to be rescaled to be between 0 and 1
    min_val : float
        worst acceptable value (rescales to 0)
    max_val : float
        best value cared about (rescales to 1)
    greater_than_max : 0 or 1
        what is returned for ``vals`` above ``max_val``. (in some cases
        anything higher than ``max_val`` should also return one,
        in some cases it should return zero)

    Returns
    -------
    array of floats between 0 and 1 inclusive rescaled so that
    ``vals`` equal to ``min_val`` equal 0 and those equal to
    ``max_val`` equal 1

    Examples
    --------
    rescale an array of altitudes to be between 0 and 1,
    with the best (60) going to 1 and worst (35) going to
    0. For values outside the range, the rescale should
    return 0 below 35 and 1 above 60.
    >>> from astroplan.constraints import max_best_rescale
    >>> import numpy as np
    >>> altitudes = np.array([20, 30, 40, 45, 55, 70])
    >>> max_best_rescale(altitudes, 35, 60)
    array([ 0. , 0. , 0.2, 0.4, 0.8, 1. ])
    """
    rescaled = (vals - min_val) / (max_val - min_val)
    below = vals < min_val
    above = vals > max_val
    rescaled[below] = 0
    rescaled[above] = greater_than_max

    return rescaled


#
#  Code reviewed 2018-01-01.  Looks OK. Need to write test cases.
#
class MaskAngleConstraint(Constraint):
    """
       MaskAngleConstraint.

    """

    def __init__(self, observer,
                 design_parang=0.0 * u.deg,
                 max_mask_angle=30.0 * u.deg,
                 grid_times_targets=False,
                 debug=False):
        """
        A subclass of the Constraint class that evaluates whether the target's current parallactic
        angle is within the allowed tolerance from the mask's design angle.

        Rather than range from 1.0 to 0.0, the lower limit is a small number, e.g., 0.1.  This allows
        targets to still be scheduled outside of the allowed limit, but with a low constraint score.
        If the target is outside of the allowed angle limit, little light will pass through the mask.

        Parameters
        ----------
        observer : Observer
            An instance of the Observer class.
        design_parang : astropy.quantity.deg
            The design (position) angle for the mask
        max_mask_angle : astropy.quantity.deg
            The maximum allowed deviation of the target's parallactic angle from the design angle.
        grid_times_targets : bool
            A boolean from the Constraint class on whether to grid the constraint over both times and targets
        debug : bool
            A boolean on whether to print debug statements.

        Returns
        -------
        self : MaskAngleConstraint
            A subclass of Constraint
        """
        if False:
            print("MaskAngleConstraint init")
        try:
            design_parang = Angle(self.design_parang)
        except:
            design_parang = Angle(0 * u.deg)

        try:
            max_mask_angle = Angle(self.max_mask_angle)
        except:
            max_mask_angle = Angle(30 * u.deg)
        self.observer = observer
        self.design_parang = design_parang
        self.max_mask_angle = max_mask_angle
        self.grid_times_targets = grid_times_targets
        self.debug = debug

    # @profile
    def get_score(self, target, time):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """
        if self.debug:
            print("target: ", repr(target), ", time: ", repr(time))

        pa = self.observer.parallactic_angle(time, target)
        da = self.design_parang
        ma = self.max_mask_angle
        if self.debug:
            print("pa: ", repr(pa))
            print("dp: ", repr(da))
            print("ma: ", repr(ma))
        ang = abs(pa - da) <= abs(ma)
        if ang == True:
            score = 1.0
        else:
            # score = 0.1
            score = 0.5

        # New code for logging comments for special caseVs.
        try:
            if score < 1.0:
                comment = 'MaskAngleConstraint fails with angle (deg): {}'.format(round((pa - da).degree, 2))
            else:
                # comment = 'MaskAngleConstraint passes'
                comment = ''

            if len(comment) > 0:
                self.process_comment(target, time, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return score

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
        Actually do the real work of computing the constraint.

        Parameters
        ----------
        times : `~astropy.time.Time`
            The times to compute the constraint
        observer : `~astroplan.Observer`
            The observaton location from which to apply the constraints
        targets : sequence of `~astroplan.Target`
            The targets on which to apply the constraints.

        Returns
        -------
        mask_numpy : 1D or 2D Numpy array of float or bool
            The constraints, with targets along the first index and times along
            the second.
        """
        if self.debug:
            print("In MaskAngleContraint compute_constraint")
        mask = []
        if targets.isscalar:
            target = targets
            for time in times:
                score = self.get_score(target, time)
                if self.debug:
                    print("Scalar MaskAngleConstraint score: ", repr(score))
                mask.append(score)
            if self.debug:
                print("Handling as scalar")
            # Turn the mask into a numpy array
            mask_numpy = np.array(mask)
            # I don't see anyway to implement the grid_times_targets here.
            # If the targets are scalar, the result will just a 1-D array.
            # It will broadcast as needed by numpy for any further calculations.
        else:
            for target in targets:
                for time in times:
                    score = self.get_score(target, time)
                    if self.debug:
                        print("Iterater MaskAngleConstraint score: ", repr(score))
                    mask.append(score)
            if self.debug:
                print("Handling as iterable")
            # Turn the mask into a numpy array and reshape.
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if self.debug:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask")
            print(repr(mask))
            print("mask_numpy")
            print(repr(mask_numpy))

        return mask_numpy

#
#  Code reviewed 2018-01-01.  Looks OK. Need to write test cases.
#


class MaskNumberConstraint(Constraint):
    """
       MaskNumberConstraint.  Determines if the maximum number of masks used as reached
       the allowed limit.  Returns True/1 if the masks_used is less than or equal to the
       number of allowed masks.  Returns False/0 otherwise.

    """

    def __init__(self, mask_id=None,
                 masks_used=None,
                 grid_times_targets=False,
                 debug=False):
        self.mask_id = mask_id
        self.masks_used = masks_used
        self.grid_times_targets = grid_times_targets
        self.debug = debug

    # @profile
    def get_score(self, time, target):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """
        if len(self.masks_used) < 10 or \
                self.mask_id in self.masks_used:
            # Just a note that we don't want to include this mask_id in
            # masks_used here since the block may or may not ultimately
            # get scheduled.  "masks_used" will be updated in the scheduler
            # with updated value passed to this class.
            score = 1.0
        else:
            score = 0.0

        # New code for logging comments for special caseVs.
        try:
            if score < 1.0:
                comment = 'MaskNumberConstraint fails'
            else:
                # comment = 'MaskNumberConstraint passes'
                comment = ''
            if len(comment) > 0:
                self.process_comment(target, time, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return score

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
             The times to compute the constraint
         observer : `~mmtqueue.Observer`
             The observaton location from which to apply the constraints
         targets : sequence of `~mmtqueue.Target`
             The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
             The constraints, with targets along the first index and times along
             the second.
         """
        mask = []
        if targets.isscalar:
            target = targets
            for time in times:
                score = self.get_score(target, time)
                mask.append(score)
            mask_numpy = np.array(mask)
        else:
            for target in targets:
                for time in times:
                    score = self.get_score(target, time)
                    mask.append(score)
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if False:
            print("targets", repr(targets))
            print("times", repr(times))
            print("mask_id", repr(self.mask_id))
            print("mask", repr(mask))

        return mask_numpy


#
#
#
class MaskReadyConstraint(Constraint):
    """
    MaskReadyConstraint.  Determines of the mask is ready (i.e., available) to be used.
    If the mask has not been cut yet, it is not available and observing blocks/targets
    that require that mask cannot be scheduled.

    Parameters
    ----------
    mask_ready : bool
        Whether the mask has been cut and is available for installation on the telescope.
    debug : bool
        Boolean flag on whether to produce verbose output.
    """

    def __init__(self, mask_ready=True,
                 debug=False):
        self.mask_ready = mask_ready
        self.debug = debug

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
             The times to compute the constraint
         observer : `~astroplan.Observer`
             The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
             The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
             The constraints, with targets along the first index and times along
             the second.
         """
        # Testing if targets is scalar
        mask = []
        if self.mask_ready:
            score = True
        else:
            score = False
        if targets.isscalar:
            for t in times:
                mask.append(score)
            mask_numpy = np.array(mask)
        else:
            for target in targets:
                for t in times:
                    mask.append(score)
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if self.debug:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask_numpy")
            print(repr(mask_numpy))

        # New code for logging comments for special caseVs.
        try:
            score = np.prod(mask_numpy)
            if score < 1.0:
                comment = 'MaskReadyConstraint fails'
            else:
                # comment = 'MaskReadyConstraint passes'
                comment = ''
            if len(comment) > 0:
                self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return mask_numpy

#
#
#

class ObjectTypeConstraint(Constraint):
    """
    ObjectTypeConstraint.  Assigns a score to an observing block based upon
    the object type.  The scheduling approach results in some object types 
    being easier to schedule that others, e.g., imaging is easier than slit-mask.
    This constraints lowers the overall score for imaging objects so that other object
    types are preferentially scheduled when they can.

    Parameters
    ----------
    object_type : 'imaging', 'longslit', 'mask'
    debug : bool
        Boolean flag on whether to produce verbose output.
    """

    def __init__(self, object_type='imaging',
                 debug=False):
        self.object_type = object_type.lower()
        self.debug = debug

        # Should have these in the configuration structure.
        self.imaging_score = 0.8
        self.longslit_score = 0.9
        self.mask_score = 1.0

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
             The times to compute the constraint
         observer : `~astroplan.Observer`
             The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
             The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
             The constraints, with targets along the first index and times along
             the second.
         """
        # Testing if targets is scalar
        if targets.isscalar:
            if self.object_type == 'mask':
                mask = [self.mask_score * len(times)]
            elif self.object_type == 'longslit':
                mask = [self.longslit_score * len(times)]
            elif self.object_type == 'imaging':
                mask = [self.imaging_score * len(times)]
            else:
                mask = [1.0 * len(times)]
            mask_numpy = np.array(mask)
        else:
            if self.object_type == 'mask':
                mask = [self.mask_score * len(targets)]
            elif self.object_type == 'longslit':
                mask = [self.longslit_score * len(targets)]
            elif self.object_type == 'imaging':
                mask = [self.imaging_score * len(targets)]
            else:
                mask = [1.0 * len(targets)]
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if False:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask_numpy")
            print(repr(mask_numpy))

        # New code for logging comments for special caseVs.
        # comment = ''
        try:
            score = np.prod(mask_numpy)
            # if score < 1.0:
            #    comment = 'ObjectTypeConstraint fails'
            # else:
            #     comment = 'ObjectTyperConstraint passes'
            # if len(comment) > 0:
            #     self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return mask_numpy



"""
Not using this.  Using MaskAngleConstraint instead.  
We may want it at some point in the future.

"""


class TransitConstraint(AltitudeConstraint):
    """
    Constrain the transit behavior of a target.

    Parameters
    ----------
    mode : 1,2, or 3.
        `1`  strictly rising
        `2`  strictly setting
        `3`  no transit (i.e., either strictly rising or setting). 
    debug : bool
        `True`  verbose output
        `False` suppress output

    Examples
    --------
        Constrain to strictly rising targets with debugging output
            TransitConstraint(mode=1,debug=True)
        Constrain to strictly setting targets
            TransitConstraint(mode=2)
        Constrain to strictly rising targets without grid_times_targets with debugging output
            TransitConstraint(mode=1,grid_times_targets=False, debug=True)

    """

    def __init__(self, mode=0, debug=False):
        self.mode = mode
        self.debug = debug

    # Checking for monotonicity
    #      "having the property either of never increasing or of never decreasing as the
    #       values of the independent variable or the subscripts of the terms increase"
    # From https://stackoverflow.com/questions/4983258/python-how-to-check-list-monotonicity

    def strictly_increasing(self, x):
        # Check if altitudes are sorted in strictly_increasing order.
        sorted_arr = np.sort(x)
        test = x == sorted_arr
        if self.debug:
            print("test for strictly_increasing")
            print("x:")
            print(x)
            print("sorted_arr:")
            print(sorted_arr)
            print("test:")
            print(test)
        return test

    def strictly_decreasing(self, x):
        # Check if altitudes are sorted in strictly_increasing order,
        # then invert the output.
        test = np.invert(self.strictly_increasing(x))
        if self.debug:
            print("test for strictly_decreasing")
            print(test)
        return test

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
             The times to compute the constraint
         observer : `~astroplan.Observer`
             The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
             The targets on which to apply the constraints.

         Returns
         -------
         mask : 1D or 2D Numpy array of boolean values.
             The constraints, with targets along the first index and times along
             the second.
         """
        cached_altaz = _get_altaz(times, observer, targets)

        alt = cached_altaz['altaz'].alt
        if self.debug:
            print("alt:")
            print(alt)

        # Case for rising targets
        if self.mode == 1:
            mask = self.strictly_increasing(alt)

        # Case for setting targets
        elif self.mode == 2:
            mask = self.strictly_decreasing(alt)

        # Case of no transit, i.e., either strictly rising or setting setting targets
        elif self.mode == 3:
            mask = self.strictly_ascending(alt) or self.strictly_decreasing(alt)

        else:
            raise ValueError("Improper altitude values in MaskConstraint.")

        if self.debug:
            print("Times:")
            print(times)
            print("Targets:")
            print(targets)
            print("Mode:")
            print(self.mode)
            print("Mask:")
            print(mask)
        return mask

#
#  Code checked against QueueSchedueler1.0/mmtscheduler.py  2018-01-31 17:20
#
# Define an astroplan constraint for the distance of the target from the meridian.
# The returned value is either a boolean [0,1] if the target is outside of an allowed time
# from meridian transit or a float from [0.0:1.0], where the value is 1.0
# when the target is o/typen the meridian to 0.0 when it is at the anti-meridian (12 hours from the meridian)


class MeridianConstraint(Constraint):
    """Constrains the time for targets from meridian transit.

    Principal investigators (PI's) are required to assigned an integer priority from 1 (highest) to 3 (lowest) to each of their targets.
    The targets should be equally divided into the three priorities (i.e., 1, 2, and 3) so that 1/3 of the requested time correspondes to each
    priority.
    This equal division into the three priorities by even time requested is needed to keep scheduling fair for all projects.
    Every effort will be made to observe all targets, but PI's should anticipate that at least part of their priority 3 targets will not be observed because of poor weather or other causes.


    """

    def __init__(self,
                 duration,
                 pi_priority,
                 daily_events,
                 mode="sunset",
                 min_alt_degrees=20 * u.deg,
                 max_solar_altitude=-12 * u.deg,
                 grid_times_targets=False,
                 debug=False,
                 verbose=True):
        """
        Parameters
        ----------
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable separation (in decimal hours) between meridian and target (inclusive).
            `None` indicates no constraint of how far the target can be from the meridian.
        boolean_constraint : bool
            If True, the constraint is treated as a boolean (True for within the
            limits and False for outside).  If False, the constraint returns a
            float on [0, 1], where 0 is when the target is on the anti-meridian and
            1 is when the target is on the meridian.
        """
        self.mode = mode
        self.duration = duration
        self.pi_priority = pi_priority
        self.daily_events = daily_events
        self.min_alt_degrees = min_alt_degrees
        self.max_solar_altitude = max_solar_altitude
        self.grid_times_targets = grid_times_targets
        self.debug = debug
        self.verbose = verbose

        try:
            # pi_priority must be an integer.  Nominally, 1, 2 or 3.
            self.pi_priority = int(self.pi_priority)
        except:
            self.pi_priority = 1

        # 12 hours: the maximum possible time for a target to be from the meridian
        self.seconds_in_12hrs = 43200     # 12 hours ==> 12 * 60 * 60 = 43200 seconds

        # Set up to TimeDelta constatnt values for future use.
        self.dt_0hrs = TimeDelta(0, format='sec')
        self.dt_1hrs = TimeDelta(3600, format='sec')
        self.dt_1_5hrs = TimeDelta(3600 * 1.5, format='sec')
        self.dt_2hrs = TimeDelta(3600 * 2.0, format='sec')
        self.dt_2_5hrs = TimeDelta(3600 * 2.5, format='sec')
        self.dt_3hrs = TimeDelta(3600 * 3.0, format='sec')
        self.dt_3_5hrs = TimeDelta(3600 * 3.5, format='sec')
        self.dt_4hrs = TimeDelta(3600 * 4.0, format='sec')

        if False:
            print("MeridianConstraint initialized")

    # @profile
    def get_score(self, target, time, observer):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """

        comment = ''
        # We do a series of tests to see if the observing block is
        # is setting early in the evening.  We want to give high
        # priority to PIPriority == 1 observing blocks that are
        # setting near sunset and that can still be observed.
        #
        # We take into account the duration of observing blocks
        # when doing this special "sunset" mode.
        #
        # Step 1:  Get the time of the previous sunset.
        #       This will be used to see if the time is close to sunset.
        if use_daily_events:
            comment = 'previous_sunset'
            key = get_sun_moon_date_key(comment, time)

            if key in self.daily_events:
                prev_sun_set_time = self.daily_events[key]
                if False:
                    print("prev_sun_set_time: ", repr(prev_sun_set_time))
            else:
                prev_sun_set_time = observer.sun_set_time(time,
                                                          which='previous',
                                                          horizon=self.max_solar_altitude)
                self.daily_events[key] = prev_sun_set_time
        else:
            prev_sun_set_time = observer.sun_set_time(time,
                                                      which='previous',
                                                      horizon=self.max_solar_altitude)

        # Step 2: Calculate the time from the previous sunset.
        #        This is an indication of how close we are to sunset.
        #        It will be a small number if we are trying to observe
        #        just after sunset.
        td1 = time - prev_sun_set_time

        # Step 3: Get the time of the next target rise.
        #       This will be used to see if the target is close to rising
        #       in the east in the morning.
        if use_daily_events:
            # def get_target_date_key(target, time, mode='set'):
            comment = 'next_target_set'
            key = get_target_date_key(target, time, comment)
            if key in self.daily_events:
                next_target_set = self.daily_events[key]
                if False:
                    print("next_target_set: ", repr(next_target_set))
            else:
                next_target_set = observer.target_set_time(time, target,
                                                           which="next",
                                                           horizon=self.min_alt_degrees)
                self.daily_events[key] = next_target_set
        else:
            next_target_set = observer.target_set_time(time, target,
                                                       which="next",
                                                       horizon=self.min_alt_degrees)

        # Step 4: Calculate the time from the next target setting.
        #        This number will be a small positive number when
        #        the target is above the western horizon.
        td2 = next_target_set - time

        # if False and verbose:
        #    print("tx1: {}, sec: {}, tx2: {}, sec: {}".format(tx1, tx1.sec, tx2, tx2.sec))

        # Note: The next three steps use a "sunset" mode where we want to give
        #       high priority to targets that will be setting within the next 2-4
        #       hours.  We use a graded approach for scoring, based on the duration
        #       of the target/observing block.

        # Step 5: This is the "2-hour-target-duration" sunset special case.
        #       Evaluate for targets/observing blocks that are more the
        #       _2_ hours in duration, and we are within _4_ hours after sunset,
        #       and the target will be setting within _4_ hours.
        #       Only do this for priority 1 targets.
        #       This is our only chance to observe them.
        #       The score is set to 1.0 to give it the maximum chance of being observed.

        if self.mode == "sunset" and \
                td1 <= self.dt_4hrs and td1 > self.dt_0hrs and \
                td2 <= self.dt_4hrs and td2 > self.dt_0hrs and \
                self.duration >= 2.0 * u.hour and \
                self.pi_priority == 1.0:
            score = 1.0
            comment = "Sunset special case (>= 2-hr duration)"
            if self.verbose:
                print("Sunset special case (>= 2-hr duration), score:", score)

        # Step 6: This is the "1-hour-target-duration" sunset special case.
        #       Evaluate for targets/observing blocks that are 1-2 hours
        #       in duration, and we are within _3_ hours after sunset,
        #       and the target will be setting within _3_ hours.
        #       Only do this for priority 1 targets.
        #       This is our only chance to observe them.
        #       The score is set to 1.0 to give it the maximum chance of being observed.
        elif self.mode == "sunset" and \
                td1 <= self.dt_3hrs and td1 > self.dt_0hrs and  \
                td2 <= self.dt_3hrs and td2 > self.dt_0hrs and  \
                self.duration >= 1.0 * u.hour and \
                self.pi_priority == 1.0:
            score = 1.0
            comment = "Sunset special case, (>= 1-hour and < 2-hour duration) score:"
            if self.verbose:
                print("Sunset special case, (>= 1-hour and < 2-hour duration) score:", score)

        # Step 7: This is the "<1-hour-target-duration" sunset special case.
        #       Evaluate for targets/observing blocks that are <1 hour
        #       in duration, and we are within _2_ hours after sunset,
        #       and the target will be setting within _2_ hours.
        #       Only do this for priority 1 targets.
        #       This is our only chance to observe them.
        #       The score is set to 1.0 to give it the maximum chance of being observed.
        elif self.mode == "sunset" and  \
                td1 <= self.dt_2hrs and td1 > self.dt_0hrs and \
                td2 <= self.dt_2hrs and td2 > self.dt_0hrs and \
                self.pi_priority == 1.0:
            score = 1.0
            comment = "Sunset special case, (any duration) score:"
            if self.verbose:
                print("Sunset special case, (any duration) score:", score)

        # Step 8:  If all of the other conditions have not been true,
        #       Determine how far the target is from the meridian in seconds
        #       and divide by 12 hours (== 43200 seconds)
        #       The target can be in either rising towards the meridian or
        #       setting away from the meridian.
        else:
            if use_daily_events:
                # def get_target_date_key(target, time, mode='set'):
                comment = 'nearest_target_meridian'
                key = get_target_date_key(target, time, comment)
                if key in self.daily_events:
                    meridian_time = self.daily_events[key]
                    if False:
                        print("meridian_time: ", repr(meridian_time))
                else:
                    meridian_time = observer.target_meridian_transit_time(time,
                                                                          target,
                                                                          which='nearest')
                    self.daily_events[key] = meridian_time
            else:
                meridian_time = observer.target_meridian_transit_time(time,
                                                                      target,
                                                                      which='nearest')
            diff = abs(time.unix - meridian_time.unix)

            # There are times when the meridian time is 24 hours off.
            # So, the math here accounts for that.
            # If the time difference is more than 24 hours (43200 seconds),
            # subtract 24 hours.
            while diff > self.seconds_in_12hrs:
                diff -= self.seconds_in_12hrs * 2  # 24 hours in seconds.
                # Recheck that we are using an absolute value.
                diff = abs(diff)

            # Here is the meridian scoring algorithm.
            # The closer to the meridian, the closer the score is to 1.0.
            # The range of scores is 1.0 (on the meridian) to
            # 0.0 (on the anti-meridian).
            score = 1.0 - (diff / self.seconds_in_12hrs)

            # The score should already range from 0.0 to 1.0.
            if score < 0.0:
                score = 0.0
            if score > 1.0:
                score = 1.0

            # New code for logging comments for special cases.
            try:
                if len(comment) > 0:
                    self.process_comment(target, time, score, comment)
            except:
                print("traceback: ", traceback.format_exc())

        if False:
            print("score: ", repr(score))
        return score

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
             The times to compute the constraint
         observer : `~astroplan.Observer`
             The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
             The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
             The constraints, with targets along the first index and times along
             the second.
         """

        """
        The MeridianConstraint is calculated by: 1) determining the number of hours the target is from the meridian, and 2) calculating a constraint using Math.abs((12. - hours_from_meridan)/12.0) for the target's position at the beginning, middle, and end of the observing block.  
        The calculated scores for these three times will be different.  
        This causes the constraint to equal 1.0 on the meridian and 0.0 on the anti-meridian (12 hours away).  
        Since the absolute value is used, it doesn't matter which direction the target is from the meridan, i.e., positive hours or negative hours.  Values will always vary from 1.0 to 0.0
        It is possible that the target passes through the meridian during the observing block, i.e., it "transits". 
        Caution should be used in cases where the target transits in that azimuth velocities can be very large if the target is close to zenith.
        The maximum AltitudeConstraint should help prevent extremely large azimuth velocities.  

        It should be remembered that constraint scores are calculated at the beginning, middle, and end of each observing block as part of score for the block.
        This causes the constraint to be multiplied by itself three times and the constraint to vary as 1/X^^3 rather than 1/X.

        """
        if False:
            print("MeridianConstraint computing")

        mask = []
        # Testing if targets is scalar
        if targets.isscalar:
            target = targets
            for time in times:
                score = self.get_score(target, time, observer)
                mask.append(score)
            mask_numpy = np.array(mask)
        else:
            for target in targets:
                for time in times:
                    score = self.get_score(target, time, observer)
                    mask.append(score)
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        return mask_numpy




#
#  Code checked against QueueScheduler1.0/mmtscheduler.py  2018-01-31 17:20
#
class PIPriorityConstraint(Constraint):
    """
       PIPriorityConstraint.

    """

    def __init__(self, pi_priority=1.0,
                 grid_times_targets=False,
                 debug=False):
        try:
            self.pi_priority = float(pi_priority)
        except:
            self.pi_priority = 1.0
        self.grid_times_targets = grid_times_targets
        self.debug = debug

    # @profile
    def compute_constraint(self, times, observer, targets):
        """
        Actually do the real work of computing the constraint.

        Parameters
        ----------
        times : `~astropy.time.Time`
         The times to compute the constraint
        observer : `~astroplan.Observer`
         The observaton location from which to apply the constraints
        targets : sequence of `~astroplan.Target`
         The targets on which to apply the constraints.

        Returns
        -------
        mask_numpy : 1D or 2D Numpy array of float or bool
        The constraints, with targets along the first index and times along
        the second.
        """
        # Testing if targets is scalar
        if targets.isscalar:
            # Shouldn't have a negative pi_priority, but someone may try.
            mask = [1.0 / abs(float(self.pi_priority))
                    for time in times]
            mask_numpy = np.array(mask)
        else:
            mask = [([1.0 / abs(float(self.pi_priority))
                      for time in times])
                    for target in targets]
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if False:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask_numpy")
            print(repr(mask_numpy))

        return mask_numpy


#
# This is not working as of 2018-06-25
# I think the problem is the found_ToO and which_ToO
# are not getting updated globally and the new
# values passed down to this constraint.
#
class TargetOfOpportunityConstraint(Constraint):
    """
       TargetOfOpportunityConstraint.

    """

    def __init__(self, ToO=None,
                 grid_times_targets=False,
                 debug=False):
        self.ToO = ToO
        self.grid_times_targets = grid_times_targets
        self.debug = debug


    def get_score(self, target):

        if self.ToO:
            found_ToO = self.ToO['found_ToO']
            observable_ToO = self.ToO['observable_ToO']
        else:
            found_ToO = False
            observable_ToO = []

        if False:
            print("954 ToO: ", repr(self.ToO))
        if found_ToO and len(observable_ToO) > 0:
            # Note observable_ToO of SkyCoord's, not FixedTarget's.
            if target.coord in observable_ToO:
                score = 1.0
            else:
                score = 0.1
        else:
            score = 1.0

        if False:
            print("954 ToO score: ", score)

        return score



    # @profile
    def compute_constraint(self, times, observer, targets):
        """
        Actually do the real work of computing the constraint.

        Parameters
        ----------
        times : `~astropy.time.Time`
         The times to compute the constraint
        observer : `~astroplan.Observer`
         The observaton location from which to apply the constraints
        targets : sequence of `~astroplan.Target`
         The targets on which to apply the constraints.

        Returns
        -------
        mask_numpy : 1D or 2D Numpy array of float or bool
        The constraints, with targets along the first index and times along
        the second.
        """
        if False:
            print ("153 ToO: ", self.ToO)

        mask = []

        try:
            if targets.isscalar:
                target = targets
                score = self.get_score(target)
                for t in times:
                    mask.append(score)
                mask_numpy = np.array(mask)
            else:
                for target in targets:
                    score = self.get_score(target)
                    for t in times:
                        mask.append(score)
                mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        except:
            print("153, traceback: ", traceback.format_exc())


        if False:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask_numpy")
            print(repr(mask_numpy))

        return mask_numpy


#
#  Code checked against QueueSchedueler1.0/mmtscheduler.py  2018-01-31 17:20
#
class RotatorConstraint(Constraint):
    # Passing in

    def __init__(self, max=None, min=None,
                 posang=None,
                 instrument_offset=None,
                 grid_times_targets=False,
                 debug=True):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable rotator angle (inclusive).
            `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable rotator angle (inclusive).
            `None` indicates no limit.
        """
        if False:
            print("RotatorConstraint max", repr(max))
            print("RotatorConstraint min", repr(min))
            print("RotatorConstraint posang", repr(posang))

        if max == None:
            self.max = 180.0 * u.deg
        else:
            self.max = max

        if min == None:
            self.min = -180.0 * u.deg
        else:
            self.min = min

        if posang == None:
            self.posang = 0.0 * u.deg
        # It can't seem to handle 0.0 degrees???
        elif abs(float(posang.value)) < 0.1:
            self.posang = 0.1 * u.deg
        else:
            self.posang = posang

        if instrument_offset == None:
            self.instrument_offset = 0.0 * u.deg
        else:
            self.instrument_offset = instrument_offset * u.deg

        self.target_posang = Angle(self.posang)
        self.instrument_offset_angle = Angle(self.instrument_offset)
        self.upper_limit = Angle(self.max)
        self.lower_limit = Angle(self.min)

        self.grid_times_targets = grid_times_targets
        self.debug = debug

    def get_score(self, target, time, observer):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """
        parang = observer.parallactic_angle(time, target)
        # rotator angle = parallactic_angle - position_angle
        # "sky offset" == "position angle"
        # rotator position == rototor angle + rotator offset
        #
        details = False

        if details:
            print("parang: ", repr(parang.degree))
            print("target_posang: ", repr(self.target_posang.degree))

        rot_ang = parang - self.target_posang - self.instrument_offset_angle

        if details:
            print("rot_ang1: ", repr(rot_ang.degree))
        if rot_ang < -180.0 * u.deg:
            rot_ang = rot_ang + 360.0 * u.deg
        if details:
            print("rot_ang2: ", repr(rot_ang.degree))

        upper_test = rot_ang <= self.upper_limit
        if details:
            print("upper_test1", repr(upper_test))
        # New code to test limits again at 360 degrees rotation.  JDG 2016-12-05
        if not upper_test:
            x = rot_ang - 360 * u.deg
            # Make sure that we are still within the rotator limits after
            # subtracting 360 degrees
            if self.lower_limit <= x and x <= self.upper_limit:
                upper_test = x <= self.upper_limit
        if details:
            print("upper_test2", repr(upper_test))

        lower_test = self.lower_limit <= rot_ang
        if details:
            print("lower_test1", repr(lower_test))
        # New code to test limits again at 360 degrees rotation.  JDG 2016-12-05
        if not lower_test:
            x = rot_ang + 360.0 * u.deg
            # Make sure that we are still within the rotator limits after
            # adding 360 degrees
            if self.lower_limit <= x and x <= self.upper_limit:
                lower_test = self.lower_limit <= x
        if details:
            print("lower_test2", repr(lower_test))

        score = lower_test & upper_test

        # New code for logging comments for special cases.
        try:
            score = np.prod(score)

            if score == 1.0:
                # RotatorConstraint passes
                comment = ''
            else:
                comment = "RotatorConstraint fails, rotator angle (deg): {}".format(round(rot_ang.degree, 2))

            if len(comment) > 0:
                self.process_comment(target, time, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        if details:
            print("RotatorConstraint score", repr(score))
        return score

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~astroplan.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
          The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """
        mask = []
        # Testing if targets is scalar
        if targets.isscalar:
            target = targets
            for time in times:
                score = self.get_score(target, time, observer)
                mask.append(score)
            mask_numpy = np.array(mask)
        else:
            for target in targets:
                for time in times:
                    score = self.get_score(target, time, observer)
                    mask.append(score)
            mask_numpy = np.reshape(np.array(mask), [len(targets), len(times)])

        if False:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask_numpy")
            print(repr(mask_numpy))

        return mask_numpy


#
#  Code reviewed 2018-01-01.  Looks OK. Need to write test cases.
#
class TimeAllocationConstraint(Constraint):
    """
       TimeAllocationConstraint.

    """

    def __init__(self, program,
                 stats,
                 grid_times_targets=False,
                 debug=False):
        self.program = program
        self.stats = stats
        # "program_hours_allocated" (in hours) from "stats" structure
        self.program_hours_allocated = float(self.stats[self.program]['program_hours_allocated'])
        # "total_hours_used" (in hours) from "stats" structure
        # This will be mutated by this call.
        self.total_hours_used = float(self.stats[self.program]['total_hours_used'])
        self.grid_times_targets = grid_times_targets
        self.debug = debug
        self.comment = ''

    def get_score(self, target, time):
        """
        Calculate the score for this constraint for a particular target and time.

        Parameters
        ----------
        target : ~mmtqueue.Target`
            An instance of the Target class.
        time : `~astropy.time.Time`
            The time to calculate the constraint.

        Returns
        -------
        score : float
            The constraint score as a float
        """
        score = 1.0 - (self.total_hours_used / self.program_hours_allocated)
        almost_zero = 0.1
        if score < almost_zero:
            score = almost_zero
            comment = 'Time not available in program.'
        elif score > 1.0:
            score = 1.0
            # comment = 'Time available in program.'
            comment = ''
        else:
            # comment = 'Time available in program.'
            comment = ''

        # New code for logging comments for special cases.
        try:
            if len(comment) > 0:
                self.process_comment(targets, times, score, comment)
        except:
            print("traceback: ", traceback.format_exc())

        return score

    def compute_constraint(self, times, observer, targets):
        """
         Actually do the real work of computing the constraint.

         Parameters
         ----------
         times : `~astropy.time.Time`
          The times to compute the constraint
         observer : `~astroplan.Observer`
          The observaton location from which to apply the constraints
         targets : sequence of `~astroplan.Target`
          The targets on which to apply the constraints.

         Returns
         -------
         mask_numpy : 1D or 2D Numpy array of float or bool
         The constraints, with targets along the first index and times along
         the second.
         """
        mask = []

        if targets.isscalar:
            target = targets
            for time in times:
                score = self.get_score(target, time)
                if False:
                    print("Scalar TimeAllocationConstraint score: ", repr(score))
                mask.append(score)
            if False:
                print("Handling as scalar")
            numpy_mask = np.array(mask)
        else:
            for target in targets:
                for time in times:
                    score = self.get_score(target, time)
                    if False:
                        print("Iterater TimeAllocationConstraint score: ", repr(score))
                    mask.append(score)
            if False:
                print("Handling as iterable")
            # Turn the mask into a numpy array and reshape.
            numpy_mask = np.reshape(np.array(mask), [len(targets), len(times)])

        if False:
            print("targets")
            print(repr(targets))
            print("times")
            print(repr(times))
            print("mask")
            print(repr(mask))
            print("numpy_mask")
            print(repr(numpy_mask))

        return numpy_mask
