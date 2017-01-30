# from vhelio calculation.ipynb
from astropy import coordinates

def helio_corr(t, loc, target=None):
    """
    Compute the heliocentric velocity correction at a given time and place.
    
    Paramters
    ---------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.
    target : SkyCoord or None
        The on-sky location at which to compute the correction.  If None,
        the function will return the cartesian vector instead of the radial
        velocity correction.
        
    Returns
    -------
    vcorr : astropy.units.Quantity with velocity units
        The heliocentric correction with a positive sign.  I.e., *add* this
        to an observed radial velocity to get the heliocentric velocity.
    """
    vsun = coordinates.get_body_barycentric_posvel('sun', t)[1]
    vearth = coordinates.get_body_barycentric_posvel('earth', t)[1]

    vsunearth = vearth - vsun
    
    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)
    
    vsuntarg = (vsunearth.xyz + gcrs_v.xyz).to(u.km/u.s)
    if target is None:
        return vsuntarg
    else:
        gtarg = target.transform_to(coordinates.GCRS(obstime=t, obsgeoloc=gcrs_p))
        targxyz = gtarg.represent_as(coordinates.UnitSphericalRepresentation).to_cartesian().xyz
        return coordinates.matrix_utilities.matrix_product(vsuntarg, targxyz)

