# from vhelio calculation.ipynb

def helio_corr(t, loc, target=None):
    """
    The thing you add to a target to get the heliocentric velocity
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

