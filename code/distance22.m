function rng = distance22(lat1, lon1, lat2, lon2)

% Calculate great circle distance between points on a sphere using the
% Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
% length and has the same units as the radius of the sphere, R.  (If R is
% 1, then RNG is effectively arc length in radians.)

a = sind((lat2-lat1)/2).^2 + cosd(lat1) .* cosd(lat2) .* sind((lon2-lon1)/2).^2;
rng = 180/pi* 2 * atan2(sqrt(a),sqrt(1 - a));


