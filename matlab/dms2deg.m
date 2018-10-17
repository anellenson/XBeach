function dms = dms2deg( deg, min, sec )

dms = abs(deg) + min/60 + sec/3600;
dms = dms * sign(deg);