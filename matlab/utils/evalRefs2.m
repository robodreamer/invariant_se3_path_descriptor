function [tRef, xRef, vRef, aRef] = evalRefs2(pp, tRate, tRatePrev, h)

persistent tLast
if isempty(tLast)
    tLast = pp.x.breaks(1);
end

tRef = tLast + tRate*h;

if tRef > pp.x.breaks(end)
    tRef = [];
    xRef = [];
    vRef = [];
    aRef = [];
else
    tRateDot = (tRate - tRatePrev)/h;
    xRef = ppval(pp.x, tRef);
    vRef = ppval(pp.v, tRef)*tRate;
    aRef = ppval(pp.v, tRef)*tRateDot + ppval(pp.a, tRef)*tRate*tRate;
    tLast = tRef;
end

end