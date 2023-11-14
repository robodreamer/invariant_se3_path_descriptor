function [tRef, xRef, vRef] = evalRefs(pp, tRate, h)

persistent tLast
if isempty(tLast)
    tLast = pp.x.breaks(1);
end

tRef = tLast + tRate*h;

if tRef > pp.x.breaks(end)
    tRef = [];
    xRef = [];
    vRef = [];
else
    xRef = ppval(pp.x, tRef);
    vRef = ppval(pp.v, tRef)*tRate;
    tLast = tRef;
end

end