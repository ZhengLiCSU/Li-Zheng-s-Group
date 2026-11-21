function layout=hsr_eight_segment_layout(rotationRad)
%HSR_EIGHT_SEGMENT_LAYOUT Eight-joint 5+2+1 layout for a reference ring.
% Angles follow the Shiziyang prototype drawing, not its diameter, water
% pressure, climate, or a claim that all rings share one orientation.
% Source: He et al., China Engineering Science 14(10), 2012, Fig. 2.
% https://www.engineering.org.cn/sscae/CN/1160098445321102171
% The image has the key-block centre at 180 degrees (left springline) in
% this model's positive-down mesh coordinates. Rotate the entire ring with
% rotationRad to study another assembly orientation.
if nargin<1, rotationRad=0; end
assert(isscalar(rotationRad)&&isreal(rotationRad)&&isfinite(rotationRad), ...
    'rotationRad must be a finite real scalar.');
keyDeg=16+21/60+49.09/3600;
otherDeg=(360-keyDeg)/7; % exact closure; published angle is rounded
keyLeftDeg=180-keyDeg/2;
baseDeg=sort(mod([keyLeftDeg,180+keyDeg/2+(0:6)*otherDeg],360));
baseJointLabels={'B3-B2','B2-B1','B1-L1','L1-F', ...
    'F-L2','L2-B5','B5-B4','B4-B3'};
baseSegmentLabels={'B2','B1','L1','F','L2','B5','B4','B3'};
[angles,order]=sort(mod(baseDeg*pi/180+rotationRad,2*pi));
layout=struct('segmentCount',8,'jointAnglesRad',angles, ...
    'jointAnglesDeg',angles*180/pi, ...
    'jointLabels',{baseJointLabels(order)}, ...
    'segmentLabels',{baseSegmentLabels(order)}, ...
    'keyBlockAngleDeg',keyDeg,'otherBlockAngleDeg',otherDeg, ...
    'keyCenterAngleRad',mod(pi+rotationRad,2*pi), ...
    'rotationRad',rotationRad, ...
    'source','https://www.engineering.org.cn/sscae/CN/1160098445321102171');
assert(abs(sum([keyDeg,7*otherDeg])-360)<1e-10);
end
