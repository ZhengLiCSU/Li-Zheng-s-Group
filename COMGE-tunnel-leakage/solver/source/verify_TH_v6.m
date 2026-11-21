function result=verify_TH_v6()
% Run invariant tests plus a nonfreezing analytical conduction benchmark.
% This verifies implementation, not field calibration or an ice-lens model.
result.invariants=TH_2D_conservative_v6('selftest');
errors=zeros(3,1); levels=[10 20 40]; dtDay=[1 .5 .25];
for k=1:3
    [y,z]=meshgrid(linspace(0,1,3),linspace(0,1,levels(k)+1));
    xy=[y(:),z(:)]; t=[];
    for j=1:2
        for i=1:levels(k)
            a=(j-1)*(levels(k)+1)+i; b=a+levels(k)+1;
            t=[t;a b a+1;b b+1 a+1]; %#ok<AGROW>
        end
    end
    S=struct('pAll',xy,'tSoil',t);
    o=struct('mesh',S,'durationDays',20,'dtInitial',dtDay(k)*86400,'dtMax',dtDay(k)*86400, ...
        'topMin',2,'topMax',2,'topHTC',inf,'geothermalFlux',0,'hydroFar',false, ...
        'topSeepage',false,'tunnelDrain',false,'Ss',1e-9,'plot',false,'saveOutput',false,'verbose',false);
    o.mechanics.enable=false;
    o.initial=struct('mode','fields','T',2+sin(pi*xy(:,2)/2),'h',xy(:,2)+1);
    r=TH_2D_conservative_v6(o);
    capacity=.65*2650*800+.35*1000*4180; conductivity=.65*2.3+.35*.57;
    exact=2+sin(pi*xy(:,2)/2).*exp(-(conductivity/capacity)*(pi/2)^2*20*86400);
    errors(k)=norm(r.T-exact,inf);
end
assert(all(diff(errors)<0),'Refinement did not reduce analytical conduction error.');
result.conduction=struct('verticalCells',levels,'dtDays',dtDay,'maxErrorC',errors);
disp(result.conduction);
save('TH_v6_verification.mat','result','-v7');
end
