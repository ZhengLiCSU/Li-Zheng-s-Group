function result=verify_TH_v7()
% Run invariant tests plus a nonfreezing analytical conduction benchmark.
% This verifies implementation, not field calibration or an ice-lens model.
result.invariants=TH_2D_conservative_v7('selftest');
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
    r=TH_2D_conservative_v7(o);
    capacity=.65*2650*800+.35*1000*4180; conductivity=.65*2.3+.35*.57;
    exact=2+sin(pi*xy(:,2)/2).*exp(-(conductivity/capacity)*(pi/2)^2*20*86400);
    errors(k)=norm(r.T-exact,inf);
end
assert(all(diff(errors)<0),'Refinement did not reduce analytical conduction error.');
result.conduction=struct('verticalCells',levels,'dtDays',dtDay,'maxErrorC',errors);
disp(result.conduction);
% Work-based time-step recovery and output-boundary truncation regression.
S=struct('pAll',[0 0;1 0;0 1;1 1],'tSoil',[1 2 3;2 4 3]);
o=struct('mesh',S,'durationDays',.02,'dtInitial',1,'dtMin',1,'dtMax',100, ...
    'outputEveryDays',.003,'topMin',2,'topMax',2,'topHTC',inf, ...
    'geothermalFlux',0,'hydroFar',false,'topSeepage',false,'tunnelDrain',false, ...
    'plot',false,'saveOutput',false,'verbose',false);
o.mechanics.enable=false;
o.initial=struct('mode','fields','T',2*ones(4,1),'h',S.pAll(:,2)+1);
q=TH_2D_conservative_v7(o);
assert(max(q.history(:,2))*86400>90, ...
    'Time step failed to recover after a small initial/output-boundary step.');
result.timeController=struct('passed',true,'acceptedSteps',size(q.history,1), ...
    'maxDtSeconds',max(q.history(:,2))*86400);
save('TH_v7_verification.mat','result','-v7');
end
