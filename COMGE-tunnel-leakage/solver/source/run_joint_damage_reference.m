function report=run_joint_damage_reference(folder)
%RUN_JOINT_DAMAGE_REFERENCE Independent radius-moment and phase quadrature.
% Checks the number-weighted r^4 and r^2 lognormal integrals directly.
% This is mathematical verification of an illustrative connected network,
% not an experiment on EPDM or a freeze-thaw damage evolution calibration.
arguments
    folder (1,1) string = fullfile(fileparts(mfilename('fullpath')),'joint_damage_20261002','joint_reference')
end
if ~isfolder(folder),mkdir(folder);end
L=struct('jointAnglesRad',0,'jointConductance',2e-10, ...
    'jointFreeze',joint_freeze_defaults(struct('enabled',true,'viscosityEnabled',false)), ...
    'jointStorage',struct('enabled',true,'tortuosity',1), ...
    'jointDamage',struct('enabled',true,'severity',0));
D=[0 .25 .5 .75 1];T=[-2 -.5 -.1 -.03 -.01 0 2];
rows=zeros(numel(D)*numel(T),8);n=0;
s=L.jointFreeze.logRadiusStd;rq=L.jointFreeze.flowWeightedRadius_m;
m=log(rq)-4*s^2;
beta=2*L.jointFreeze.interfacialEnergy_J_m2*L.jointFreeze.meltingTemperature_K/ ...
    (L.jointFreeze.iceDensity_kg_m3*L.jointFreeze.latentHeat_J_kg);
M4=integral(@(z)exp(4*s*z-z.^2/2)/sqrt(2*pi),-Inf,Inf,'RelTol',1e-11);
M2=integral(@(z)exp(2*s*z-z.^2/2)/sqrt(2*pi),-Inf,Inf,'RelTol',1e-11);
for j=1:numel(D)
    L.jointDamage.severity=D(j);[law,geometry]=joint_microchannel_damage(L);
    lambda=geometry.radiusMultiplier;
    for k=1:numel(T)
        upper=Inf;
        if T(k)<0,upper=(log(beta/(-T(k)))-log(lambda)-m)/s;end
        fraction4=integral(@(z)exp(4*s*z-z.^2/2)/sqrt(2*pi),-Inf,upper, ...
            'RelTol',1e-10,'AbsTol',1e-25)/M4;
        fraction2=integral(@(z)exp(2*s*z-z.^2/2)/sqrt(2*pi),-Inf,upper, ...
            'RelTol',1e-10,'AbsTol',1e-25)/M2;
        actual=law.jointConductance/L.jointConductance*joint_freeze_conductance(T(k),law.jointFreeze);
        expected=lambda^4*fraction4;
        [~,~,v]=joint_channel_inventory(T(k),law.jointFreeze,law.jointStorage, ...
            law.jointConductance,.45,1000,3.34e5,4180,2100,2500*900,9.81);
        n=n+1;rows(n,:)=[D(j),T(k),actual,expected, ...
            abs(actual-expected)/max(expected,realmin),v.liquidVolumeFraction, ...
            fraction2,abs(v.liquidVolumeFraction-fraction2)];
    end
end
comparison=array2table(rows,'VariableNames',{'damageLevel','temperature_C', ...
    'productionConductanceRatio','independentConductanceRatio','relativeConductanceError', ...
    'productionLiquidVolumeFraction','independentLiquidVolumeFraction','absoluteLiquidFractionError'});
report=struct('passed',max(rows(:,5))<1e-7&&max(rows(:,8))<1e-10, ...
    'comparison',comparison,'maxRelativeConductanceError',max(rows(:,5)), ...
    'maxAbsoluteLiquidFractionError',max(rows(:,8)), ...
    'scope','Independent number-weighted radius quadrature; illustrative fixed hydraulic geometry, no damage evolution or experimental calibration.');
save(fullfile(folder,'joint_damage_reference.mat'),'report','-v7');
writetable(comparison,fullfile(folder,'joint_damage_reference.csv'));
assert(report.passed,'JointDamageReference:Failed','Independent physical-moment comparison failed.');
end
