function report=run_joint_leakage_reference_benchmark(folder)
%RUN_JOINT_LEAKAGE_REFERENCE_BENCHMARK Independent quadrature of two laws.
% Integrates a flow-weighted pore-radius distribution directly (not the
% production normal CDF) and variable-aperture Poiseuille resistance by
% adaptive numerical integration (not the production endpoint formula).
% This verifies the stated reduced laws; all parameters are illustrative.
% Local mixing/harmonic resistance is assumed, not persistent-radius tubes.
root=fileparts(mfilename('fullpath'));
if nargin<1,folder=fullfile(root,'leakage_research_20261001','joint_reference');end
if ~isfolder(folder),mkdir(folder);end
settings=joint_freeze_defaults(struct('enabled',true,'viscosityEnabled',false));
temperatures=[-8;-2;-.5;-.1;-.03;-.01;0;2];
actual=joint_freeze_conductance(temperatures,settings);reference=ones(size(actual));
beta=2*settings.interfacialEnergy_J_m2*settings.meltingTemperature_K/ ...
    (settings.iceDensity_kg_m3*settings.latentHeat_J_kg);
for j=1:numel(temperatures)
    if temperatures(j)<settings.bulkFreezingTemperature_C
        limit=log(beta/(settings.bulkFreezingTemperature_C-temperatures(j))/ ...
            settings.flowWeightedRadius_m)/settings.logRadiusStd;
        reference(j)=integral(@(z)exp(-z.^2/2)/sqrt(2*pi),-Inf,limit, ...
            'AbsTol',1e-25,'RelTol',1e-10);
    end
end
relative=abs(actual-reference)./max(reference,realmin);
capillary=table(temperatures,actual,reference,relative, ...
    'VariableNames',{'temperature_C','productionMobility','independentMobility','relativeDifference'});
weights=[.3 .7];radii=[6;5.9;5.55];T=[3 12];fraction=.001;
gaps=[.001 .004 .002;.002 .002 .002;.000001 .003 .002;0 .003 .002];
slit=zeros(size(gaps,1),4);
for j=1:size(gaps,1)
    gap=gaps(j,:);B=struct('nj',1,'nr',3,'nq',2,'nodalRadius_m',radii, ...
        'nodalGap0',gap','nodalGapU',speye(3),'radiusOuter_m',6, ...
        'radiusInner_m',5.55,'sampleWeights',weights);
    options=struct('openingTransmissionFraction',fraction,'baseMobilityFactor',0, ...
        'viscosityEnabled',true,'freezeEnabled',false);
    [computed,~,~,D]=joint_opening_p1_hydraulics(2e-10,B,zeros(3,1),T,options);
    expected=0;
    if all(gap>0)
        cuts=6-.45*[0 cumsum(weights)];resistance=0;
        for k=1:2
            mu=water_viscosity_iapws(T(k));
            integrand=@(r)mu./(fraction*interp1(flipud(radii),fliplr(gap)',r,'linear')).^3;
            resistance=resistance+integral(integrand,cuts(k+1),cuts(k), ...
                'Waypoints',radii(radii>cuts(k+1)&radii<cuts(k)), ...
                'RelTol',1e-10,'AbsTol',1e-3);
        end
        expected=1000*9.81/(12*resistance);
    end
    slit(j,:)=[computed expected abs(computed-expected)/max(expected,realmin) D.newslitGeometricVolume_m2];
end
slit=table((1:size(gaps,1))',slit(:,1),slit(:,2),slit(:,3),slit(:,4), ...
    'VariableNames',{'caseID','productionConductance_m_s','independentConductance_m_s','relativeDifference','geometricVolume_m2'});
passed=max(relative)<1e-7&&max(slit.relativeDifference)<1e-7;
report=struct('passed',passed,'capillary',capillary,'slit',slit, ...
    'settings',settings,'scope','Independent mathematical verification of illustrative capillary/slit resistance laws; NOT tunnel-gasket experiments or calibration.');
save(fullfile(folder,'joint_reference_report.mat'),'report','-v7');
writetable(capillary,fullfile(folder,'capillary_reference.csv'));
writetable(slit,fullfile(folder,'slit_reference.csv'));
assert(passed,'JointReference:Mismatch','Production laws disagree with independent resistance quadrature.');
end
