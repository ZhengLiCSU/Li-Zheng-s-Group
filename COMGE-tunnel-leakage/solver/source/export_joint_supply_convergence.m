function result=export_joint_supply_convergence(studyFolder)
%EXPORT_JOINT_SUPPLY_CONVERGENCE Compare independent complete mesh histories.
% Refinement changes the local polygon geometry slightly. Each level has
% its own twenty-year history from the same prescribed initial profile.
% This reports two-grid sensitivity, not an asymptotic spatial error bound.
root=fileparts(mfilename('fullpath'));
if nargin<1, studyFolder=fullfile(root,'joint_supply_study_20260930'); end
folder=fullfile(studyFolder,'convergence_exports');
if ~isfolder(folder), mkdir(folder); end
names={'base','supplySeasonal'}; rows=table(); data=cell(2,2);
for k=1:2
    d=load(fullfile(studyFolder,names{k},'TH_segment_spinup_results.mat'),'result');
    f=load(fullfile(studyFolder,[names{k} 'Fine'],'TH_segment_spinup_results.mat'),'result');
    coarse=d.result; fine=f.result;
    assert(coarse.completedYears==20&&fine.completedYears==20);
    a=coarse.lastYear; b=fine.lastYear;
    assert(strcmp(a.version,b.version)&&strcmp(a.version,'TH_v14_20260930_joint_supply'));
    assert(a.parameters.Ksat==b.parameters.Ksat&& ...
        isequaln(a.parameters.forcing,b.parameters.forcing)&& ...
        isequaln(a.parameters.leakage,b.parameters.leakage));
    assert(a.mesh.jointRefinementReport.level==1&&b.mesh.jointRefinementReport.level==2);
    qa=1e3*a.leakage.totalVolume_m3_per_m; qb=1e3*b.leakage.totalVolume_m3_per_m;
    da=1e3*a.leakage.dailyTotalVolume_m3_per_m(:);
    db=1e3*b.leakage.dailyTotalVolume_m3_per_m(:);
    sa=1e3*sum(a.history(:,2)*86400.*a.leakage.hydraulics.unrestrictedRate_m2_per_s(:,7));
    sb=1e3*sum(b.history(:,2)*86400.*b.leakage.hydraulics.unrestrictedRate_m2_per_s(:,7));
    row=table(string(names{k}),qa,qb,100*abs(qa-qb)/qb, ...
        100*sum(abs(da-db))/qb,sa,sb,100*(1-qa/sa),100*(1-qb/sb), ...
        a.solverWork.rejectedSteps,b.solverWork.rejectedSteps, ...
        min(a.history(:,2)),min(b.history(:,2)), ...
        max(abs(a.history(:,8))),max(abs(b.history(:,8))), ...
        max(abs(a.history(:,9))),max(abs(b.history(:,9))), ...
        'VariableNames',{'caseName','Q_level1_L_per_m','Q_level2_L_per_m', ...
        'annualDifferenceRelativeToLevel2_percent','dailyL1RelativeToLevel2_percent', ...
        'sameT_reference_level1_L_per_m','sameT_reference_level2_L_per_m', ...
        'referenceReduction_level1_percent','referenceReduction_level2_percent', ...
        'year20Rejected_level1','year20Rejected_level2', ...
        'year20MinStep_level1_day','year20MinStep_level2_day', ...
        'year20MaxWaterDefect_level1_m2','year20MaxWaterDefect_level2_m2', ...
        'year20MaxEnergyDefect_level1_J_per_m','year20MaxEnergyDefect_level2_J_per_m'});
    if isempty(rows), rows=row; else, rows=[rows;row]; end %#ok<AGROW>
    data(k,:)={a,b};
end
result=struct('comparison',rows,'method', ...
    'Independent fresh twenty-year histories; two-grid sensitivity, not an error bound');
writetable(rows,fullfile(folder,'full_history_mesh_comparison.csv'));
save(fullfile(folder,'full_history_mesh_comparison.mat'),'result','-v7');
fig=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
clean=onCleanup(@()close(fig)); %#ok<NASGU>
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for k=1:2
    ax=nexttile(tl,k); hold(ax,'on');
    for j=1:2
        o=data{k,j};
        plot(ax,1:365,1e3*o.leakage.dailyTotalVolume_m3_per_m,'LineWidth',1.5);
    end
    grid(ax,'on'); xlabel(ax,'Phase day'); ylabel(ax,'Leakage (L day^{-1} m^{-1})');
    title(ax,[names{k} ': full-history local mesh sensitivity'],'Interpreter','none');
    legend(ax,{'Level 1: 1 cm','Level 2: 5 mm'},'Location','best');
    ax=nexttile(tl,k+2); hold(ax,'on');
    for j=1:2
        o=data{k,j}; h=o.history; hy=o.leakage.hydraulics;
        plot(ax,h(:,1)-19*365,hy.drivingHead_m(:,7),'LineWidth',1.5);
    end
    grid(ax,'on'); xlabel(ax,'Phase day'); ylabel(ax,'Actual joint driving head (m)');
    title(ax,sprintf('Annual difference %.3f%%; daily L1 %.3f%%', ...
        rows.annualDifferenceRelativeToLevel2_percent(k),rows.dailyL1RelativeToLevel2_percent(k)));
end
set(findall(fig,'Type','axes'),'FontSize',11);
title(tl,'Independent twenty-year histories on two locally refined meshes','FontSize',14);
exportgraphics(fig,fullfile(folder,'figure_full_history_mesh.png'),'Resolution',300);
exportgraphics(fig,fullfile(folder,'figure_full_history_mesh.pdf'),'ContentType','vector');
disp(rows);
end
