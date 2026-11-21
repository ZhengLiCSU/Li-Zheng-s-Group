function files=export_joint_supply_benchmark(inputFolder,outputFolder)
%EXPORT_JOINT_SUPPLY_BENCHMARK Plot actual warm steady benchmark points.
% No TH solves run here. Read the 27 mesh/K/C linear benchmark calculations
% and their pairwise mesh changes; display computed points without an
% interpolated contour. These warm steady rates are not year-20 leakage.
root=fileparts(mfilename('fullpath'));
if nargin<1||isempty(inputFolder)
 inputFolder=fullfile(root,'joint_supply_study_20260930','benchmark');
 if ~isfile(fullfile(inputFolder,'warm_linear_supply.csv'))
  inputFolder=fullfile(root,'joint_supply_study_20260930','warm_linear_benchmark');
 end
end
if nargin<2||isempty(outputFolder), outputFolder=inputFolder; end
R=readtable(fullfile(inputFolder,'warm_linear_supply.csv'));
V=readtable(fullfile(inputFolder,'warm_linear_mesh_convergence.csv'));
assert(height(R)==27&&isequal(unique(R.meshLevel),(0:2)'), ...
 'Expected 27 computed cases on mesh levels 0, 1 and 2.');
assert(height(V)==18&&all(isfinite(R.Q_m2_per_s))&& ...
 all(R.C_times_staticHead_m2_per_s>0),'Benchmark tables are incomplete.');
if ~isfolder(outputFolder), mkdir(outputFolder); end
files=struct('PNG',fullfile(outputFolder,'figure_warm_supply_benchmark.png'), ...
 'PDF',fullfile(outputFolder,'figure_warm_supply_benchmark.pdf'), ...
 'pointCSV',fullfile(outputFolder,'warm_benchmark_plot_points.csv'), ...
 'notesTXT',fullfile(outputFolder,'READ_BENCHMARK_FIGURE.txt'));
P=R(:,{'meshLevel','K_m_per_s','C_m_per_s','Lambda', ...
 'Q_m2_per_s','C_times_staticHead_m2_per_s','negativeTransmissivityEdges'});
P.actualOverStatic=R.Q_m2_per_s./R.C_times_staticHead_m2_per_s;
P.uniformLoadSeriesApprox=1./(1+R.Lambda);
writetable(P,files.pointCSV,'Encoding','UTF-8');
colors=[.50 .50 .50;.10 .40 .75;.25 .60 .35]; markers={'s','o','^'};
labels={'Original mesh','Active-joint refinement 1','Active-joint refinement 2'};
fig=figure('Visible','off','Color','w','Position',[80 80 1450 950]);
cleanFig=onCleanup(@() close(fig)); %#ok<NASGU>
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(t); F=sortrows(R(R.meshLevel==2,:),{'K_m_per_s','C_m_per_s'});
ratio=F.Q_m2_per_s./F.C_times_staticHead_m2_per_s;
scatter(ax,log10(F.C_m_per_s),log10(F.K_m_per_s),165,ratio,'filled', ...
 'MarkerEdgeColor',[.15 .15 .15]);
for k=1:height(F)
 text(ax,log10(F.C_m_per_s(k)),log10(F.K_m_per_s(k))+.09, ...
  sprintf('%.3f',ratio(k)),'HorizontalAlignment','center','FontSize',10, ...
  'FontName','Arial');
end
colormap(ax,parula(256)); clim(ax,[0 1]); cb=colorbar(ax);
cb.Label.String='Q_{steady}/(C H_{reference})'; cb.FontSize=10;
xlim(ax,[min(log10(F.C_m_per_s))-.45,max(log10(F.C_m_per_s))+.45]);
ylim(ax,[min(log10(F.K_m_per_s))-.35,max(log10(F.K_m_per_s))+.35]);
xticks(ax,log10(unique(F.C_m_per_s))); yticks(ax,log10(unique(F.K_m_per_s)));
xlabel(ax,'log_{10} joint reference C (m s^{-1})');
ylabel(ax,'log_{10} uniform saturated soil K (m s^{-1})');
title(ax,'(a) Nine actual calculations on refinement level 2'); formatAxes(ax);

ax=nexttile(t); hold(ax,'on');
for level=0:2
 mask=R.meshLevel==level;
 scatter(ax,R.Lambda(mask),R.Q_m2_per_s(mask)./R.C_times_staticHead_m2_per_s(mask), ...
  46,colors(level+1,:),markers{level+1},'filled','DisplayName',labels{level+1});
end
xx=logspace(log10(min(R.Lambda))-.2,log10(max(R.Lambda))+.2,200);
plot(ax,xx,1./(1+xx),'k--','LineWidth',1.3, ...
 'DisplayName','Uniform-band load approximation: 1/(1+\Lambda)');
set(ax,'XScale','log'); ylim(ax,[0 1.05]);
xlabel(ax,'\Lambda=C/G_{soil}, G_{soil} from independent Laplace operator');
ylabel(ax,'Q_{steady}/(C H_{reference})');
title(ax,'(b) All 27 points and the resistance approximation');
legend(ax,'Location','southwest','FontSize',9); formatAxes(ax);

ax=nexttile(t);
pairs=sortrows(unique([V.K_m_per_s,V.C_m_per_s],'rows'),[1 2]);
change=zeros(size(pairs,1),2); pairLabels=cell(size(pairs,1),1);
for k=1:size(pairs,1)
 for fine=1:2
  mask=V.K_m_per_s==pairs(k,1)&V.C_m_per_s==pairs(k,2)&V.fineLevel==fine;
  assert(nnz(mask)==1,'Mesh comparison is missing a K/C pair.');
  change(k,fine)=V.rateChangeRelativeToFine_percent(mask);
 end
 % XTICKLABELS may split newline-containing cell labels into separate
 % ticks. Keep each K/C pair as one complete single-line tick label.
 pairLabels{k}=sprintf('K%g/C%g',pairs(k,1),pairs(k,2));
end
bars=bar(ax,1:size(pairs,1),change,'grouped');
bars(1).FaceColor=[.55 .55 .55]; bars(2).FaceColor=[.10 .40 .75];
set(ax,'YScale','log');
xticks(ax,1:size(pairs,1)); xticklabels(ax,pairLabels); xtickangle(ax,40);
xlabel(ax,'Like-for-like soil and joint parameter pairs (m s^{-1})');
ylabel(ax,'|Q_{fine}-Q_{coarse}|/Q_{fine} (%)');
title(ax,'(c)  Computed leakage changes between meshes','Interpreter','none');
legend(ax,bars,{'Original \rightarrow refinement 1','Refinement 1 \rightarrow 2'}, ...
 'Location','best','FontSize',9); formatAxes(ax); ax.FontSize=9;

ax=nexttile(t); hold(ax,'on');
% G/K is geometric and independent of the permeability scale. Show all
% observed mesh-level values, with their exact min/max rather than averaging.
gMin=zeros(3,1); gMax=gMin; edge=zeros(3,1);
for level=0:2
 m=R.meshLevel==level;
 shape=R.uniformBandSoilConductance_m_per_s(m)./R.K_m_per_s(m);
 gMin(level+1)=min(shape); gMax(level+1)=max(shape);
 edge(level+1)=max(R.maximumSupportedInterfaceEdge_m(m));
end
yyaxis(ax,'left');
plot(ax,0:2,gMin,'-o','Color',[.10 .40 .75],'LineWidth',1.6,'MarkerSize',7);
plot(ax,0:2,gMax,':','Color',[.10 .40 .75],'LineWidth',1,'HandleVisibility','off');
ylabel(ax,'Geometric soil conductance G_{soil}/K');
yyaxis(ax,'right');
semilogy(ax,0:2,edge,'-s','Color',[.90 .50 .10],'LineWidth',1.6,'MarkerSize',6);
ylabel(ax,'Maximum interface edge touching joint band (m)');
xticks(ax,0:2); xlim(ax,[-.15 2.15]); xlabel(ax,'Local mesh refinement level');
title(ax,'(d) Resolution of the 0.02 m numerical joint support'); formatAxes(ax);
title(t,'Warm saturated steady benchmark: soil supply, joint resistance and mesh resolution', ...
 'FontName','Arial','FontSize',15,'FontWeight','normal');
exportgraphics(fig,files.PNG,'Resolution',300);
exportgraphics(fig,files.PDF,'ContentType','vector');
fid=fopen(files.notesTXT,'w'); assert(fid>=0,'Cannot write benchmark figure notes.');
cleanFile=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,['The figure contains 27 actual warm saturated steady linear hydraulic solves.\n' ...
 'Panel (a) shows nine computed level-2 points, without an interpolated phase contour.\n' ...
 'Gsoil is independently calculated from a normalized uniform-band load and the soil Laplace matrix.\n' ...
 '1/(1+Lambda) is a uniform-load scalar resistance approximation, not the exact finite-band Robin solution.\n' ...
 'Panel (c) uses |Qfine-Qcoarse|/Qfine; refinement level 2 is a reference, not an exact solution.\n' ...
 'Original negative cotangent weights were not clipped. Refined geometry is slightly changed by circular node insertion.\n' ...
 'These are warm steady verification cases, not frozen TH histories or year-20 results.\n']);
end

function formatAxes(ax)
set(ax,'FontName','Arial','FontSize',11,'LineWidth',.8,'Box','on','TickDir','out');
grid(ax,'on'); ax.GridAlpha=.15;
end
