function files=plot_v16_verification(caseName)
%PLOT_V16_VERIFICATION Actual short-run THM fields; no annual extrapolation.
if nargin<1, caseName='real_mesh_profile'; end
root=fileparts(mfilename('fullpath'));
folder=fullfile(root,'v16_verification_20260930',caseName);
d=load(fullfile(folder,'thm_results.mat'),'out'); o=d.out;
xy=o.pAll; tri=o.mesh.tSoil; M=o.mechanics.mesh; c=o.mechanics.contact;
f=figure('Visible','off','Color','w','Position',[80 80 1450 900]);
cleaner=onCleanup(@()close(f)); %#ok<NASGU>
lay=tiledlayout(f,2,3,'TileSpacing','loose','Padding','loose');
lay.OuterPosition=[0 0 1 .93];
annotation(f,'textbox',[.01 .95 .98 .04], ...
    'String',sprintf('Contact THM verification: day %.2f (illustrative parameters)',o.timeDay), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',14,'EdgeColor','none');
a=nexttile(lay);
patch(a,'Faces',tri,'Vertices',xy,'FaceVertexCData',o.T, ...
    'FaceColor','interp','EdgeColor','none'); hold(a,'on');
for name={'qTunnelGlobal','qInvertGlobal','tInvertTriGlobal'}
    if isfield(o.mesh,name{1})&&~isempty(o.mesh.(name{1}))
        patch(a,'Faces',o.mesh.(name{1}),'Vertices',xy,'FaceVertexCData',o.T, ...
            'FaceColor','interp','EdgeColor','none');
    end
end
domainAxes(a,xy); title(a,'(a) Temperature (^{\circ}C)'); colorbar(a);
a=nexttile(lay); dh=o.h-o.initialX(numel(o.T)+(1:numel(o.h)));
nodal=zeros(size(xy,1),1); nodal(o.soilNodes)=dh;
patch(a,'Faces',tri,'Vertices',xy,'FaceVertexCData',nodal, ...
    'FaceColor','interp','EdgeColor','none');
domainAxes(a,xy); title(a,'(b) Change in liquid pressure head (m)'); colorbar(a);
a=nexttile(lay);
patch(a,'Faces',tri,'Vertices',xy,'FaceVertexCData',100*o.state.elementVolumetricStrain, ...
    'FaceColor','flat','EdgeColor','none');
domainAxes(a,xy); title(a,'(c) Incremental volumetric strain (%)'); colorbar(a);
a=nexttile(lay); colors=lines(8); hold(a,'on');
for j=1:8
    pair=M.jointNodePairs(:,1,j);
    radius=hypot(M.p(pair,1)-M.center(1),M.p(pair,2)-M.center(2));
    depth=1000*(max(radius)-radius);
    width=1.1; if j==7, width=2.8; end
    style='-o'; if j==8, style='--s'; end
    plot(a,depth,1e6*c.jointGap_m(:,j),style,'Color',colors(j,:), ...
        'LineWidth',width,'MarkerSize',3,'DisplayName',sprintf('Joint %d',j));
end
yline(a,0,':k','HandleVisibility','off'); grid(a,'on');
xlabel(a,'Distance from outer face (mm)'); ylabel(a,'Signed mechanical gap (\mum)');
title(a,'(d) Gap profile; positive means separation');
legend(a,'Location','best','NumColumns',2,'FontSize',8);
a=nexttile(lay);
bar(a,1:8,max(c.jointContactPressure_Pa,[],1)/1e6, ...
    'FaceColor',[.25 .42 .64]); grid(a,'on'); xticks(a,1:8);
xlabel(a,'Joint ID'); ylabel(a,'Maximum normal contact pressure (MPa)');
title(a,'(e) Compression-only contact');
a=nexttile(lay);
s=squeeze(mean(c.stressGauss_Pa,2));
vm=sqrt(.5*((s(:,1)-s(:,2)).^2+(s(:,2)-s(:,4)).^2+ ...
    (s(:,4)-s(:,1)).^2)+3*s(:,3).^2)/1e6;
patch(a,'Faces',M.q,'Vertices',M.p,'FaceVertexCData',vm, ...
    'FaceColor','flat','EdgeColor',[.4 .4 .4],'LineWidth',.2);
axis(a,'equal'); set(a,'YDir','reverse'); box(a,'on');
xlabel(a,'x (m)'); ylabel(a,'Depth (m)'); colorbar(a);
title(a,'(f) Ring von Mises stress (MPa)');
set(findall(f,'Type','axes'),'FontName','Arial','FontSize',10);
colormap(f,turbo(256));
files=struct('png',fullfile(folder,'thm_verification_fields.png'), ...
    'pdf',fullfile(folder,'thm_verification_fields.pdf'));
exportgraphics(f,files.png,'Resolution',180);
exportgraphics(f,files.pdf,'ContentType','vector');
disp(files);
end

function domainAxes(a,xy)
axis(a,'equal'); set(a,'YDir','reverse'); box(a,'on');
xlim(a,[min(xy(:,1)),max(xy(:,1))]); ylim(a,[min(xy(:,2)),max(xy(:,2))]);
xlabel(a,'x (m)'); ylabel(a,'Depth (m)');
end
