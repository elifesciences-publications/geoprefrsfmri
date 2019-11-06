function batch_connSocAffectcorr
%

%%  IC02_IC10
x_var = 'IC02_IC10';

nboot = 100000;
PLOTHIST = 0;
SAVERESID = 1;
YLIM = [0 15000];

% ADOSatscan
y_var = 'ADOS_CoSoTot';
Grp = 'GeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);

Grp = 'nonGeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);


%%  IC05_IC10
x_var = 'IC05_IC10';

nboot = 100000;
PLOTHIST = 0;
SAVERESID = 1;
YLIM = [0 15000];

% ADOSatscan
y_var = 'ADOS_CoSoTot';
Grp = 'GeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);

Grp = 'nonGeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);

%%  IC09_IC10
x_var = 'IC09_IC10';

nboot = 100000;
PLOTHIST = 0;
SAVERESID = 1;
YLIM = [0 15000];

% ADOSatscan
y_var = 'ADOS_CoSoTot';
Grp = 'GeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);

Grp = 'nonGeoASD';
result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM);


end % function




%%
function result = connSocAffectcorr(x_var,y_var,Grp,nboot,PLOTHIST,SAVERESID,YLIM)
warning off;

% initial stuff to set
fontSize = 14;
fontWeight = 'b';
lineColor = zeros(1,3);
lineWidth = 2;
cov_var = 'scan_age';
ci_type = 'per';

% read in data
rootpath = '/Users/mvlombardo/Dropbox/ACE_rsfMRI_Geo/reproAnalysis';
resultpath = fullfile(rootpath,'data','tidy');
datapath = fullfile(rootpath,'data','tidy');
datafile = fullfile(datapath,'data4connSocOrientcorr.txt');
D = readtable(datafile);

% make mask for specific group
if strcmp(Grp,'GeoASD')
    mask2use = ismember(D.subgrp2,'GeoASD');
    faceColor2use = [248 118 109]./255;
elseif strcmp(Grp,'nonGeoASD')
    mask2use = ismember(D.subgrp2,'nonGeoASD');
    faceColor2use = [167 143 6]./255;
end % if

% construct X and Y variables for analysis
X = [D.(cov_var) D.(x_var)]; X = X(mask2use,:);
Y = D.(y_var); Y = Y(mask2use);
idx2use = size(X,2);

% find correlation in sample
[x_tmp,y_tmp,r_tmp,p_tmp,se,meany,stats] = partialcor(X,Y,idx2use,0,1);

% make bootstrapping result reproducible
s=RandStream('mlfg6331_64');
options=statset('UseSubstreams', true, 'Streams', s);
reset(s);
% make bootstrapping result reproducible

% run bootstrapping
[ci,bootstat]  = bootci(nboot,{@robcorADOS,X,Y,idx2use}, 'options',options,'type',ci_type);

result.Grp = Grp;
result.x_var = x_var;
result.y_var = y_var;
result.cov_var = cov_var;
result.nboot = nboot;
result.ci_type = ci_type;
result.r = r_tmp;
result.pval = p_tmp;
result.ci = ci;
result.bootstat = bootstat;
warning on;

if PLOTHIST
    figure; set(gcf,'color','w');
    h = histogram(bootstat,100);
    ylim(YLIM);
    h.FaceColor = faceColor2use;
    h.EdgeColor = faceColor2use;
    grid on;
    hold on;
    l = line([r_tmp,r_tmp],ylim);
    l.Color = lineColor;
    l.LineWidth = lineWidth;
    ll = line([ci(1),ci(1)],ylim);
    ll.Color = zeros(1,3);
    ll.LineWidth = lineWidth;
    ll.LineStyle = ':';
    ul = line([ci(2),ci(2)],ylim);
    ul.Color = zeros(1,3);
    ul.LineWidth = lineWidth;
    ul.LineStyle = ':';
    xlabel('Correlation');
    ylabel('Frequency');
    xlim([-1 1]);

    set(gca,'fontsize',fontSize,'fontweight',fontWeight,'XTick',-1:0.2:1);
    fname2save = fullfile(plotpath,sprintf('%s_%s_%s_bootstrapCorr.pdf',Grp,x_var,y_var));
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end

fname2save = fullfile(resultpath,sprintf('%s_%s_%s_bootstrapCorr.mat',Grp,x_var,y_var));
save(fname2save,'result');

writetable(cell2table(num2cell([r_tmp, p_tmp, ci']),'VariableNames',{'r','p','ci95lo','ci95hi'}), ...
    fullfile(resultpath,sprintf('%s_%s_%s_corrRes.txt',Grp,x_var,y_var)),'FileType','text','delimiter','\t');

if SAVERESID
    writetable(cell2table([D.subjectId(mask2use), D.subgrp2(mask2use), num2cell([x_tmp y_tmp])],'VariableNames',{'subjectId','subgrp','x_adj','y_adj'}), ...
        fullfile(datapath,sprintf('%s_%s_%s_xy_adjdata.txt',Grp,x_var,y_var)), 'FileType','text','delimiter','\t');
end
end % function



%%
function r_tmp = robcorADOS(X,y,idx2use)
% rng(123);
warning off;
[x_tmp,y_tmp,r_tmp,p,se,meany,stats] = partialcor(X,y,idx2use,0,1);
warning on;
end
