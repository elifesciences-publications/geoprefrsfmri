function estimateConnectivity
%
%

%% initial stuff to set
addpath ~/Dropbox/matlab/FSLNets;
addpath ~/Dropbox/matlab/L1precision;
addpath ~/Dropbox/matlab/pwling;
addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')));
ncomps = 30;
rho = 1;

rootpath = '/Users/mvlombardo/Dropbox/ACE_rsfMRI_Geo/';
ica_dir = fullfile(rootpath,'reproAnalysis','ica30_ASDTDLDDDSIB.ica');
dr_dir = fullfile(ica_dir,'dualreg');
ts_dir = fullfile(dr_dir,'dr_stage1');
result_dir = fullfile(rootpath,'reproAnalysis','data','tidy');
scanInfoFile = fullfile(result_dir,'final_allETrsfMRIsubs_phenodata04_ASDTDLDDDSIB.txt');

%% Set components to use
comps2use = [2 4 5 6 9 10 11 21 26 28]; 

%% load scan information
S = readtable(scanInfoFile,'delimiter','\t');
p2f2 = S.p2f2_scannum;
subID = S.subjectId;
subgrpDx = S.subgrp2;
subgrpDx_geomidsoctd = S.subgrp;
Dx = S.Dx;
CaseControl = S.CaseControl;
sex = S.sex;
FixSoc = S.Percent_Fixation_Social;
FixGeo = S.Percent_Fixation_Geometric;
ETage = S.ET1_Age;
scan_age = S.scan_age;

% Component Names
for icomp = 1:length(comps2use)
    compNames{icomp} = sprintf('IC%02d',comps2use(icomp));
end % for icomp


%% grab timeseries data
d = dir(fullfile(ts_dir,'*.txt'));
nframeCensored = zeros(length(d),1);
for isub = 1:length(d)
    % filename
	tsfile{isub,1} = fullfile(ts_dir,d(isub).name);
    % read in data and extract
	a = readtable(tsfile{isub},'delimiter','space','HeaderLines',0, ...
        'ReadVariableNames',false,'FORMAT',repmat('%f',1,ncomps*2));
	tmp = table2array(a(:,1:2:end-1));

    ts_data{isub,1} = tmp(:,comps2use);

    tmp_r = nets_netmats(ts_data{isub},1,'ridgep',rho);

    tmp_z = r2z_transform(tmp_r);
	r(:,:,isub) = tmp_r;
    z(:,:,isub) = tmp_z;
    % grab lower triangle of correlation matrix
    lower_tri_mask = tril(tmp_r,-1)~=0;
    rts_lowertri(isub,:) = tmp_r(lower_tri_mask);
    ts_lowertri(isub,:) = tmp_z(lower_tri_mask);
end % for isub

% make connection names to label the columns of the rts_lowertri matrix
idx2use = find(lower_tri_mask);
[x,y] = ind2sub(size(tmp_r),idx2use);
for i = 1:length(idx2use)
    colLabels{i} = sprintf('%s_%s',compNames{y(i)},compNames{x(i)});
end % for i

fname2save = fullfile(result_dir,sprintf('partialCorDataASDTDLDDDSIB_ridge_lambda%d.txt',rho));

% make table for rts_lowertri and export to text file
tab2write=cell2table([subID CaseControl subgrpDx sex num2cell([scan_age,FixGeo,ETage]) num2cell(rts_lowertri)], ...
    'VariableNames',[{'subjectId','CaseControl','subgrp','sex','scan_age','FixGeo','ETage'},colLabels]);
writetable(tab2write,fname2save,'FileType','text','delimiter','\t');

end % function estimateConnectivity

%%
function Z = r2z_transform(R)

Z = zeros(size(R));

% Loop over columns in R
for icolumn = 1:size(R,2)
    % Loop over rows in R
    for irow = 1:size(R,1)
        % this converts r-values to z-scores
        Z(irow,icolumn) = 0.5*log((1+R(irow,icolumn))/(1-R(irow,icolumn)));
    end % for irow
end % for icolumn

end % function r2z_transform
