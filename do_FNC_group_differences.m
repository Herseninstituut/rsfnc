%% do_FNC_group_differences.m
%
% Companion code to:
% Cerliani L., Mennes M., Thomas RM, DI Martino A, Thioux M. Keysers C. (2015)
% "Increased functional connectivity between subcortical and cortical 
% resting-state networks in autism spectrum disorder"; JAMA Psychiatry
%
% PLEASE READ CAREFULLY THE FOLLOWING LINES.
% EVERY USE OF THIS SOFTWARE IS CONDITIONED TO THE ACCEPTANCE OF THE
% DISCLAIMER REPORTED BELOW.
%
% This script performs the estimation of group differences in Functional 
% Network Connectivity (FNC) between participants with ASD and typically 
% developing (TD) participants which is reported in the manuscript. 
%
% The FSL software library (http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) must
% be installed and its binary must be in the current path in order to run 
% the script. In particular, the Matlab(r) binaries provided in the 
% $FSLDIR/etc/matlab directory must be in the path, as well as all the
% binaries in $FSLDIR/bin 
%
% IT IS STRONGLY ADVISED TO EVALUATE ONE CELL AT THE TIME (using CMD+ENTER 
% on Mac), AND TO INSPECT THE RESULTS ACCORDING TO THE COMMENTS REPORTED 
% IN THE SCRIPT, RATHER THAN RUNNING THE WHOLE SCRIPT AT ONCE
%
% Input: subject-specific network-specific time courses derived from the 
% dual-regression stage 1 (i.e. spatial regression) performed on the 
% results of the metaICA. These time courses are provided in the 'data'
% directory, and refer to the 359 participants from the ABIDE database
% (http://fcon_1000.projects.nitrc.org/indi/abide/) whose rs-fmri data were
% selected for the analyses described in the manuscript.
%
% Output: matrices displaying significant (p<.05 FDR corrected) differences
% in FNC between the ASD and the TD sample. Positive differences indicate
% higher FNC in the ASD group.
%
% This script requires additional software to run. This additional software
% is provided in the 'addons' directory. Specifically:
%
% do_randomise_wrapper.m is a simple wrapped for the 'randomise' software
% provided with FMRIB FSL software (http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).
% This also means that FSL must be installed and its binaries must be in
% the path in order to run the script.
%
% FDR_nichols.m is the FDR.m software developed by Tom Nichols to estimate
% the FDR correction for multiple comparisons. It can be found at 
% http://www-personal.umich.edu/~nichols/FDR/FDR.m
%
% xticklabel_rotate.m was developed by Brian Katz, and can be downloaded
% from:
% http://fr.mathworks.com/matlabcentral/fileexchange/3486-xticklabel-rotate
% The relative copyright notive, list of conditions and disclaimer are to
% be found in the 'addons' directory
%
% notBoxPlot.m (and the associated SEM_calc.m and tInterval_Calc.m) was
% developed by Rob Campbell, and can be downloaded from:
% http://www.mathworks.com/matlabcentral/fileexchange/26508-notboxplot-alternative-to-box-plots
% The relative copyright notive, list of conditions and disclaimer are to
% be found in the 'addons' directory
%
% The software is provided 'AS IS'. The public distribution of this script
% is UNIQUELY intended to provide an accurate documentation of the analyses 
% performed in the published manuscript. The author declines any
% responsibility regarding machine malfunctioning or data loss due to the
% use of this script in the version provided here or in any modification
% of it.
%
% Author: Leonardo Cerliani
% leonardo.cerliani@gmail.com
% June 7th, 2015


%% setup environmental variables

clear all; close all; clc

basedir=pwd;
cd(basedir)

% add additional software to the path
addpath(genpath([basedir '/addons']));

% directory where the output of dr_stage1 are stored
dr_dir=[basedir '/data/dr_stage1'];


% Load the phenotypical information for the subset selected for the
% analysis.
% These data are extracted from the Composite Phenotypic File 
% (Phenotypic_V1_0b.csv) available on the ABIDE database
load([basedir '/data/pheno_359.mat']);

% read the participants' IDs. This list matches the order of phenomat(:,2)
IDs = dlmread([basedir '/data/list_359_IDs']);


% Define which of the 52 ICs are resting state networks,
% here indexed by the variable 'fn_relevant'
fn_relevant = [1 3 5 8 9 10 13 15 16 17 19 21 23 24 25 27 29 30 33];

nIC_fn = length(fn_relevant);

 
labels = {...
'IC1 - V1' ...
'IC3 - pCRB' ...
'IC5 - dSI+dM1+mPMC' ...
'IC8 - pfus+V1' ...
'IC9 - BA32sg+bOFC' ...
'IC10 - BA32+9+10'...
'IC13 - aCRB' ...
'IC15 - dPcun'...
'IC16 - TE+pSTG+PARop+pIC'...
'IC17 - Bg+Th'...
'IC19 - OCC pole'...
'IC21 - IPS+vPMC(RH)'...
'IC23 - Saliency'...
'IC24 - STS+IFG(LH)'...
'IC25 - LH FTP'...
'IC27 - DMN'...
'IC29 - vSI+vM1+pIC'...
'IC30 - IFG+pIPS(LH)'...
'IC33 - TP+Amy+Pons'...
};


% Assign the values in phenomat to single variables
SITE_ID = phenomat(:, 1);
SUB_ID = phenomat(:, 2);
DSM_IV_TR = phenomat(:, 4);
AGE_AT_SCAN = phenomat(:, 5);
SEX = phenomat(:, 6);
FIQ = phenomat(:, 7);
EYE_STATUS_AT_SCAN = phenomat(:,16);
TRs = phenomat(:,17);
FD_mean = phenomat(:,18);



%% Create one index for ASD and one for TD
idx_ASD = find(DSM_IV_TR == 1 | DSM_IV_TR == 2); % 1=Autism, 2=Asperger
idx_TD = find(DSM_IV_TR == 0);

clc
disp(['There are ' num2str(length(idx_ASD)) ' ASD']);
disp(['There are ' num2str(length(idx_TD)) ' TD']);
disp(['In total ' num2str(length(idx_ASD) + length(idx_TD)) ' participants'])

nsubj = length(idx_ASD) + length(idx_TD);




%% Estimate the FNC matrix for each participant
%
%  !!! IMPORTANT !!!
%  
%  We noticed some very slight variation (2nd-3rd decimal) in the estimated 
%  FNC values in different machines. 
%
%  These variations do not influence the results published in the paper.
%
%  For this reason we provide here the values which were calculated on the
%  machine where all the analyses published in the paper were performed, in
%  order to reproduce the exact values reported in the manuscript.
%
%  If you wish to perform again all the calculations by yourself, simply
%  comment the following line, loading the FNC values published in the
%  paper, and uncomment all the other lines in this cell.
%
load CC_bptf_published.mat



% % prepare a 3D matrix to store the FNC scores for each participant
% CC = zeros(nIC_fn,nIC_fn,nsubj);
% 
% 
% % NB: the subject number output by dualreg are from 0..358, therefore 
% % while the counter starts at 1, the first participant is (subj-1)
% for subj=1:nsubj
%     
%     
%     % load the dr_stage1 time courses for each subject, and demean
%     clear a    
%     a = dlmread([dr_dir '/dr_stage1_subject' sprintf('%05d', subj-1) '.txt']);
% 
%     meana = mean(a);
%     a_demean =  a - repmat(meana, size(a,1), 1);
%     a_demean_fn = a_demean(:,fn_relevant);
% 
%     % do the bandpass filtering between 0.009 and 0.08 Hz    
%     TR_subj_i = TRs(subj);
%     
%     HP = 1/(2 * 0.009 * TR_subj_i);
%     LP = 1/(2 * 0.08 * TR_subj_i);
%     
%     Nvox = size(a_demean_fn,2);
%     Ntimepoints = size(a_demean_fn,1);
%     
%     nii_tmp_bptf = reshape(a_demean_fn', [1 1 Nvox Ntimepoints]);
%     
%     save_avw(nii_tmp_bptf, 'tmp_nii.nii.gz','f', [1 1 1 3]);
% 
%     comando = ['fslmaths tmp_nii.nii.gz -bptf ' num2str(HP) ' ' num2str(LP) ' tmp_nii_bptf.nii.gz '];
%     system(comando)
%     
%     tmp = squeeze(read_avw('tmp_nii_bptf.nii.gz'));
%     
%     a_demean_fn_bptf = tmp';
%     
%     !rm tmp_nii.nii.gz
%     !rm tmp_nii_bptf.nii.gz
% 
%     % calculate FNC
%     CC(:,:,subj) = corrcoef(a_demean_fn_bptf);
%         
%     clc
%     disp(['loading dr_stage1 time courses... ' num2str(subj/nsubj*100) ' % done']);
% 
% 
% end
% 
% !rm -rf randomise
% 
% % remove tmp variables
% clear a a_demean a_demean_fn a_demean_fn_bptf





%% prepare the nifti image that will be used for randomise GLM

dim = size(CC,1);
nii_CC = reshape(CC, [dim dim 1 nsubj]);
save_avw(nii_CC,'nii_CC.nii.gz','f',[1 1 1 3]); 




%% define covariates of no interest

% Create a matrix made of one column of ones for each site of
% acquisition, which will serve to model the mean for that site.

site_indices = unique(SITE_ID);
howmanysites = length(unique(SITE_ID));

covariates_site = zeros(nsubj,howmanysites);

for ith_site = 1:howmanysites
    
    site_rows = find(SITE_ID==site_indices(ith_site));
    covariates_site(site_rows,ith_site) = 1;
    
end


% retrieve the covariates from phenomat
%  SEX is not used as a covariate since the participants in the 359 sample 
%  are all males.
%  Covariates controlling for site are added below
covariates_pheno = [AGE_AT_SCAN FIQ EYE_STATUS_AT_SCAN FD_mean];
covariates_pheno_demean = covariates_pheno - repmat(mean(covariates_pheno), [nsubj 1]);

%  site is modeled as one column of 1's (i.e. mean) for each site
% we add it here so that it does not get demeaned in the previous step
covariates = [covariates_pheno_demean covariates_site];

clear covariates_pheno covariates_pheno_demean




 

%% FSL randomise to estimate groups differences
%  Performed using the do_randomise_wrapper.m in the /addons directory


% Build a grouping variable
groups = zeros(nsubj,1);
groups(idx_ASD) = 1;
groups(idx_TD) = 2;


% Define the model matrix for randomise. 
model = zeros(nsubj, max(groups));
model(idx_ASD,1) = 1;
model(idx_TD,2) = 1;

nperms = 20000
disp(['Number of permutations = ' num2str(nperms)]);

% Define the contrasts and their names
contrasts = [1 -1    % ASD > TD
            -1  1];  % TD > ASD
         
contrast_names={'ASD > TD', 'TD > ASD'};


clc

reply = input(['Do you want to add covariates to the model? [Y/N]'], 's');

    if strcmp(reply, 'Y')

        % attach the covariates of no interest to the model matrix
        model = [model covariates];

        cov_contrasts = zeros(size(contrasts,1), size(covariates,2));    
        contrasts = [contrasts cov_contrasts];     
    end
    % end of attach the covariates

    
    

% run randomise
do_randomise_wrapper('nii_CC.nii.gz', model, contrasts, contrast_names, nperms)



p_ASD_maj_TD_whole = 1 - read_avw('randomise/ttest_vox_p_tstat1.nii.gz');
p_TD_maj_ASD_whole = 1 - read_avw('randomise/ttest_vox_p_tstat2.nii.gz');

t_ASD_maj_TD_whole = read_avw('randomise/ttest_tstat1.nii.gz');
t_TD_maj_ASD_whole = read_avw('randomise/ttest_tstat2.nii.gz');



% estimate FDR-corrected alpha value
p_vec_ASD_maj_TD_whole = [triu(p_ASD_maj_TD_whole,1)];
p_vec_ASD_maj_TD_whole = p_vec_ASD_maj_TD_whole(find(p_vec_ASD_maj_TD_whole));

p_vec_TD_maj_ASD_whole = [triu(p_TD_maj_ASD_whole,1)];
p_vec_TD_maj_ASD_whole = p_vec_TD_maj_ASD_whole(find(p_vec_TD_maj_ASD_whole));

p_vec = [p_vec_ASD_maj_TD_whole' p_vec_TD_maj_ASD_whole'];

pID_whole = FDR_nichols(p_vec, 0.05)


alpha_whole = pID_whole;
clc
disp(['Using alpha = ' num2str(alpha_whole)]);

save('FSL_RANDOMISE_20k_BPTF.mat')




%% Find indices of FDR-corrected significant differences

p_ASD_maj_TD_whole_sig = p_ASD_maj_TD_whole<=alpha_whole;
p_TD_maj_ASD_whole_sig = p_TD_maj_ASD_whole<=alpha_whole;


[isig1 jsig1] = find(triu(p_ASD_maj_TD_whole_sig,1));
[isig2 jsig2] = find(triu(p_TD_maj_ASD_whole_sig,1));


isig = [isig1 ; isig2];
jsig = [jsig1 ; jsig2];

clc
disp([' '])
disp(['There are ' num2str(length(isig)) ' group differences']);
disp([' '])




%% Plot the two matrices showing group differences

close all


figure

n_rows_cols = length(fn_relevant);

subplot(1,2,1)
imagesc(t_ASD_maj_TD_whole.*(p_ASD_maj_TD_whole <= alpha_whole));
title(['ASD > TD p(FDR)<' num2str(alpha_whole)])
axis square; colormap summer; colorbar
grid on


set(gca,'YTick',[1:n_rows_cols]);
set(gca,'XTick',[1:n_rows_cols]);
set(gca,'YTickLabel', labels);
set(gca, 'XTickLabel', labels);
xticklabel_rotate([],90)


subplot(1,2,2)
imagesc(t_TD_maj_ASD_whole.*(p_TD_maj_ASD_whole <= alpha_whole));
title(['TD > ASD p(FDR)<' num2str(alpha_whole)])
axis square; colormap summer; colorbar
grid on


set(gca,'YTick',[1:n_rows_cols]);
set(gca,'XTick',[1:n_rows_cols]);
set(gca,'YTickLabel', labels);
set(gca, 'XTickLabel', labels);
xticklabel_rotate([],90)



%% figure with RED = ASD>TD and BLUE = TD>ASD

figure

n_rows_cols = length(fn_relevant);

ASD_maj_TD = t_ASD_maj_TD_whole.*(p_ASD_maj_TD_whole <= alpha_whole);
TD_maj_ASD = t_TD_maj_ASD_whole.*(p_TD_maj_ASD_whole <= alpha_whole);

% the sign of TD_maj_ASD is inverted for purposes of visualization
results_matrix = ASD_maj_TD - TD_maj_ASD;
imagesc(results_matrix)

h = colormap(jet(256));
h(100:156,:)=1;
caxis([-5 5]);

axis square; colormap(h); colorbar
grid on

set(gca,'YTick',[1:n_rows_cols]);
set(gca,'XTick',[1:n_rows_cols]);
set(gca,'YTickLabel', labels);
set(gca, 'XTickLabel', labels);
xticklabel_rotate([],90)

title(['Presented results p(FDR)<' num2str(alpha_whole)]);





%% notBoxPlot with standard error of the mean and standard deviation
%  red = median
%  light red = standard error of the mean
%  blue = standard deviation


figure

% isig1 = ASD > TD
for i = 1:length(isig1)
    
   T_temp = num2str(t_ASD_maj_TD_whole(isig1(i), jsig1(i))); 
   p_temp = num2str(p_ASD_maj_TD_whole(isig1(i), jsig1(i)));
    
   subplot(3,3,i)

   y = squeeze(CC(isig1(i), jsig1(i), :));
   
   h=notBoxPlot(y, groups, 0.3);
   d = [h.data];
   set(d, 'markersize', 5, 'markerfacecolor', 'none', 'markeredgecolor', 'none');
   
   % boxplot(y, groups);
   axis([0.5 2.5  -1 1]);
   set(gca, 'XTick', [1 2]);
   set(gca, 'YTick', [-1 -0.5 0 0.5 1]);
   set(gca, 'XTickLabel', {'ASD','TD'})
   title(char(labels{isig1(i)} , labels{jsig1(i)}, ['T = ' T_temp ' p=' p_temp]));
   grid on
    
end


istart_2nd_plot = i;


% isig2 = TD > ASD
for i = 1:length(isig2)
    
   T_temp = num2str(t_TD_maj_ASD_whole(isig2(i), jsig2(i))); 
   p_temp = num2str(p_TD_maj_ASD_whole(isig2(i), jsig2(i)));
    
   subplot(3,3,i+istart_2nd_plot)

   y = squeeze(CC(isig2(i), jsig2(i), :));
   
   h=notBoxPlot(y, groups, 0.3);
   d = [h.data];
   set(d, 'markersize', 5, 'markerfacecolor', 'none', 'markeredgecolor', 'none');
   
   % boxplot(y, groups);
   axis([0.5 2.5  -1 1]);
   set(gca, 'XTick', [1 2]);
   set(gca, 'YTick', [-1 -0.5 0 0.5 1]);
   set(gca, 'XTickLabel', {'ASD','TD'})
   title(char(labels{isig2(i)} , labels{jsig2(i)}, ['T = ' T_temp ' p=' p_temp]));
   grid on
    
end




























