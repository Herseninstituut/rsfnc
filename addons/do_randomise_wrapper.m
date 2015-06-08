function do_randomise_wrapper(data_filename, model, contrasts, contrast_names, nperms)
%
% This is a very simple and very INflexible wrapper for FSL randomise. 
%
% Creates the model for randomise and run it, storing everything in a
% directory `pwd`/randomise
%
% data_filename = 4D file with the matrix, format nIC x nIC x 1 x nsubj
% model = matrix nsubj x predictors
% contrast_names = cell array of strings containg the names of the
%                  contrasts
% nperms = number of permutations
%
% e.g.
% do_randomise('nii_CC_all.nii.gz', model, contrasts, contrast_names, 20000)
%
%
% Leonardo Cerliani, 2014


% creation of the model

system(['mkdir ' pwd '/randomise'])
system(['rm ' pwd '/randomise/*'])
system(['rm ' pwd '/randomise/design.mat'])
system(['rm ' pwd 'randomise/design.con'])


% create the mask, since randomise is not happy when the nifti image is
% made of only one voxel.
system([ 'fslroi  ' pwd '/' data_filename ' ' pwd '/randomise/mask.nii.gz 0 1'])
system([ 'fslmaths ' pwd '/randomise/mask.nii.gz -mul 0 -add 1 ' pwd '/randomise/mask.nii.gz']);


% design.mat
fid = fopen([pwd '/randomise/design.mat'],'wt');

fprintf(fid, '/NumWaves   %i \n', size(model,2));
fprintf(fid, '/NumPoints  %i \n', size(model,1));
fprintf(fid, '/PPheights  ');

for i = 1:size(model,2)

    fprintf(fid, '  %e  ', max(model(:,i)));

end

fprintf(fid, ' \n\n');
fprintf(fid, '/Matrix \n');




for i = 1:size(model,1)
   
   fprintf(fid, '%s \n', num2str(model(i,:))); 
    
end

fclose(fid)




% design.con
fid = fopen([pwd '/randomise/design.con'],'wt');

for i = 1:length(contrast_names)
   
    fprintf(fid, '/ContrastName%i   "%s" \n', i, contrast_names{i});

end

fprintf(fid, '/NumWaves   %i \n', size(model,2));
fprintf(fid, '/NumContrasts   %i \n', length(contrast_names));

fprintf(fid, '/PPheights  ');
for i = 1:length(contrast_names)
    fprintf(fid, '  %i  ', 1);
end
fprintf(fid, ' \n');


fprintf(fid, '/RequiredEffect  ');
for i = 1:length(contrast_names)
    fprintf(fid, '  %i  ', 1);
end
fprintf(fid, ' \n');

fprintf(fid, ' \n\n');
fprintf(fid, '/Matrix \n');


for i = 1:length(contrast_names)
   
   fprintf(fid, '%s \n', num2str(contrasts(i,:))); 
    
end

clc
system(['cat ' pwd '/randomise/design.con']);


%% run randomise

com = ['randomise -i ' pwd '/' data_filename ' -m ' pwd '/randomise/mask -o ' pwd '/randomise/ttest -d ' pwd '/randomise/design.mat -t ' pwd '/randomise/design.con -x -n ' num2str(nperms) ' -V']
disp(com)
system(com)

clc
system(['cat ' pwd '/randomise/design.con']);
ls([pwd '/randomise/*.nii.gz']);
disp([ num2str(nperms) 'permutations']);
sprintf('\n')
sprintf(com)









