% % set the cellfie directory
% CELLFIEDIR = fileparts(mfilename('fullpath'));
% 
% % initialize the COBRA Toolbox
% if exist('CBTDIR', 'var')
%     cd(CBTDIR);
%     initCobraToolbox
% else
%     error('Please initialize the COBRA Toolbox or set the main directory of the COBRA Toolbox CBTDIR');
% end
% 
% % change to the cellfie main directory
% cd(CELLFIEDIR);
% 
% % add the entire cellfie repository to the MATLAB path
% addpath(genpath(CELLFIEDIR));
