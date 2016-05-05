
% disable JxBrowser to prevent idle cpu use
com.mathworks.mlwidgets.html.HtmlComponentFactory.setDefaultType('HTMLRENDERER');

addpath(genpath(pwd));
rmpath(genpath([pwd,'/drake']))  % removes the bundled (source) version of drake
rmpath([pwd,'/lib/ellipsoids/solvers/SeDuMi_1_1']);
rmpath([pwd,'/lib/mpt/solvers/SeDuMi_1_3']);
rmpath(genpath([pwd,'/drake-v0.9.11-win64/drake/examples/']));

% The following is preliminary, and only for testing
run([pwd,'/drake/drake/addpath_drake']);
%addpath(genpath('/home/jon/Software/mosek'));
%addpath(genpath('/home/jon/Software/drake-distro/spotless'));
%addpath('/drake/examples/DubinsCar');
%addpath('/home/jon/Software/drake-distro/drake');
%rmpath([pwd,'/drake-distro-lcm-win64/spotless']);
%rmpath([pwd,'/lib/spotless-master']);

