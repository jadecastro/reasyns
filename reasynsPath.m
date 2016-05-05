
% disable JxBrowser to prevent idle cpu use
com.mathworks.mlwidgets.html.HtmlComponentFactory.setDefaultType('HTMLRENDERER');

addpath(genpath(pwd));
rmpath(genpath([pwd,'/drake/drake/examples/']));
rmpath([pwd,'/lib/ellipsoids/solvers/SeDuMi_1_1']);
rmpath([pwd,'/lib/mpt/solvers/SeDuMi_1_3']);

run([pwd,'/drake/drake/addpath_drake']);

examplesPath = [pwd,'/examples/'];
