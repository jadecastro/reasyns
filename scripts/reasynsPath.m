
% disable JxBrowser to prevent idle cpu use
com.mathworks.mlwidgets.html.HtmlComponentFactory.setDefaultType('HTMLRENDERER');

addpath(genpath([pwd,'/../../reasyns']));
rmpath([pwd,'/../../reasyns/lib/ellipsoids/solvers/SeDuMi_1_1']);
rmpath([pwd,'/../../reasyns/lib/mpt/solvers/SeDuMi_1_3']);
rmpath([pwd,'/../../reasyns/lib/mpt/solvers/SeDuMi_1_3_old']);
rmpath([pwd,'/../../reasyns/lib/drake-distro-lcm-win64/drake/examples/']);
addpath(genpath([pwd,'/../../reasyns/lib/drake-distro-lcm-win64/drake/examples/DubinsCar']));

% The following is preliminary, and only for testing
rmpath(genpath([pwd,'/../../reasyns/lib/drake-distro-lcm-win64']));
run('/home/jon/Software/drake-distro/drake/addpath_drake');
addpath(genpath('/home/jon/Software/mosek'));
addpath(genpath('/home/jon/Software/drake-distro/spotless'));
addpath('/home/jon/Software/drake-distro/drake/examples/DubinsCar');
addpath('/home/jon/Software/drake-distro/drake');
rmpath([pwd,'/../../reasyns/lib/drake-distro-lcm-win64/spotless']);
rmpath([pwd,'/../../reasyns/lib/spotless-master']);

