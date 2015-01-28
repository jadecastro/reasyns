
% Automaton definition
clear aut
aut.q{1} = 1;
aut.q{2} = 2;
aut.q{3} = 3;
aut.q{4} = 2;
aut.trans{1} = [1 2];
aut.trans{2} = [2 3];
aut.trans{3} = [2 1];
aut.trans{4} = [3 4];
aut.trans{5} = [4 1];
%aut.trans{6} = [1 1];

aut.isReactive(1) = false;
aut.isReactive(2) = true;
aut.isReactive(3) = true;
aut.isReactive(4) = false;
aut.isReactive(5) = false;