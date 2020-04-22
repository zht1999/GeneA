function [bestindividual,bestfit] = best(POP)
[px,py] = size(POP.y);
bestindividual = [POP.x1(1) , POP.x2(1)];
bestfit = POP.y(1);
for i = 2:px
    if POP.y(i)>bestfit
        bestindividual = [POP.x1(i) , POP.x2(i)];
        bestfit = POP.y(i);
    end
end
