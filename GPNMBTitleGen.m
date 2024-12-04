function [titleAtPlay] = GPNMBTitleGen(pAtPlay, pArray, numReps, numPerms )


if pAtPlay >= .05/numReps || isnan(pAtPlay)
titleAtPlay='q > .05';
elseif pAtPlay < .05/numReps && pAtPlay >0
    titleAtPlay=['q=', num2str(  sum(pAtPlay > pArray)/numPerms  )];
elseif pAtPlay==0 %this can be misleading in perm tests
titleAtPlay= ['q< ', num2str(1/numPerms)];
end
