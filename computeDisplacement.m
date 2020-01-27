%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the displcement of the CoM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes the average displacement of the CoM during one
% entire trial. This can be useful to see that the forces are not so off
% that we'd expect the person to fall off the treadmill. It CANNOT be
% greater than the treadmill dimensions in the x-y plane. It also takes the
% product of this displacement with the average force measured from a trial
% to get an estiamte of the average net work on the center of mass.
% Ideally, this should be close to equal to the average net work performed
% on the center of mass as obtained using the individual limbs method.
% INPUTS:
% dFastStep: Displacement of CoM during fast step. It's the output from
%   computePower in the main code.
% dSlowStep: Displacement of CoM during slow step. It's the output from
%   computePower in the main code.
% fBias: The forces averaged over the entire duration of trial after
%   remvoing baseline drift/ offset. It's the output from forcePerStride in
%   the main code.
%
% OUTPUTS:
% avgdCoM: this is the average dispalcement of the center of mass over the
%   duration of the trial
% avgBiasWork: this is the avearge net work performed on the center of mass
%   if the average displacement of the center of mass is avgdCoM and an
%   avearge force of fBias acts on the center of mass during the trial.

function [avgdCoM,avgBiasWork] = computeDisplacement(dFastStep,dSlowStep,fBias)
dCoMx_all=[];dCoMy_all=[];dCoMz_all=[];
for m=1:length(dFastStep)
        dCoMx_all = [dCoMx_all;dFastStep{1,m}.com(:,1);dSlowStep{1,m}.com(:,1)];
        dCoMy_all = [dCoMy_all;dFastStep{1,m}.com(:,2);dSlowStep{1,m}.com(:,2)];
        dCoMz_all = [dCoMz_all;dFastStep{1,m}.com(:,3);dSlowStep{1,m}.com(:,3)];
        netdCoMx(m) = mean(dCoMx_all); % average displacement per stride
        netdCoMy(m) = mean(dCoMy_all);
        netdCoMz(m) = mean(dCoMz_all);
end
avgdCoM.x = mean(netdCoMx); % average displacement per trial
avgdCoM.y = mean(netdCoMy);
avgdCoM.z = mean(netdCoMz);

% this is the expected net work given the displacement of the CoM. Ideally,
% the displacement would be zero. Here we're just trying to see if the
% displacement is too much for some reason, such as maybe signal dependent
% noise in force.
avgBiasWork = sqrt((fBias.meanfxBiasStride*avgdCoM.x)^2 + (fBias.meanfyBiasStride*avgdCoM.y)^2 ...
    + (fBias.meanfzBiasStride*avgdCoM.z)^2);
end