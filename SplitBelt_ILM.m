% This code computes the mechanical work performed by a person walking on a
% split-belt treadmill, using the individual limbs method described in
% Donelan et al 2002.
% This code expects that each main folder contains a list of folders
% representing data from each participant. All participant folder names have
% a common key word to use to select them, and contain a mat file called
% Data. This code reads that mat file which is a structure containing
% 1. A structure called Fs which contains two fields- Force and Kinematics,
%   which contain the value of the frequency (Hz) at which force and kinematics
%   data were collected in the experiment.
% 2. A structure called Demographics which contains the Age, Sex,
%   Leg_length, and Mass as fields. Leg length is in mm and mass is in kg.
% 3. A structure called Trials whose fields are the various trial names
%   from the experiment, and each of them is a structure. Ex: for trials
%   are Standing, Walking baseline, Walking Split -15SLA etc.
% 
% References: ? 
% 1. Donelan, J. M. M., Kram, R., & Kuo, A. D. (2002). Simultaneous
% positive and negative external mechanical work in human walking. Journal
% of Biomechanics, 35(1), 117?124.
% https://doi.org/10.1016/S0021-9290(01)00169-5
% 2. ?Selgrade, B. P., Thajchayapong, M., Lee, G. E., Toney, M. E., &
% Chang, Y.-H. (2017). Changes in mechanical work during neural adaptation
% to asymmetric locomotion. The Journal of Experimental Biology, 220(16),
% 2993?3000. https://doi.org/10.1242/jeb.149450 
%
% This code was last updated by Surabhi Simha on Jan 25th 2020

clear
close all

fp = uigetdir(cd,'Select the Directory Where the Individual Data are Stored');
fns = dir(fullfile(fp,'S*')); % Here "2019" is the repeating term in every participant's folder name
folderInd = [fns.isdir]'; % This just gets a number for each file in the folder, to easily refernece the file
fns = fns(folderInd);

varIN.g = 9.81; % The constant for acceleration due to gravity

% We detect steps using ground reaction forces. In that algorithm, this is
% used as the required minimum duration of a step, to eliminate false
% detections
varIN.minStepDur = 0.4;

%%
% preallocation
M = zeros(length(fns),1); %mass
legLength = zeros(length(fns),1); %leg length
subj = cell(length(fns),1); %each cell contains results for one subject

for i=1:length(fns) %loop to go through each participant
    load(fullfile(fp,fns(i).name,'Data.mat'));
    varIN.Fs = Data.Fs.Force; % frequency at which force data were collected. Note that this assumes all participants were collected at same frequency
    M(i,1) = Data.Demographics.Mass; %participant's mass
    legLength(i,1) = str2num(Data.Demographics.Leg_Length)./1000; % converting to metres
    trials = Data.Trials;
    trialData = fieldnames(trials);

    % Analysis of all walking trials for one participant. All of the
    % analysed data for one participant is stored in their subj{i}
    % structure. 
    [subj{i}] = analyzeSubject(Data,varIN,trialData,M(i));
    
    display(strcat(num2str(i),' participants completed'))
end

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to analyse all trials of a single subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This combines all the processed data for a participant and stores them
% in two structures that are given as the output of this function.
% OUTPUTS:
% trial: This contains all the results for each trial. This is a cell array
% where each cell contains the analysed data of the corresponding trial.
function [trial] = analyzeSubject(Data,varIN,trialData,M)

for j=2:length(trialData) % this cycles through the differnet walking trials for each subject
    
    time.force = Data.Trials.(trialData{j}).Time_Vector_Force(1:end-1); %time vector for GRF
    
    speed.L = -1*Data.Trials.(trialData{j}).SpeedLeft; % multiply by -1 since the belts are moving in the
    speed.R = -1*Data.Trials.(trialData{j}).SpeedRight; % -y direction according to treadmill cordinate frame
    
    % storing the ground reaction forces into intuitive variable names
    grf.xl = Data.Trials.(trialData{j}).xGRF_L(1:end-1);
    grf.yl = Data.Trials.(trialData{j}).yGRF_L(1:end-1);
    grf.zl = Data.Trials.(trialData{j}).zGRF_L(1:end-1);
    grf.xr = Data.Trials.(trialData{j}).xGRF_R(1:end-1);
    grf.yr = Data.Trials.(trialData{j}).yGRF_R(1:end-1);
    grf.zr = Data.Trials.(trialData{j}).zGRF_R(1:end-1);
    
    [trial{j-1}] = analyzeTrial(varIN,time,speed,grf,M,trialData{j});
    
    display(strcat(trialData{j},' completed'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to analyse a single walking trial for a single subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% varOut: This is a structure whose fields are all the analysed data such
% as average work, work rate, average force etc.
function [varOut] = analyzeTrial(varIN,time,speed,grf,M,trialData)

t_force = time.force;
Fs = varIN.Fs;
g = varIN.g;
minStepNo = varIN.minStepDur;

% This breaks up the force vectors into strides, given by the left heel
% strike.
% fStride: the forces broken into strides
% tStride: time broken into strides
% rStepIndperStride: Index no. of the right heel strike in a given stride
% fBias: the forces averaged over the entire duration of trial after
%   remvoing baseline drift/ offset
[varOut.fStride,varOut.tStride,varOut.rStepIndperStride, varOut.fBias] =...
    forcePerStride(grf,g,minStepNo,t_force,M,trialData);

% This is the function that calculates velocity from force and then
% computes power. Here fast step refers to the step where the left foot
% (the leg on the fast belt) is in stance.
% dFastStep: displacement of CoM during fast step
% dSlowStep: displacement of CoM during slow step
% vStride: velocity of the CoM during one stride (begin and end with left
%   heels trike
% pStride: power generated due to each leg on the CoM, the belts,and total
%   power generated by each leg
[varOut.dFastStep, varOut.dSlowStep, varOut.vStride, varOut.pStride] =...
    computePower(varOut.rStepIndperStride,varOut.fStride,Fs,speed);

% This computes work from power 
% avgwork: average of the work performed throughout the trial 
% avg100workRate: average work rate of the last 100 strides. Work rate is
%   calculated as the work done in a stride divided by the duration of that
%   stride 
% wnetRate: this is a vector of net work rate done by the legs at each
%   stride in that trial 
% wnet: this is vector of net work done by the legs at each stride in that
%   trial
[varOut.avgwork,varOut.avg100workRate,varOut.wnetRate,varOut.wnet] =...
    workPerStride(varOut.pStride,varOut.tStride,Fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to break up the force vectors into strides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This breaks up the force vectors into strides, given by the left heel
% strike.
% OUTPUTS:
% fStride: the forces broken into strides
% tStride: time broken into strides
% rStepIndperStride: Index no. of the right heel strike in a given stride
% fBias: the forces averaged over the entire duration of trial after
%   remvoing baseline drift/ offset
function [fStride,tStride,rStrideEndInd,fBias] = forcePerStride(grf,g,minStepNo,t_force,M,trialData)

% this combines all left forces into one variable and all right
% forces into one.
flMatched = [grf.xl, grf.yl, grf.zl];
frMatched = [grf.xr, grf.yr, grf.zr];

% This part is in case there is some systematic noise in the forces
% measured i.e. there is a drift in the forces or there is a constant
% offset, independent of the forces on interest.
dfldt = diff(flMatched(:,3)); %derivative of force to see when its unchanging.
dfrdt = diff(frMatched(:,3));

dLZerodt = abs(dfldt)<1e-1; %testing when the derivative of force is zero
flZero = abs(flMatched(1:end-1,3))<0.05*max(flMatched(:,3)); %testing when the value of force isless than 5% body weight. to ensure that there is no foot on the belt.
flBaseRqd = dLZerodt.*flZero; %finding points in time when derivative of force is constant && force <5% body weight

% now we have a vector of the forces measured from the left force plate
% during the entire experiment when we expect the forces to be zero. using
% this, we calculate the slope and offset of the baseline force so that we
% can subtract it from the total forces.
xflBaseDrift = [ones(length(t_force(logical(flBaseRqd))),1),t_force(logical(flBaseRqd))]\flMatched(logical(flBaseRqd),1);
yflBaseDrift = [ones(length(t_force(logical(flBaseRqd))),1),t_force(logical(flBaseRqd))]\flMatched(logical(flBaseRqd),2);
zflBaseDrift = [ones(length(t_force(logical(flBaseRqd))),1),t_force(logical(flBaseRqd))]\flMatched(logical(flBaseRqd),3);

flMatched(:,1) = (flMatched(:,1)- (t_force*xflBaseDrift(2))) - xflBaseDrift(1);
flMatched(:,2) = (flMatched(:,2)- (t_force*yflBaseDrift(2))) - yflBaseDrift(1);
flMatched(:,3) = (flMatched(:,3)- (t_force*zflBaseDrift(2))) - zflBaseDrift(1);

% we repeat the same process for the right belt
dRZerodt = abs(dfrdt)<1e-1;
frZero = abs(frMatched(1:end-1,3))<0.05*max(frMatched(:,3));
frBaseRqd = dRZerodt.*frZero;

% calculate the slope and offset of the baseline force
xfrBaseDrift = [ones(length(t_force(logical(frBaseRqd))),1),t_force(logical(frBaseRqd))]\frMatched(logical(frBaseRqd),1);
yfrBaseDrift = [ones(length(t_force(logical(frBaseRqd))),1),t_force(logical(frBaseRqd))]\frMatched(logical(frBaseRqd),2);
zfrBaseDrift = [ones(length(t_force(logical(frBaseRqd))),1),t_force(logical(frBaseRqd))]\frMatched(logical(frBaseRqd),3);

frMatched(:,1) = (frMatched(:,1)- (t_force*xfrBaseDrift(2))) - xfrBaseDrift(1);
frMatched(:,2) = (frMatched(:,2)- (t_force*yfrBaseDrift(2))) - yfrBaseDrift(1);
frMatched(:,3) = (frMatched(:,3)- (t_force*zfrBaseDrift(2))) - zfrBaseDrift(1);

% Choose the end of the stride by finding the time when the
% vertical force crosses 32N(Young-Hui Chang 2017) in the positive
% direction i.e. right after heel strike. Ignore steps detected
% that are less than 400ms long.
lstepEndInd = (flMatched(:,3)-32)<0.5;
lstepEndInd = lstepEndInd*1000;
changeInd = find(diff(lstepEndInd)~=0);
for p=1:length(changeInd)-1
    if (changeInd(p+1)-changeInd(p))<minStepNo
        lstepEndInd(changeInd(p):changeInd(p+1))=1000;
    end
end
lStrideEndInd = find(diff(lstepEndInd)<0);
lStrideEndInd = lStrideEndInd+1;

% Repeat the same to find the right heel strike. This is used to obtain the
% center of mass displacement
for i=2:length(lStrideEndInd)
    rStepEndInd = (frMatched(lStrideEndInd(i-1):lStrideEndInd(i),3)-32)<0.5;
    if ~rStepEndInd, rStrideEndInd(i-1)=nan; continue; end
    rStepEndInd = rStepEndInd*1000;
    rChangeInd = find(diff(rStepEndInd)~=0);
    for p=1:length(rChangeInd)-1
        if (rChangeInd(p+1)-rChangeInd(p))<minStepNo
            rStepEndInd(rChangeInd(p):rChangeInd(p+1))=1000;
        end
    end
    rStrideEndInd(i-1) = find(diff(rStepEndInd)<0,1);
    rStrideEndInd(i-1) = rStrideEndInd(i-1)+1;
    if ~rStrideEndInd(i-1), rStrideEndInd(i-1)=nan; end
end

tStride = diff(t_force(lStrideEndInd));

% break up force data of left and right leg into strides
fxBiasStride=zeros(length(lStrideEndInd)-1,1);
fyBiasStride=zeros(length(lStrideEndInd)-1,1);
fzBiasStride=zeros(length(lStrideEndInd)-1,1);
flyBiasStride=zeros(length(lStrideEndInd)-1,1);
fryBiasStride=zeros(length(lStrideEndInd)-1,1);

for m=2:length(lStrideEndInd)
    fStride.l{m-1,:} = flMatched(lStrideEndInd(m-1):lStrideEndInd(m),:);
    fStride.r{m-1,:} = frMatched(lStrideEndInd(m-1):lStrideEndInd(m),:);
    fStride.total{m-1,:} = [fStride.l{m-1,1}(:,1:3)+fStride.r{m-1,1}(:,1:3)];
    
    % i can sum up the total force over each stride. ideally this
    % should be zero for each stride for the x and y directions.
    % this is to test if there is something causing non-zero forces
    % during the trial.
    fxBiasStride(m-1) = mean(fStride.total{m-1,1}(:,1));
    fyBiasStride(m-1) = mean(fStride.total{m-1,1}(:,2));
    fzBiasStride(m-1) = mean(fStride.total{m-1,1}(:,3));
    flyBiasStride(m-1) = mean(fStride.l{m-1,1}(:,2));
    fryBiasStride(m-1) = mean(fStride.r{m-1,1}(:,2));
end

fBias.meanfxBiasStride = mean(fxBiasStride);
fBias.meanfyBiasStride = mean(fyBiasStride);
fBias.meanfzBiasStride = mean(fzBiasStride)-(M*g);

if abs(fBias.meanfxBiasStride)>5 || abs(fBias.meanfyBiasStride)>5 || abs(fBias.meanfzBiasStride)>10
    warning(strcat('Average forces over the trial ',trialData,...
        ' are larger than expected. Use the function computeDisplacement to verify the errors are within acceptable ranges.'))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that calculates velocity from force and then
% computes power. Here fast step refers to the step where the left foot
% (the leg on the fast belt) is in stance.
% OUTPUTS:
% dFastStep: displacement of CoM during fast step
% dSlowStep: displacement of CoM during slow step
% vStride: velocity of the CoM during one stride (begin and end with left
%   heels trike
% pStride: power generated due to each leg on the CoM, the belts,and total
%   power generated by each leg
function [dFastStep, dSlowStep, vStride, pStride]...
    = computePower(rStepIndperStride,fStride,Fs,speed)

ftot = fStride.total;
fl = fStride.l;
fr = fStride.r;
speedL = speed.L;
speedR = speed.R;
g = 9.81;

% Integrates total forces over a stride to get center of mass velocity.
for i=1:length(ftot)
    % estimate mass for each stride. this is important since we are
    % intergrating force and then dividing by mass to get velocity. Ideally,
    % we would have very precise and accurate estimates of force and mass
    % so that we can just divide the integral by the same exact constant
    % in each stride. But since we do not, we estimate mass at each stride
    % to account for the noise in force from that stride, and to get a
    % continuous looking velocity. i.e. the velocity at the end of one
    % stride matches the velocity at the beginning of the next stride.
    mass = mean(ftot{i}(:,3))/g;

    % velocity of the center of mass.
    % medio-lateral velocity
    vcom{i}(:,1) = (cumtrapz(ftot{i}(:,1)) / (mass*Fs));
    vcom{i}(:,1) = vcom{i}(:,1) - mean(vcom{i}(:,1));
    % fore-aft velcoity
    vcom{i}(:,2) = (cumtrapz(ftot{i}(:,2)) / (mass*Fs));
    vcom{i}(:,2) = vcom{i}(:,2) - mean(vcom{i}(:,2));
    % vertical velocity
    vcom{i}(:,3) = (cumtrapz(ftot{i}(:,3)-(mass*g)) / (mass*Fs));
    vcom{i}(:,3) = vcom{i}(:,3) - mean(vcom{i}(:,3));
    
    % CoM displacement
    dcom{i}(:,1) = cumtrapz(vcom{i}(:,1))/Fs; %medio-lateral
    dcom{i}(:,2) = cumtrapz(vcom{i}(:,2))/Fs; %fore-aft
    dcom{i}(:,3) = cumtrapz(vcom{i}(:,3))/Fs; %vertical
    
    % velocity of the left belt
    vltread{i}(:,1) = zeros(length(ftot{i}(:,1)),1); %medio-lateral
    vltread{i}(:,2) = ones(length(ftot{i}(:,1)),1) * speedL; %fore-aft
    vltread{i}(:,3) = zeros(length(ftot{i}(:,1)),1); %vertical
    
    % velocity of the right belt
    vrtread{i}(:,1) = zeros(length(ftot{i}(:,1)),1); %medio-lateral
    vrtread{i}(:,2) = ones(length(ftot{i}(:,1)),1) * speedR; %fore-aft
    vrtread{i}(:,3) = zeros(length(ftot{i}(:,1)),1); %vertical
    
    % compute mechanical power
    plcom{i} = dot(fl{i}(:,1:3),vcom{i}(:,1:3),2); %power generated on the CoM due to the left leg
    prcom{i} = dot(fr{i}(:,1:3),vcom{i}(:,1:3),2); %power generated on the CoM due to the right leg
    pltread{i} = dot(fl{i}(:,1:3),vltread{i}(:,1:3),2); %power generated by left leg on left belt
    prtread{i} = dot(fr{i}(:,1:3),vrtread{i}(:,1:3),2); %power generated by right leg on right belt
    pl{i} = plcom{i} + pltread{i}; %total power generated by left leg
    pr{i} = prcom{i} + prtread{i}; %total power generated by the right leg
    
end

%separate displacement of com into that during the fast step and slow step
n=1;
for m=1:length(dcom)
    if isnan(rStepIndperStride(m))|| rStepIndperStride(m)==0, continue; end
    dFastStep{n}.com = dcom{m}(1:rStepIndperStride(m),1:3);
    dSlowStep{n}.com = dcom{m}(rStepIndperStride(m):end,1:3);
    n=n+1;
end

% combine all related varaibles
vStride.com = vcom;
vStride.ltread = vltread;
vStride.rtread = vrtread;
pStride.lcom = plcom;
pStride.rcom = prcom;
pStride.ltread = pltread;
pStride.rtread = prtread;
pStride.l = pl;
pStride.r = pr;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes work from power 
% OUTPUTS:
% avgwork: average of the work performed throughout the trial 
% avg100workRate: average work rate of the last 100 strides. Work rate is
%   calculated as the work done in a stride divided by the duration of that
%   stride 
% wnetRate: this is a vector of net work rate done by the legs at each
%   stride in that trial 
% wnet: this is vector of net work done by the legs at each stride in that
%   trial
function [avgwork,avg100workRate,wnetRate,wnet] = workPerStride(pStride,tStride,Fs)
   
%% Compute work done per stride
    for q=1:length(pStride.lcom)
        plcomstride = pStride.lcom{q};
        prcomstride = pStride.rcom{q};
        pltreadstride = pStride.ltread{q};
        prtreadstride = pStride.rtread{q};
        plstride = pStride.l{q};
        prstride = pStride.r{q};

        %% work and work rate: performed on the center of mass
        %work performed on the center of mass through the left leg. note that
        %this does not necessarily mean that the work was done by the leg, some
        %of it is done by the treadmill but tranferred through the leg. so,
        %this will not correspond to the energy spent by the leg to do work on
        %the center of mass.
        wlcompos(q) = trapz(plcomstride(plcomstride>0))*(1/Fs);
        wlcomneg(q) = trapz(plcomstride(plcomstride<0))*(1/Fs);

        %work performed on the center of mass through the right leg.
        wrcompos(q) = trapz(prcomstride(prcomstride>0))*(1/Fs);
        wrcomneg(q) = trapz(prcomstride(prcomstride<0))*(1/Fs);  
        
        %work rate performed on the center of mass through the right leg.
        wlcomposRate(q) = wlcompos(q)/tStride(q);
        wlcomnegRate(q) = wlcomneg(q)/tStride(q);
        wrcomposRate(q) = wrcompos(q)/tStride(q);
        wrcomnegRate(q) = wrcomneg(q)/tStride(q);
        
        %work performed on the center of mass through both legs
        wcompos(q) = wlcompos(q)+wrcompos(q);
        wcomneg(q) = wlcomneg(q)+wrcomneg(q);
        wnetcom(q) = wcompos(q)+wcomneg(q);

        wcomposRate(q) = wcompos(q)/tStride(q);
        wcomnegRate(q) = wcomneg(q)/tStride(q);
        wnetcomRate(q) = wnetcom(q)/tStride(q);
        %% work and work rate: performed on the belts by the legs
        %work performed on the left belt by the left leg.
        wltreadpos(q) = trapz(pltreadstride(pltreadstride>0))*(1/Fs);
        wltreadneg(q) = trapz(pltreadstride(pltreadstride<0))*(1/Fs);

        %work performed on the right belt by the right leg.
        wrtreadpos(q) = trapz(prtreadstride(prtreadstride>0))*(1/Fs);
        wrtreadneg(q) = trapz(prtreadstride(prtreadstride<0))*(1/Fs);

        %work rate performed by the legs on the respective belts
        wltreadposRate(q) = wltreadpos(q)/tStride(q);
        wltreadnegRate(q) = wltreadneg(q)/tStride(q);
        wrtreadposRate(q) = wrtreadpos(q)/tStride(q);
        wrtreadnegRate(q) = wrtreadneg(q)/tStride(q);

        %total work performed on the belts by the legs
        wtreadpos(q) = wltreadpos(q)+wrtreadpos(q);
        wtreadneg(q) = wltreadneg(q)+wrtreadneg(q);
        wnettread(q) = wtreadpos(q)+wtreadneg(q);

        wtreadposRate(q) = wtreadpos(q)/tStride(q);
        wtreadnegRate(q) = wtreadneg(q)/tStride(q);
        wnettreadRate(q) = wnettread(q)/tStride(q);
        %% total work and work rate: performed by the legs
        %total work performed by the left leg
        wlpos(q) = trapz(plstride(plstride>0))*(1/Fs);
        wlneg(q) = trapz(plstride(plstride<0))*(1/Fs);

        %total work performed by the right leg
        wrpos(q) = trapz(prstride(prstride>0))*(1/Fs);
        wrneg(q) = trapz(prstride(prstride<0))*(1/Fs);

        %total positive work performed by both legs
        wtotpos(q) = wlpos(q)+wrpos(q);
        wtotneg(q) = wlneg(q)+wrneg(q);
        wtotposRate(q) = wtotpos(q)/tStride(q);
        wtotnegRate(q) = wtotneg(q)/tStride(q);

        %net work performed by both legs
        wnet(q) = wtotpos(q) + wtotneg(q);
        wnetRate(q) = wnet(q)/tStride(q);

        %total work rate by individual legs
        wlposRate(q) = wlpos(q)/tStride(q);
        wlnegRate(q) = wlneg(q)/tStride(q);
        wrposRate(q) = wrpos(q)/tStride(q);
        wrnegRate(q) = wrneg(q)/tStride(q);
    end

%% average work and work rate: done across all strides.
avgwork.net = mean(wnet);
avgwork.totpos = mean(wtotpos);
avgwork.totneg = mean(wtotneg);
avgwork.lpos = mean(wlpos);
avgwork.lneg = mean(wlneg);
avgwork.rpos = mean(wrpos);
avgwork.rneg = mean(wrneg);

avgwork.nettread = mean(wnettread);
avgwork.treadpos = mean(wtreadpos);
avgwork.treadneg = mean(wtreadneg);
avgwork.ltreadpos = mean(wltreadpos);
avgwork.ltreadneg = mean(wltreadneg);
avgwork.rtreadpos = mean(wrtreadpos);
avgwork.rtreadneg = mean(wrtreadneg);

avgwork.netcom = mean(wnetcom);
avgwork.compos = mean(wcompos);
avgwork.comneg = mean(wcomneg);
avgwork.lcompos = mean(wlcompos);
avgwork.lcomneg = mean(wlcomneg);
avgwork.rcompos = mean(wrcompos);
avgwork.rcomneg = mean(wrcomneg);

avgworkRate.net = mean(wnetRate);
avgworkRate.totpos = mean(wtotposRate);
avgworkRate.totneg = mean(wtotnegRate);
avgworkRate.lpos = mean(wlposRate);
avgworkRate.lneg = mean(wlnegRate);
avgworkRate.rpos = mean(wrposRate);
avgworkRate.rneg = mean(wrnegRate);

avgworkRate.nettread = mean(wnettreadRate);
avgworkRate.treadpos = mean(wtreadposRate);
avgworkRate.treadneg = mean(wtreadnegRate);
avgworkRate.ltreadpos = mean(wltreadposRate);
avgworkRate.ltreadneg = mean(wltreadnegRate);
avgworkRate.rtreadpos = mean(wrtreadposRate);
avgworkRate.rtreadneg = mean(wrtreadnegRate);

avgworkRate.netcom = mean(wnetcomRate);
avgworkRate.compos = mean(wcomposRate);
avgworkRate.comneg = mean(wcomnegRate);
avgworkRate.lcompos = mean(wlcomposRate);
avgworkRate.lcomneg = mean(wlcomnegRate);
avgworkRate.rcompos = mean(wrcomposRate);
avgworkRate.rcomneg = mean(wrcomnegRate);

%% average work rate: done across only the last 100 strides of that trial i.e. steady state
avg100workRate.net = mean(wnetRate(end-100:end));
avg100workRate.totpos = mean(wtotposRate(end-100:end));
avg100workRate.totneg = mean(wtotnegRate(end-100:end));
avg100workRate.lpos = mean(wlposRate(end-100:end));
avg100workRate.lneg = mean(wlnegRate(end-100:end));
avg100workRate.rpos = mean(wrposRate(end-100:end));
avg100workRate.rneg = mean(wrnegRate(end-100:end));

avg100workRate.nettread = mean(wnettreadRate(end-100:end));
avg100workRate.treadpos = mean(wtreadposRate(end-100:end));
avg100workRate.treadneg = mean(wtreadnegRate(end-100:end));
avg100workRate.ltreadpos = mean(wltreadposRate(end-100:end));
avg100workRate.ltreadneg = mean(wltreadnegRate(end-100:end));
avg100workRate.rtreadpos = mean(wrtreadposRate(end-100:end));
avg100workRate.rtreadneg = mean(wrtreadnegRate(end-100:end));

avg100workRate.netcom = mean(wnetcomRate(end-100:end));
avg100workRate.compos = mean(wcomposRate(end-100:end));
avg100workRate.comneg = mean(wcomnegRate(end-100:end));
avg100workRate.lcompos = mean(wlcomposRate(end-100:end));
avg100workRate.lcomneg = mean(wlcomnegRate(end-100:end));
avg100workRate.rcompos = mean(wrcomposRate(end-100:end));
avg100workRate.rcomneg = mean(wrcomnegRate(end-100:end));

end