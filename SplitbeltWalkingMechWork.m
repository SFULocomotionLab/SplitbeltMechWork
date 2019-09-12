% This code computes the mechanical work performed by a person walking on a
% split-belt treadmill, using the individual limbs method decribed in
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
% This code was last updated by Surabhi Simha on July 13th 2019

clear
close all

fp = uigetdir(cd,'Select the Directory Where the Individual Data are Stored');
fns = dir(fullfile(fp,'2019*')); % Here "2019" is the repeating term in every participant's folder name
folderInd = [fns.isdir]'; % This just gets a number for each file in the folder, to easily refernece the file
fns = fns(folderInd);

varIN.g = 9.81; % The constant for acceleration due to gravity

% We detect steps using ground reaction forces. In that algorithm, this is
% used as the required minimum duration of a step, to eliminate false
% detections
varIN.minStepDur = 0.4;

%% Some tests to help debug/ verify that the analyses are correct.
% Select "1" for the following depending on the test you want to run.

% This interpolates every stride to be equal length and computes the
% average force, velocity and power (FVP) across the last 100 stride. Use this to
% see the average FVP plots for the last 100 strides. This may not help as
% much with debugging since it is an average but is useful to visualize
% and maybe present to others to explain FVP plots. Since this computes the
% average of the last 100 strides for each participant, the "plot" option
% below produces a figure for each participant. It pauses for the user to
% tap a key after producing each figure.
test.fvpInterp = 0; % This interpolates every stride to be equal length and computes the average force, velocity and power across the last 100 strides.
test.fvpInterpPlot = 0; % This plots the interpolated data. If this is 0, no figures are generated. test.fvpInterp has to be 1 for this to work.

% This generates an FVP plot where the x-axis is actual experiment time.
% This is very useful for debugging. It gives the time series data for each
% participant, one by one. Using this, we can check that: 
% 1. The forces look correct. For ex: there should not be any drift in the
%   forces over time, the vertical force should not go below zero, and that
%   the general profile looks right i.e. there is double peak for walking.
%   Note that the force is raw data collected from the force plates. So,
%   errors here are indicative of errors in data collection or filtering of
%   forces. So, this is a logical first step in debugging.
% 2. The velocity looks correct. Main thing to check is that the velocity
%   at the end of a stride is equal to the velcocity at the start of the
%   next stride. Also, compare the profile and magnitude to literature.
% 3. The power looks correct. It should be the dot product of the force and
%   velocity for every time step. Compare profile with literature.
test.fvpTimeSeries = 0;

% This outputs the force per stride averaged across all strides. This is
% also a very useful tool for debugging. The average force in the x and y
% axes should be zero, and in the z-axis it should be equal to body weight.
% These are the ideal conditions. But I've never seen this to be exactly
% true. In my experience, upto 5N away from the ideal is still fine. When
% comparing the vertical force with body weight, make sure that the manner
% in which body weight itself was measured is reliable and not without
% errors/ noise/ drift etc. It is possible that this trial is fine but the
% recorded body weight is wrong.
test.forceperstride = 0;

% This is similarly useful as checking the total force. This displays the
% average displacement of the center of mass over a trial, which would
% ideally be zero and CANNOT be greater than the treadmill dimensions in
% the x-y plane. It also takes the product of this displacement with the
% average force from the above test to get an estiamte of the average net
% work on the center of mass. Ideally, this should be close to equal to the
% average net work performed on the center of mass as obtained using the
% individual limbs method.
test.checkBias=0;
%%
% Select this to save data for further analysis. It will be saved in a mat
% file named ProcessedData<date>.m
% This happens in the section starting at line 164.
savedata = 0;


% preallocation
M = zeros(length(fns),1); %mass
legLength = zeros(length(fns),1); %leg length
standMet = struct; %results from the standing trial
subj = cell(length(fns),1); %each cell contains results for one subject
rawSLA = cell(length(fns),1); %Step Length Asymmetry measured using kinematics (This is calculated beforehand, not done here)

for i=1:length(fns) %loop to go through each participant
    load(fullfile(fp,fns(i).name,'Data.mat'));
    varIN.Fs = Data.Fs.Force; % frequency at which force data were collected. Note that this assumes all participants were collected at same frequency
    M(i,1) = Data.Demographics.Mass; %participant's mass
    legLength(i,1) = str2num(Data.Demographics.Leg_Length)./1000; % converting to metres
    trials = Data.Trials;
    trialData = fieldnames(trials);
    % this assigns a trial order number based on the trial names. These are
    % the trials we had in our second project.
    for k = 1:length(trialData)
        switch(trialData{k})
            case('StandingBaseline'); trialOrder(i,k) = -1;
            case('Base'); trialOrder(i,k) = 1;
            case('Split_1'); trialOrder(i,k) = 2;
            case('Split_2'); trialOrder(i,k) = 3;
            case('Split_3'); trialOrder(i,k) = 4;
            case('Split_4'); trialOrder(i,k) = 5;
            case('Split_5'); trialOrder(i,k) = 6;
            case('PostTied'); trialOrder(i,k) = 7;
        end
    end
    % Standing Trial  - to obtain mass and compare; and obtain metabolic
    % power for standing.
    zfl_stand = Data.Trials.(trialData{1}).zGRF_L;
    zfr_stand = Data.Trials.(trialData{1}).zGRF_R;
    mass = mean(zfl_stand+zfr_stand)/varIN.g;
    VO2 = Data.Trials.StandingBaseline.AcumVO2_L;
    VCO2 = Data.Trials.StandingBaseline.AcumVCO2_L;
    t_met = Data.Trials.StandingBaseline.Time_Vector_MetCost;
    [standMet(i).RER,standMet(i).Emet,standMet(i).Pmet] = computeMetPower(t_met, VO2, VCO2, 2);
    
    % Analysis of all walking trials for one participant. All of the
    % analysed data for one participant is stored in their subj{i}
    % structure. 
    [subj{i},rawSLA{i}] = analyzeSubject(Data,varIN,trialData,test,M(i));
    
    for q=1:length(subj{i})
        % Interpolates to average FVP of last 100 strides and plots FVP. Refer
        % to tests above. Note that these values are not saved in the subj
        % structure. So, these need to be separately added the list of
        % variables to be saved, if required.
        % avgforce: time average of the forces from the last 100 strides
        % avgvel: time average of the velocities from the last 100 strides
        % avgpower: time average of the power from the last 100 strides
        if test.fvpInterp
            [avgforce(i,q),avgvel(i,q),avgpower(i,q)] = computeStrideAvg(i,subj{i}{q}.fStride,...
                subj{i}{q}.vStride,subj{i}{q}.pStride,varIN.Fs,test,fns(i).name,trialData{q+1});
        end
        % plots FVP time series for each participant
        if test.fvpTimeSeries
            fvpTimeSeries(subj{i}{q}.fStride,subj{i}{q}.vStride,subj{i}{q}.pStride,trialData{q+1});
        end
    end
    display(strcat(num2str(i),' participants completed'))
    if test.checkBias==1, pause; end
end
%% This creates a 2D array of only the variables that we want to save.
% Add variables you want to save under the second for loop.
if savedata
    clear p q
    RER = zeros(length(subj),1); Emet = zeros(length(subj),1); Pmet = zeros(length(subj),1);
    % avg100workRate = struct;
    for p=1:length(subj) % participant loop
        for q=1:length(subj{p,1}) % trial loop
            RER(p,q) = subj{p,1}{1,q}.RER;
            Emet(p,q) = subj{p,1}{1,q}.Emet;
            Pmet(p,q) = subj{p,1}{1,q}.Pmet;
            avg100workRate(p,q) = subj{p,1}{1,q}.avg100workRate;
            wnetRate{p,q} = subj{p,1}{1,q}.wnetRate;
        end
    end
    dateInfo = datetime();
    sYear=num2str(dTemp.Year);
    sMonth = num2str(dTemp.Month);
    sDay=num2str(dTemp.Day);
    sep ='-';
    save(['ProcessedData',sYear,sep,sMonth,sep,sDay],'RER','Emet','Pmet',...
        'trialOrder','avg100workRate','wnetRate');
end

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to analyse all trials of a single subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This combines all the processed data for a participant and stores them
% in two structures that are given as the output of this function.
% OUTPUTS:
% trial: This contains all the results of analysing forces (ex: mechanical
% work) and metabolic cost for each trial. This is a cell array where each
% cell contains the analysed data of the corresponding trial.
% rawSLA: This is a structure that contains the Step Length Asymmetry data
% as provided  in the original Data structure (obtained using kinematics).
% It just now organized in a manner similar to how the force data is
% analyzed.
function [trial,rawSLA] = analyzeSubject(Data,varIN,trialData,test,M)

for j=2:length(trialData) % this cycles through the differnet walking trials for each subject
    
    %These are contained in the Data structure
    SLA = Data.Trials.(trialData{j}).SLA(1:end-1)'; %measured step length asymmetry. It is obtained through kinematics.
    rawSLA(j,1:length(SLA)) = SLA;
    VO2 = Data.Trials.(trialData{j}).AcumVO2_L(1:end-1); %measured volume of oxygen in liters
    VCO2 = Data.Trials.(trialData{j}).AcumVCO2_L(1:end-1); %measured volume of CO2 in litres
    
    time.force = Data.Trials.(trialData{j}).Time_Vector_Force(1:end-1); %time vector for GRF
    time.lHS_kinematics = Data.Trials.(trialData{j}).HS_Time_Left(1:end-1)'; %time of each left heel strike based on video capture
    time.met = Data.Trials.(trialData{j}).Time_Vector_MetCost(1:end-1); %time vector for metabolics
    
    speed.L = -1*Data.Trials.(trialData{j}).SpeedLeft; % multiply by -1 since the belts are moving in the
    speed.R = -1*Data.Trials.(trialData{j}).SpeedRight; % -y direction according to treadmill cordinate frame
    
    % storing the ground reaction forces into intuitive variable names
    grf.xl = Data.Trials.(trialData{j}).xGRF_L(1:end-1);
    grf.yl = Data.Trials.(trialData{j}).yGRF_L(1:end-1);
    grf.zl = Data.Trials.(trialData{j}).zGRF_L(1:end-1);
    grf.xr = Data.Trials.(trialData{j}).xGRF_R(1:end-1);
    grf.yr = Data.Trials.(trialData{j}).yGRF_R(1:end-1);
    grf.zr = Data.Trials.(trialData{j}).zGRF_R(1:end-1);
    
    cop.L = Data.Trials.(trialData{j}).Filtered_COP_L_FA;
    cop.R = Data.Trials.(trialData{j}).Filtered_COP_R_FA;
    
    [trial{j-1}] = analyzeTrial(varIN,VO2,VCO2,time,speed,grf,test,M,trialData{j});
    
    if test.checkBias==1
        figure(88);hold on;
        plot(j,trial{j-1}.avgwork.netcom,'ro');
        plot(j,trial{j-1}.avgBiasWork,'bo')
        display(strcat('Average medio-lateral displacement of CoM during',...
            num2str(trialData{j}),' is', num2str(round(trial{j-1}.avgdCoM.x,3)),' metres'));
        display(strcat('Average fore-aft displacement of CoM during',...
            num2str(trialData{j}),' is', num2str(round(trial{j-1}.avgdCoM.y,3)),' metres'));
        display(strcat('Average vertical displacement of CoM during',...
            num2str(trialData{j}),' is', num2str(round(trial{j-1}.avgdCoM.z,3)),' metres'));
    end
    display(strcat(trialData{j},' completed'))
end

%This just labels the figure if the test is selected
if test.checkBias==1
    figure(88);
    title({'Comparing avg net work on CoM:','Product of avg force and avg displacement of CoM (blue)','Individual limbs method (red)'});
    ylabel('Work done on Center of Mass (J)')
    xlabel('Trial Number')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to analyse a single walking trial for a single subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% varOut: This is a structure whose fields are all the analysed data such
% as average work, work rate, average force etc.
function [varOut] = analyzeTrial(varIN,VO2,VCO2,time,speed,grf,test,M,trialData)

t_force = time.force;

t_met = time.met;
Fs = varIN.Fs;
g = varIN.g;
minStepNo = varIN.minStepDur;

% This gives the metabolic energy and power for this trial. The output
% values are 
% RER: average Respiratory Exchange Ratio from the last 3min of each trial.
% Emet: average energy consumed during the last 3min of the trial (Joule)
% Pmet: average metabolic power from the last 3min (Watt)
[varOut.RER,varOut.Emet,varOut.Pmet] = computeMetPower(t_met, VO2, VCO2);

% This breaks up the force vectors into strides, given by the left heel
% strike.
% rStepIndperStride: Index no. of the right heel strike in a given stride
% fStride: the forces broken into strides
% tStride: time broken into strides
[varOut.fStride,varOut.tStride,varOut.rStepIndperStride, varOut.fBias] =...
    forcePerStride(grf,test,g,minStepNo,t_force,M,trialData);

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

% This computes the expected net displacement of the CoM during the entire
% trial. This can be useful to see that the forces are not so off that we'd
% expect the person to fall off the treadmill
% avgdCoM: this is the average dispalcement of the center of mass over the
%   duration of the trial
% avgBiasWork: this is the avearge net work performed on the center of mass
%   if the average displacement of the center of mass is avgdCoM and an
%   avearge force of fBias acts on the center of mass during the trial.
[varOut.avgdCoM,varOut.avgBiasWork] = ...
    computeDisplacement(varOut.dFastStep,varOut.dSlowStep,varOut.fBias);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to break up the force vectors into strides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This breaks up the force vectors into strides, given by the left heel
% strike.
% OUTPUTS:
% rStepIndperStride: Index no. of the right heel strike in a given stride
% fStride: the forces broken into strides
% tStride: time broken into strides
% fBias: the forces averaged over the entire duration of trial after
%   remvoing baseline drift/ offset
function [fStride,tStride,rStrideEndInd,fBias] = forcePerStride(grf,test,g,minStepNo,t_force,M,trialData)

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

if abs(fBias.meanfxBiasStride)>5 || abs(fBias.meanfyBiasStride)>5 || abs(fBias.meanfzBiasStride)>5
    warning(strcat('Average forces over the trial ',trialData,' are larger than expected. Use the forceperstride and CheckBias tests to verify the errors are within acceptable ranges.'))
end

if test.forceperstride==1
    % i can sum up the total force over each stride. here i just
    % average this force across all strides. ideally this should be
    % zero for each stride for the x and y directions. this is to test
    % if there is something causing non-zero forces during the trial.
%     figure(1);hold on;plot(fxBiasStride);pause
%     figure(2);hold on;plot(fyBiasStride);grid on;pause
    fBias.meanflyBiasStride = mean(flyBiasStride);
%     figure(3);hold on;plot(flyBiasStride);pause
    fBias.meanfryBiasStride = mean(fryBiasStride);
%     figure(4);hold on;plot(fryBiasStride);pause
    display(strcat('Forces averaged over the entire duration of trial: ',trialData))
    disp(fBias)
    pause
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

        %work performed on the center of mass through both leg
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
avgwork.lcompos = mean(wlcompos);
avgwork.lcomneg = mean(wlcomneg);

avgwork.rcompos = mean(wrcompos);
avgwork.rcomneg = mean(wrcomneg);

avgwork.ltreadpos = mean(wltreadpos);
avgwork.ltreadneg = mean(wltreadneg);

avgwork.rtreadpos = mean(wrtreadpos);
avgwork.rtreadneg = mean(wrtreadneg);

avgwork.lpos = mean(wlpos);
avgwork.lneg = mean(wlneg);

avgwork.rpos = mean(wrpos);
avgwork.rneg = mean(wrneg);

avgwork.totpos = mean(wtotpos);
avgwork.totneg = mean(wtotneg);

avgwork.net = mean(wnet);
avgwork.nettread = mean(wnettread);
avgwork.netcom = mean(wnetcom);
avgwork.treadpos = mean(wtreadpos);
avgwork.compos = mean(wcompos);
avgwork.treadneg = mean(wtreadneg);
avgwork.comneg = mean(wcomneg);

avgworkRate.net = mean(wnetRate);
avgworkRate.totpos = mean(wtotposRate);
avgworkRate.totneg = mean(wtotnegRate);
avgworkRate.nettread = mean(wnettreadRate);
avgworkRate.netcom = mean(wnetcomRate);
avgworkRate.treadpos = mean(wtreadposRate);
avgworkRate.compos = mean(wcomposRate);
avgworkRate.treadneg = mean(wtreadnegRate);
avgworkRate.comneg = mean(wcomnegRate);
avgworkRate.lpos = mean(wlposRate);
avgworkRate.lneg = mean(wlnegRate);
avgworkRate.rpos = mean(wrposRate);
avgworkRate.rneg = mean(wrnegRate);
avgworkRate.ltreadpos = mean(wltreadposRate);
avgworkRate.ltreadneg = mean(wltreadnegRate);
avgworkRate.rtreadpos = mean(wrtreadposRate);
avgworkRate.rtreadneg = mean(wrtreadnegRate);

%% average work and work rate: done across only the last 100 strides of that trial i.e. steady state SLA
avg100workRate.net = mean(wnetRate(end-100:end));
avg100workRate.totpos = mean(wtotposRate(end-100:end));
avg100workRate.totneg = mean(wtotnegRate(end-100:end));
avg100workRate.nettread = mean(wnettreadRate(end-100:end));
avg100workRate.netcom = mean(wnetcomRate(end-100:end));
avg100workRate.treadpos = mean(wtreadposRate(end-100:end));
avg100workRate.compos = mean(wcomposRate(end-100:end));
avg100workRate.treadneg = mean(wtreadnegRate(end-100:end));
avg100workRate.comneg = mean(wcomnegRate(end-100:end));
avg100workRate.lpos = mean(wlposRate(end-100:end));
avg100workRate.lneg = mean(wlnegRate(end-100:end));
avg100workRate.rpos = mean(wrposRate(end-100:end));
avg100workRate.rneg = mean(wrnegRate(end-100:end));
avg100workRate.ltreadpos = mean(wltreadposRate(end-100:end));
avg100workRate.ltreadneg = mean(wltreadnegRate(end-100:end));
avg100workRate.rtreadpos = mean(wrtreadposRate(end-100:end));
avg100workRate.rtreadneg = mean(wrtreadnegRate(end-100:end));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the metabolic power for each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This gives the metabolic energy and power for this trial.
% OUTPUTS: 
% RER: average Respiratory Exchange Ratio from the last 3min of each trial.
% Emet: average energy consumed during the last 3min of the trial (Joule)
% Pmet: average metabolic power from the last 3min (Watt)
function [RER,Emet, Pmet] = computeMetPower(t_met, VO2, VCO2,Dur)
if ~exist('Dur','var'), Dur = 3; end
[~,metAvgBegin] = min(abs(t_met-(t_met(end)-Dur)));
VO2tot = VO2(end)-VO2(metAvgBegin); %Total volume of oxygen consumed in the duration of our interest
VCO2tot = VCO2(end)-VCO2(metAvgBegin); %total volume of carbon dioxide produced in the duration of our interest 
Emet = (VO2tot*1000*16.48) + (VCO2tot*1000*4.48); % applying the Brockway equation
Tavg = ((t_met(end)-t_met(metAvgBegin))*60); % Measure the actual duration over which this energy is measured
Pmet = Emet/Tavg; % Divide by the duration to get power
RER = VCO2tot/VO2tot;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the displcement of the CoM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes the expected net displacement of the CoM during the entire
% trial. This can be useful to see that the forces are not so off that we'd
% expect the person to fall off the treadmill
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to interpolate and compute average (across last 100 strides) force, velocity and power within a stride
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% avgforce: time average of the forces from the last 100 strides
% avgvel: time average of the velocities from the last 100 strides
% avgpower: time average of the power from the last 100 strides
function [avgforce,avgvel,avgpower] = computeStrideAvg(m,fStride,vStride,pStride,Fs,test,fn,trialData)

fl = fStride.l;
fr = fStride.r;
vcom = vStride.com;
vltread = vStride.ltread;
vrtread = vStride.rtread;
plcom = pStride.lcom;
prcom = pStride.rcom;
pltread = pStride.ltread;
prtread = pStride.rtread;
pl = pStride.l;
pr = pStride.r;

% Preallocate
flxinterp=zeros(1000,length(plcom)); frxinterp=zeros(1000,length(plcom));
flyinterp=zeros(1000,length(plcom)); fryinterp=zeros(1000,length(plcom));
flzinterp=zeros(1000,length(plcom)); frzinterp=zeros(1000,length(plcom));
vxcominterp=zeros(1000,length(plcom)); vlxtreadinterp=zeros(1000,length(plcom));
vrxtreadinterp=zeros(1000,length(plcom)); vycominterp=zeros(1000,length(plcom));
vlytreadinterp=zeros(1000,length(plcom)); vrytreadinterp=zeros(1000,length(plcom));
vzcominterp=zeros(1000,length(plcom)); vlztreadinterp=zeros(1000,length(plcom));
vrztreadinterp=zeros(1000,length(plcom)); plcominterp=zeros(1000,length(plcom));
prcominterp=zeros(1000,length(plcom)); pltreadinterp=zeros(1000,length(plcom));
prtreadinterp=zeros(1000,length(plcom)); plinterp=zeros(1000,length(plcom));
printerp=zeros(1000,length(plcom)); tvals=zeros(1000,length(plcom));


for i=1:length(plcom)
    % This stores the power over a stride from a cell array to a double to
    % make calculations easier.
    plcomstride = plcom{i};
    prcomstride = prcom{i};
    pltreadstride = pltread{i};
    prtreadstride = prtread{i};
    plstride = pl{i};
    prstride = pr{i};
    
    flxstride = fl{i}(:,1);
    frxstride = fr{i}(:,1);
    flystride = fl{i}(:,2);
    frystride = fr{i}(:,2);
    flzstride = fl{i}(:,3);
    frzstride = fr{i}(:,3);
    
    vxcomstride = vcom{i}(:,1);
    vlxtreadstride = vltread{i}(:,1);
    vrxtreadstride = vrtread{i}(:,1);
    vycomstride = vcom{i}(:,2);
    vlytreadstride = vltread{i}(:,2);
    vrytreadstride = vrtread{i}(:,2);
    vzcomstride = vcom{i}(:,3);
    vlztreadstride = vltread{i}(:,3);
    vrztreadstride = vrtread{i}(:,3);
    
    strideLen = length(frzstride)-1;
    tvals(:,i) = linspace(0,strideLen,1000);
    
    % This is interpolating across a stride to contain 1000 points.
    flxinterp(:,i) = interp1(flxstride,tvals(:,i));
    frxinterp(:,i) = interp1(frxstride,tvals(:,i));
    flyinterp(:,i) = interp1(flystride,tvals(:,i));
    fryinterp(:,i) = interp1(frystride,tvals(:,i));
    flzinterp(:,i) = interp1(flzstride,tvals(:,i));
    frzinterp(:,i) = interp1(frzstride,tvals(:,i));
    
    vxcominterp(:,i) = interp1(vxcomstride,tvals(:,i));
    vlxtreadinterp(:,i) = interp1(vlxtreadstride,tvals(:,i));
    vrxtreadinterp(:,i) = interp1(vrxtreadstride,tvals(:,i));
    vycominterp(:,i) = interp1(vycomstride,tvals(:,i));
    vlytreadinterp(:,i) = interp1(vlytreadstride,tvals(:,i));
    vrytreadinterp(:,i) = interp1(vrytreadstride,tvals(:,i));
    vzcominterp(:,i) = interp1(vzcomstride,tvals(:,i));
    vlztreadinterp(:,i) = interp1(vlztreadstride,tvals(:,i));
    vrztreadinterp(:,i) = interp1(vrztreadstride,tvals(:,i));
    
    plcominterp(:,i) = interp1(plcomstride,tvals(:,i));
    prcominterp(:,i) = interp1(prcomstride,tvals(:,i));
    pltreadinterp(:,i) = interp1(pltreadstride,tvals(:,i));
    prtreadinterp(:,i) = interp1(prtreadstride,tvals(:,i));
    plinterp(:,i) = interp1(plstride,tvals(:,i));
    printerp(:,i) = interp1(prstride,tvals(:,i));
end

%% Calculating average values of power over a stride

avgforce.lx = mean(flxinterp(:,end-100:end),2);
avgforce.rx = mean(frxinterp(:,end-100:end),2);
avgforce.ly = mean(flyinterp(:,end-100:end),2);
avgforce.ry = mean(fryinterp(:,end-100:end),2);
avgforce.lz = mean(flzinterp(:,end-100:end),2);
avgforce.rz = mean(frzinterp(:,end-100:end),2);

avgvel.xcom = mean(vxcominterp(:,end-100:end),2);
avgvel.lxtread = mean(vlxtreadinterp(:,end-100:end),2);
avgvel.rxtread = mean(vrxtreadinterp(:,end-100:end),2);
avgvel.ycom = mean(vycominterp(:,end-100:end),2);
avgvel.lytread = mean(vlytreadinterp(:,end-100:end),2);
avgvel.rytread = mean(vrytreadinterp(:,end-100:end),2);
avgvel.zcom = mean(vzcominterp(:,end-100:end),2);
avgvel.lztread = mean(vlztreadinterp(:,end-100:end),2);
avgvel.rztread = mean(vrztreadinterp(:,end-100:end),2);

avgpower.lcom = mean(plcominterp(:,end-100:end),2);
avgpower.rcom = mean(prcominterp(:,end-100:end),2);
avgpower.ltread = mean(pltreadinterp(:,end-100:end),2);
avgpower.rtread = mean(prtreadinterp(:,end-100:end),2);
avgpower.l = mean(plinterp(:,end-100:end),2);
avgpower.r = mean(printerp(:,end-100:end),2);

if test.fvpInterpPlot==1
    % Plot the positive, negative, and net power from slow and fast belts
    xvals = linspace(0,1,1000);
    figure(m); hold on;
    subplot(3,2,1); hold on; grid on; title(strcat(trialData,' Left belt'))
    plot(xvals,avgforce.lx,'r-');
    plot(xvals,avgforce.ly,'g-');
    plot(xvals,avgforce.lz,'b-');
    ylabel('Force (N)'); ylim([-200 900])
    legend('x','y','z')
    
    subplot(3,2,2); hold on; grid on; title(strcat(trialData,' Right belt'))
    plot(xvals,avgforce.rx,'r-');
    plot(xvals,avgforce.ry,'g-');
    plot(xvals,avgforce.rz,'b-');
    ylabel('Force (N)'); ylim([-200 900])
    legend('x','y','z')
    
    subplot(3,2,3); hold on; grid on;
    plot(xvals,avgvel.xcom,'r-');
    plot(xvals,avgvel.ycom,'g-');
    plot(xvals,avgvel.zcom,'b-');
    plot(xvals,avgvel.lytread,'g--');
    ylabel('Velocity (m/s)'); ylim([-1.5 0.5])
    legend('x','y','z')
    
    subplot(3,2,4); hold on; grid on;
    plot(xvals,avgvel.xcom,'r-');
    plot(xvals,avgvel.ycom,'g-');
    plot(xvals,avgvel.zcom,'b-');
    plot(xvals,avgvel.rytread,'g--');
    ylabel('Velocity (m/s)'); ylim([-1.5 0.5])
    legend('x','y','z')
    
    subplot(3,2,5); hold on; grid on;
    plot(xvals,avgpower.l,'r-');
%     plot(xvals,avgpower.lcom,'g-');
%     plot(xvals,avgpower.ltread,'b-');
    ylabel('Power (W)'); ylim([-300 300])
    xlabel('Fraction of Stride')
    
    subplot(3,2,6); hold on; grid on;
    plot(xvals,avgpower.r,'r-')
%     plot(xvals,avgpower.rcom,'g-')
%     plot(xvals,avgpower.rtread,'b-')
    ylabel('Power (W)'); ylim([-300 300]); 
    xlabel('Fraction of Stride'); pause
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the fvp time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fvpTimeSeries(fStride,vStride,pStride,trialData)
% This can be used to trouble shoot by looking at individual FVP
% plots for each subject. Choose fvpTimeSeries = 1 at the beginning if
% you want to look at this.

% Combining all stride by stride values into a vector for easy
% plotting.
fl_all=[]; fr_all=[];vcom_all=[];vltread_all=[];vrtread_all=[];plcom_all=[];
prcom_all=[];pltread_all=[];prtread_all=[];pl_all=[];pr_all=[];
meanPowerLPerStride=[];meanPowerRPerStride=[];
for n=1:length(pStride.l)
    fl_all = [fl_all;fStride.l{n}];
    fr_all = [fr_all;fStride.r{n}];
    vcom_all = [vcom_all;vStride.com{n}];
    vltread_all = [vltread_all;vStride.ltread{n}];
    vrtread_all = [vrtread_all;vStride.rtread{n}];
    plcom_all = [plcom_all;pStride.lcom{n}];
    prcom_all = [prcom_all;pStride.rcom{n}];
    pltread_all = [pltread_all;pStride.ltread{n}];
    prtread_all = [prtread_all;pStride.rtread{n}];
    pl_all = [pl_all;pStride.l{n}];
    pr_all = [pr_all;pStride.r{n}];
    meanPowerLPerStride = [meanPowerLPerStride;mean((pStride.l{n}))];
    meanPowerRPerStride = [meanPowerRPerStride;mean((pStride.r{n}))];
end
ftot_all = fl_all+fr_all;
plcom_all = [plcom_all];
prcom_all = [prcom_all];
pltread_all = [pltread_all];
prtread_all = [prtread_all];
pl_all = [pl_all];
pr_all = [pr_all];

figure(1); clf; hold on;
subplot(3,2,1); hold on; grid on; ylim([-200 900]); grid on; 
title(strcat(trialData,' Left belt'))
plot(fl_all(:,1),'r-')
plot(fl_all(:,2),'g-')
plot(fl_all(:,3),'b-')
ylabel('Force (N)');
legend('x','y','z')

subplot(3,2,2); hold on; grid on; ylim([-200 900]); grid on;
title(strcat(trialData,' Right belt'))
plot(fr_all(:,1),'r-')
plot(fr_all(:,2),'g-')
plot(fr_all(:,3),'b-')
ylabel('Force (N)');

subplot(3,2,3); hold on; grid on; ylim([-1.5 0.5]); grid on; 
plot(vcom_all(:,1),'r-')
plot(vcom_all(:,2),'g-')
plot(vcom_all(:,3),'b-')
plot(vltread_all(:,2),'g--')
ylabel('Velocity (m/s)')

subplot(3,2,4); hold on; grid on; ylim([-1.5 0.5]); grid on; 
plot(vcom_all(:,1),'r-')
plot(vcom_all(:,2),'g-')
plot(vcom_all(:,3),'b-')
plot(vrtread_all(:,2),'g--')
ylabel('Velocity (m/s)')

subplot(3,2,5); hold on; grid on; ylim([-300 300]); grid on; 
plot(pl_all(:,1),'r-');
ylabel('Velocity (m/s)')
ylabel('Power')
xlabel('Sample time');

subplot(3,2,6); hold on; grid on; ylim([-300 300]); grid on; 
plot(pr_all(:,1),'r-');
ylabel('Power')
xlabel('Sample time');pause

end