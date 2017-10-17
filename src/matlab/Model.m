clc
clear all
tic

%% PARAMETERS

% Parameters that we can change
aggregationKernelConstant=1e-6; %1e-7
breakageKernelConstant = 0;%0.8; %0.868
breakageProbability=0;
consolidationConstant=0;%1e-5
initialPorosity=0.5;
minimumPorosity=0.2;
granuleSaturationFactor=0; %0.2035*1e-3
numberOfCompartments=3;
timeStep=0.5; % seconds
% particleAverageVelocity=0; % m/s
% liquidAdditionRate=0; % m3/s
sieveGrid=[38, 63, 90, 125, 250, 355, 500, 710, 850, 1000, 1400, 2000, 2380, 4000]; %Sieve (in micron)
numberOfSieveGrid=numel(sieveGrid); % Grid size

% Geometry
granulatorLength=0.38; % meter
impellerDiameter=0.114; % meter

% Material
solidDensity=476; % kg/m3
liquidDensity=1000; % kg/m3

% Process
premixingTime=5; % seconds
liquidAdditionTime =25;  % seconds
throughput=15; % kg/hr
liquidToSolidRatio=0.35;
particleResidenceTime=20; % seconds %11.07
impellerSpeed=1000; % RPM

% PBM Bins
numberOfFirstSolidBins=16; %Number of first solid bins
s_coef=5d-16; % m^3
s_base=3; %Volume increament between two bins

numberOfSecondSolidBins=16; %Number of second solid bins
ss_coef=5d-16; % m^3
ss_base=3; %Volume increament between twobins

% DEM
collisionEfficiencyConstant=0.01;
timeStepDEM=1e-7; %%%%%%WILL COME FROM DEM%%%%%%
numberOfDEMBins=16; %%%%%%WILL COME FROM DEM%%%%%%

% Particle Inflow--------------------------
fIn(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
fIn(1,1)=726657587;
fIn(2,2)=286654401;
fIn(3,3)=118218011;
fIn(4,4)=50319795;
fIn(5,5)=20954036;
fIn(6,6)=7345998;
fIn(7,7)=1500147;
fIn(8,8)=76518;
fIn(9,9)=149;
fIn(10,10)=0;
fIn(11,11)=0;
fIn(12,12)=0;
fIn(13,13)=0;
fIn(14,14)=0;
fIn(15,15)=0;
% Particle Inflow__________________

%% INITIAL CALCULATION
% Process Parameter Calculation----------
finalTime=premixingTime+liquidAdditionTime; % seconds
distanceBetweenCompartments= granulatorLength/numberOfCompartments; % m
particleAverageVelocity=granulatorLength/particleResidenceTime; % m/s
liquidAdditionRateAllCompartments(1:numberOfCompartments)=0;
liquidAdditionRate=((liquidToSolidRatio*throughput)/liquidDensity)/3600; % m3/s
liquidAdditionRateAllCompartments(1)=liquidAdditionRate/numberOfCompartments;
% Process Parameter Calculation_______

% PBM Bin-------------------------------------
vs(1:numberOfFirstSolidBins)=0;
vss(1:numberOfSecondSolidBins)=0;
diameter(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
for s=1:numberOfFirstSolidBins
    vs(s)=s_coef*s_base^(s-1); % m^3
end
for ss=1:numberOfSecondSolidBins
    vss(ss)=ss_coef*ss_base^(ss-1); % m^3
end
for s=1:numberOfFirstSolidBins
    for ss=1:numberOfSecondSolidBins
        diameter(s,ss)=((6/pi)*(vs(s)+vss(ss)))^(1/3);
    end
end
% PBM Bin______________________

% DEM
DEMDiameter(1:numberOfDEMBins)=1e3*rand(1,numberOfDEMBins); % micron %%%%%%WILL COME FROM DEM%%%%%%
numberOfCollisions(1:numberOfDEMBins,1:numberOfDEMBins)=1; %%%%%%WILL COME FROM DEM%%%%%%
numberOfImpacts(1:numberOfDEMBins)=1; %%%%%%WILL COME FROM DEM%%%%%%
% numberOfCollisions(1:numberOfDEMBins,1:numberOfDEMBins)=1*rand(numberOfDEMBins); %%%%%%WILL COME FROM DEM%%%%%%
% numberOfImpacts(1:numberOfDEMBins)=1*rand(1,numberOfDEMBins); %%%%%%WILL COME FROM DEM%%%%%%

% Cell Average----------------------------
[s_meshxy, ss_meshxy] =ndgrid(vs,vss);

s_agg=bsxfun(@plus,vs,vs');
ss_agg=bsxfun(@plus,vss,vss');
sAggregationCheck=repmat(reshape(s_agg<=vs(numberOfFirstSolidBins),[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]),[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]);
ssAggregationCheck=repmat(reshape(ss_agg<=vss(numberOfSecondSolidBins),[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]),[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]);

[s_mesh, ss_mesh]=ndgrid(vs,vss);
s_low(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=s_mesh(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins);
s_high(1:numberOfFirstSolidBins-1,1:numberOfSecondSolidBins)=s_mesh(2:numberOfFirstSolidBins,1:numberOfSecondSolidBins);
s_high(numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
ss_low(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=ss_mesh(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins);
ss_high(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins-1)=ss_mesh(1:numberOfFirstSolidBins,2:numberOfSecondSolidBins);
ss_high(1:numberOfFirstSolidBins,numberOfSecondSolidBins)=0;

sloc(1:numberOfFirstSolidBins,1:numberOfFirstSolidBins)=floor(log(s_agg/s_coef)/log(s_base)+1);
ssloc(1:numberOfSecondSolidBins,1:numberOfSecondSolidBins)=floor(log(ss_agg/ss_coef)/log(ss_base)+1);
sind=repmat(reshape(sloc,[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]),[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]);
ssind=repmat(reshape(ssloc,[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]),[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]);

s_break=bsxfun(@minus,vs,vs');
s_break(s_break<0)=0;
ss_break=bsxfun(@minus,vss,vss');
ss_break(ss_break<0)=0;

sloc_break(1:numberOfFirstSolidBins,1:numberOfFirstSolidBins)=(floor(log(s_break/s_coef)/log(s_base)+1))';
ssloc_break(1:numberOfSecondSolidBins,1:numberOfSecondSolidBins)=(floor(log(ss_break/ss_coef)/log(ss_base)+1))';
s_check_b=repmat(reshape(sloc_break>=1,[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]),[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]);
ss_check_b=repmat(reshape(ssloc_break>=1,[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]),[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]);

s_mesh_break=repmat(reshape(vs,[numberOfFirstSolidBins 1 1 1]),[1 numberOfSecondSolidBins numberOfFirstSolidBins numberOfSecondSolidBins])-repmat(reshape(vs,[1 1 numberOfFirstSolidBins 1]),[numberOfFirstSolidBins numberOfSecondSolidBins 1 numberOfSecondSolidBins]);
ss_mesh_break=repmat(reshape(vss,[1 numberOfSecondSolidBins 1 1]),[numberOfFirstSolidBins 1 numberOfFirstSolidBins numberOfSecondSolidBins])-repmat(reshape(vss,[1 1 1 numberOfSecondSolidBins]),[numberOfFirstSolidBins numberOfSecondSolidBins numberOfFirstSolidBins 1]);

sind_b=repmat(reshape(sloc_break,[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]),[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]);
ssind_b=repmat(reshape(ssloc_break,[1 numberOfSecondSolidBins 1 numberOfSecondSolidBins]),[numberOfFirstSolidBins 1 numberOfFirstSolidBins 1]);
sind_b(sind_b<1)=numberOfFirstSolidBins+1;
ssind_b(ssind_b<1)=numberOfSecondSolidBins+1;
% Cell Average___________________

%% INITIALIZATION
fAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
flAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
fgAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
dfdtAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
dfldtAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
dfgdtAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
liquidBinsAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
internalVolumeBinsAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
externalVolumeBinsAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
gasBinsAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
totalVolumeBinsAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
internalLiquidAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
externalLiquidAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
aggregationKernelAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
collisionEfficiencyAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
collisionFrequencyAllCompartments(1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
time=0; %Process Time initialization
timeIndex=1;
% load fInitial.mat

%% TIME LOOP
while time<=finalTime
    Time(timeIndex)=time;
    
    for c=1:numberOfCompartments
        f=reshape(fAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        fl=reshape(flAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        fg=reshape(fgAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        liquidAdditionRate=liquidAdditionRateAllCompartments(c);
        
        if c==1
            fPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
            flPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
            fgPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
            fComingIn=fIn;
            fgComingIn=reshape((fIn.*(repmat(vs,[numberOfSecondSolidBins 1])+repmat(vss,[numberOfFirstSolidBins 1]))),[numberOfFirstSolidBins numberOfSecondSolidBins])*initialPorosity*timeStep;
        else
            fPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=fAllCompartments(c-1,:,:);
            flPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=flAllCompartments(c-1,:,:);
            fgPreviousCompartment(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=fgAllCompartments(c-1,:,:);
            fComingIn(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
            fgComingIn(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        end
        
        [dfdt,dfldt,dfgdt,liquidBins,gasBins,internalVolumeBins,aggregationKernel,collisionEfficiency,collisionFrequency,breakageKernel]=Compartment(fComingIn,fgComingIn,fPreviousCompartment,flPreviousCompartment,fgPreviousCompartment,particleAverageVelocity,distanceBetweenCompartments,premixingTime,liquidAdditionTime,liquidAdditionRate,fl,fg,f,time,timeStep, aggregationKernelConstant,numberOfCollisions,breakageKernelConstant,consolidationConstant,minimumPorosity,granuleSaturationFactor,numberOfFirstSolidBins,numberOfSecondSolidBins,vs,vss,s_meshxy,ss_meshxy,sAggregationCheck,ssAggregationCheck,sind,ssind,sind_b,ssind_b,s_high,s_low,ss_high,ss_low,s_check_b,ss_check_b,s_mesh_break,ss_mesh_break,numberOfDEMBins,timeStepDEM,collisionEfficiencyConstant,DEMDiameter,diameter,numberOfImpacts,breakageProbability);
        
        dfdtAllCompartments(c,:,:)=dfdt;
        dfldtAllCompartments(c,:,:)=dfldt;
        dfgdtAllCompartments(c,:,:)=dfgdt;
        liquidBinsAllCompartments(c,:,:)=liquidBins;
        gasBinsAllCompartments(c,:,:)=gasBins;
        internalVolumeBinsAllCompartments(c,:,:)=internalVolumeBins;
        aggregationKernelAllCompartments(c,:,:,:,:)=aggregationKernel;
        collisionEfficiencyAllCompartments(c,:,:,:,:)=collisionEfficiency;
        collisionFrequencyAllCompartments(c,:,:,:,:)=collisionFrequency;
        
    end
    
    int_check=[-dfdtAllCompartments(fAllCompartments(:)~=0)./fAllCompartments(fAllCompartments(:)~=0); -dfldtAllCompartments(flAllCompartments(:)~=0)./flAllCompartments(flAllCompartments(:)~=0); -dfgdtAllCompartments(fgAllCompartments(:)~=0)./fgAllCompartments(fgAllCompartments(:)~=0)];
    
    maxAll = -Inf;
    maxLiquid = -Inf;
    maxGas = -Inf;
    if(sum(sum(sum(fAllCompartments~=0))) > 0)
        maxAll = max(-dfdtAllCompartments(fAllCompartments~=0)./fAllCompartments(fAllCompartments(:)~=0))
    end
    if(sum(sum(sum(flAllCompartments~=0))) > 0)
        maxLiquid = max(-dfldtAllCompartments(flAllCompartments(:)~=0)./flAllCompartments(flAllCompartments(:)~=0))
    end
    if(sum(sum(sum(fgAllCompartments~=0))) > 0)
        maxGas = max(-dfgdtAllCompartments(fgAllCompartments(:)~=0)./fgAllCompartments(fgAllCompartments(:)~=0))
    end
    
    maxofthree = max(maxAll, max(maxLiquid, maxGas))
    
    while and(max(max(max(int_check(:))))<0.1/timeStep,timeStep<0.25);
        timeStep=timeStep*2;
    end
    while and(max(max(max(int_check(:))))>=0.1/timeStep,timeStep>=5e-5);
        timeStep=timeStep/2;
    end
    
    fAllCompartments=fAllCompartments+dfdtAllCompartments*timeStep; %Calculation of PSD for next time
    flAllCompartments=flAllCompartments+dfldtAllCompartments*timeStep; %Calculation of liquid PSD for next time
    fgAllCompartments=fgAllCompartments+dfgdtAllCompartments*timeStep; %Calculation of gas PSD for next time
    
    flAllCompartments(flAllCompartments<0)=0;
    fgAllCompartments(fgAllCompartments<0)=0;
    
    if sum(sum(sum(isnan(fAllCompartments))))>0
        keyboard
    end
    
    if min(min(min(fAllCompartments)))<0
        keyboard
    end
    
    % Bin Calculation-----------------------
    for c=1:numberOfCompartments
        f=reshape(fAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        fl=reshape(flAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        fg=reshape(fgAllCompartments(c,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
        
        for s=1:numberOfFirstSolidBins
            for ss=1:numberOfSecondSolidBins
                if f(s,ss)~=0 % If there are particles in the bin...
                    liquidBins(s,ss)=fl(s,ss)/f(s,ss); % New liquid bins are calculated as total amount liquid in that size class divided by the number of particle in that size class
                    gasBins(s,ss)=fg(s,ss)/f(s,ss); % New gas bins are calculated as total amount liquid in that size class divided by the number of particle in that size class
                else
                    liquidBins(s,ss)=0; % No particles in the bin -> particles in that size class has no liquid
                    gasBins(s,ss)=0; % No particles in the bin -> particles in that size class has no gas
                end
            end
        end
        
        internalLiquid(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        externalLiquid(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        internalVolumeBins(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        externalVolumeBins(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        totalVolumeBins(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        for s=1:numberOfFirstSolidBins
            for ss=1:numberOfSecondSolidBins
                if f(s,ss)~=0 % If there are particles in the bin...
                    internalLiquid(s,ss)=min(granuleSaturationFactor*gasBins(s,ss),liquidBins(s,ss));
                    externalLiquid(s,ss)=max(0,liquidBins(s,ss)-internalLiquid(s,ss));
                    internalVolumeBins(s,ss)=vs(s)+vss(ss)+internalLiquid(s,ss)+gasBins(s,ss);
                    externalVolumeBins(s,ss)=vs(s)+vss(ss)+externalLiquid(s,ss)+gasBins(s,ss);
                    totalVolumeBins(s,ss)=vs(s)+vss(ss)+liquidBins(s,ss)+gasBins(s,ss);
                else
                    internalLiquid(s,ss)=0;
                    externalLiquid(s,ss)=0;
                    internalVolumeBins(s,ss)=0;
                    externalVolumeBins(s,ss)=0;
                    totalVolumeBins(s,ss)=0;
                end
            end
        end
        
        %         internalVolumeBins=s_meshxy+ss_meshxy+internalLiquid+gasBins;
        %         externalVolumeBins=s_meshxy+ss_meshxy+externalLiquid+gasBins;
        %         totalVolumeBins=s_meshxy+ss_meshxy+liquidBins+gasBins;
        % Bin Calculation______________
        
        liquidBinsAllCompartments(c,:,:)=liquidBins;
        gasBinsAllCompartments(c,:,:)=gasBins;
        internalVolumeBinsAllCompartments(c,:,:)=internalVolumeBins;
        externalVolumeBinsAllCompartments(c,:,:)=externalVolumeBins;
        totalVolumeBinsAllCompartments(c,:,:)=totalVolumeBins;
    end
    
    % Saving Over Time--------------------
    fAllCompartmentsOverTime(timeIndex,:,:,:)=fAllCompartments(:,:,:);
    flAllCompartmentsOverTime(timeIndex,:,:,:)=flAllCompartments(:,:,:);
    fgAllCompartmentsOverTime(timeIndex,:,:,:)=fgAllCompartments(:,:,:);
    liquidBinsAllCompartmentsOverTime(timeIndex,:,:,:)=liquidBinsAllCompartments(:,:,:);
    gasBinsAllCompartmentsOverTime(timeIndex,:,:,:)=gasBinsAllCompartments(:,:,:);
    internalVolumeBinsAllCompartmentsOverTime(timeIndex,:,:,:)=internalVolumeBinsAllCompartments(:,:,:);
    externalVolumeBinsAllCompartmentsOverTime(timeIndex,:,:,:)=externalVolumeBinsAllCompartments(:,:,:);
    totalVolumeBinsAllCompartmentsOverTime(timeIndex,:,:,:)=totalVolumeBinsAllCompartments(:,:,:);
    % Saving Over Time______________
    
    time
    time=time+timeStep;
    timeIndex=timeIndex+1;
end

%% POST TIME LOOP CALCULATION
d10OverTime(1:timeIndex,1:numberOfCompartments)=0;
d50OverTime(1:timeIndex,1:numberOfCompartments)=0;
d90OverTime(1:timeIndex,1:numberOfCompartments)=0;
totalVolumeAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments)=0;
totalSolidVolumeAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments)=0;
totalPoreVolumeAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments)=0;
totalLiquidVolumeAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments)=0;
totalGasVolumeAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments)=0;
volumeDistributionAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments,1:numberOfSieveGrid)=0;
totalSolidLeavingOverTime(1:timeIndex)=0;
totalStuffLeavingOverTime(1:timeIndex)=0;
totalLiquidLeavingOverTime(1:timeIndex)=0;
cumulativeVolumeDistributionAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments,1:numberOfSieveGrid)=0;
totalVolumeGridAllCompartmentsOverTime(1:timeIndex,1:numberOfCompartments,1:numberOfSieveGrid)=0;

for n=1:timeIndex-1
    for s=1:numberOfFirstSolidBins
        for ss=1:numberOfSecondSolidBins
            totalStuffLeavingOverTime(n)=totalStuffLeavingOverTime(n)+fAllCompartmentsOverTime(n,numberOfCompartments,s,ss)*((particleAverageVelocity)/distanceBetweenCompartments)*totalVolumeBinsAllCompartmentsOverTime(n,numberOfCompartments,s,ss)*timeStep;
            totalSolidLeavingOverTime(n)=totalSolidLeavingOverTime(n)+fAllCompartmentsOverTime(n,numberOfCompartments,s,ss)*((particleAverageVelocity)/distanceBetweenCompartments)*(vs(s)+vss(ss))*timeStep;
            totalLiquidLeavingOverTime(n)=totalLiquidLeavingOverTime(n)+fAllCompartmentsOverTime(n,numberOfCompartments,s,ss)*((particleAverageVelocity)/distanceBetweenCompartments)*(totalVolumeBinsAllCompartmentsOverTime(n,numberOfCompartments,s,ss)-vs(s)-vss(ss))*timeStep;
        end
    end
    
    for c=1:numberOfCompartments
        for s=1:numberOfFirstSolidBins
            for ss=1:numberOfSecondSolidBins
                totalVolumeAllCompartmentsOverTime(n,c)=totalVolumeAllCompartmentsOverTime(n,c)+fAllCompartmentsOverTime(n,c,s,ss)*totalVolumeBinsAllCompartmentsOverTime(n,c,s,ss);
                totalSolidVolumeAllCompartmentsOverTime(n,c)=totalSolidVolumeAllCompartmentsOverTime(n,c)+fAllCompartmentsOverTime(n,c,s,ss)*(vs(s)+vss(ss));
                totalPoreVolumeAllCompartmentsOverTime(n,c)=totalPoreVolumeAllCompartmentsOverTime(n,c)+fAllCompartmentsOverTime(n,c,s,ss)*(totalVolumeBinsAllCompartmentsOverTime(n,c,s,ss)-vs(s)-vss(ss));
                totalLiquidVolumeAllCompartmentsOverTime(n,c)=totalLiquidVolumeAllCompartmentsOverTime(n,c)+fAllCompartmentsOverTime(n,c,s,ss)*(liquidBinsAllCompartmentsOverTime(n,c,s,ss));
                totalGasVolumeAllCompartmentsOverTime(n,c)=totalGasVolumeAllCompartmentsOverTime(n,c)+fAllCompartmentsOverTime(n,c,s,ss)*(gasBinsAllCompartmentsOverTime(n,c,s,ss));
            end
        end
        f(1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
        f(:,:)=fAllCompartmentsOverTime(n,c,:,:);
        totalVolumeBins(:,:)=totalVolumeBinsAllCompartmentsOverTime(n,c,:,:);
        internalVolumeBins(:,:)=internalVolumeBinsAllCompartmentsOverTime(n,c,:,:);
        
        diameter=((6/pi)*totalVolumeBins).^(1/3)*1e6;
        
        totalVolumeGrid(1:numberOfSieveGrid)=0;
        
        for d=1:numberOfSieveGrid-1
            for s=1:numberOfFirstSolidBins
                for ss=1:numberOfSecondSolidBins
                    if and(diameter(s,ss)<sieveGrid(d+1),diameter(s,ss)>=sieveGrid(d))
                        totalVolumeGrid(d)=totalVolumeGrid(d)+f(s,ss)*totalVolumeBins(s,ss);
                    end
                end
            end
        end
        totalVolumeGridAllCompartmentsOverTime(n,c,:)=totalVolumeGrid;
        
        volumeDistribution(1:numberOfSieveGrid)=0;
        for d=1:numberOfSieveGrid
            volumeDistribution(d)=totalVolumeGrid(d)/sum(totalVolumeGrid);
        end
        cumulativeVolumeDistribution=cumsum(volumeDistribution);
        cumulativeVolumeDistributionAllCompartmentsOverTime(n,c,:)=cumulativeVolumeDistribution;
        volumeDistributionAllCompartmentsOverTime(n,c,:)=volumeDistribution;
        
        d10=0.1*(sieveGrid(2)/cumulativeVolumeDistribution(1));
        d50=0.5*(sieveGrid(2)/cumulativeVolumeDistribution(1));
        d90=0.9*(sieveGrid(2)/cumulativeVolumeDistribution(1));
        for i=2:numberOfSieveGrid
            if and(cumulativeVolumeDistribution(i-1)<0.5,cumulativeVolumeDistribution(i)>=0.5)
                loc=i;
                d50=((sieveGrid(loc)-sieveGrid(loc-1))/(cumulativeVolumeDistribution(loc)-cumulativeVolumeDistribution(loc-1)))*(0.5-cumulativeVolumeDistribution(loc-1))+sieveGrid(loc-1);%
            end
            if and(cumulativeVolumeDistribution(i-1)<0.1,cumulativeVolumeDistribution(i)>=0.1)
                loc10=i;
                d10=((sieveGrid(loc10)-sieveGrid(loc10-1))/(cumulativeVolumeDistribution(loc10)-cumulativeVolumeDistribution(loc10-1)))*(0.1-cumulativeVolumeDistribution(loc10-1))+sieveGrid(loc10-1);%
            end
            if and(cumulativeVolumeDistribution(i-1)<0.9,cumulativeVolumeDistribution(i)>=0.9)
                loc90=i;
                d90=((sieveGrid(loc90)-sieveGrid(loc90-1))/(cumulativeVolumeDistribution(loc90)-cumulativeVolumeDistribution(loc90-1)))*(0.9-cumulativeVolumeDistribution(loc90-1))+sieveGrid(loc90-1);%
            end
        end
        
        d10OverTime(n,c)=d10;
        d50OverTime(n,c)=d50;
        d90OverTime(n,c)=d90;
    end
end
toc

%% PLOTS
%  Diameter------------------------------
figure;
plot(Time', d50OverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
hold on
plot(Time', d50OverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
hold on
plot(Time', d50OverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('D_{50} (\mu\it{m})');
title({'Dynamic Particle Size'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
%  Diameter_____________________

% Total Volume---------------------------
figure;
plot(Time', totalVolumeAllCompartmentsOverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
hold on
plot(Time', totalVolumeAllCompartmentsOverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
hold on
dplot(Time', totalVolumeAllCompartmentsOverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('Volume (m^{3})');
title({'Total Volume'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
% Total Volume___________________

% Solid Volume---------------------------
figure;
plot(Time', totalSolidVolumeAllCompartmentsOverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
hold on
plot(Time', totalSolidVolumeAllCompartmentsOverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
hold on
plot(Time', totalSolidVolumeAllCompartmentsOverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('Volume (m^{3})');
title({'Total Solid Volume'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
% Solid Volume___________________

% Pore Volume
% figure;
% plot(Time', totalPoreVolumeAllCompartmentsOverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
% hold on
% plot(Time', totalPoreVolumeAllCompartmentsOverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
% hold on
% plot(Time', totalPoreVolumeAllCompartmentsOverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
% hold on
%
% hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
% hXLabel = xlabel('Time (seconds)');
% hYLabel = ylabel('Volume (m^{3})');
% % title({'Total Pore Volume'},'FontName','Georgia','FontSize',14)
% set( gca,'FontName' , 'garamond' ); % Legend Font
% set([hLegend, gca],'FontSize',12); % Legend font-size
% set(hLegend,'box','off') % Hiding legend box
% set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
% set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)

% Liquid Volume--------------------------
figure;
plot(Time', totalLiquidVolumeAllCompartmentsOverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
hold on
plot(Time', totalLiquidVolumeAllCompartmentsOverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
hold on
plot(Time', totalLiquidVolumeAllCompartmentsOverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('Volume (m^{3})');
title({'Total Liquid Volume'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
% Liquid Volume__________________

% Gas Volume----------------------------
figure;
plot(Time', totalGasVolumeAllCompartmentsOverTime(1:timeIndex-1,1),'LineWidth',2,'Color','k');
hold on
plot(Time', totalGasVolumeAllCompartmentsOverTime(1:timeIndex-1,2),'LineWidth',2,'Color','r');
hold on
plot(Time', totalGasVolumeAllCompartmentsOverTime(1:timeIndex-1,3),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Compartment_{1}','Compartment_{2}','Compartment_{3}','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('Volume (m^{3})');
title({'Total Gas Volume'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
% Gas Volume___________________

% Leaving--------------------------------
figure;
plot(Time', totalStuffLeavingOverTime(1:timeIndex-1),'LineWidth',2,'Color','k');
hold on
plot(Time', totalSolidLeavingOverTime(1:timeIndex-1),'LineWidth',2,'Color','r');
hold on
plot(Time', totalLiquidLeavingOverTime(1:timeIndex-1),'LineWidth',2,'Color','b');
hold on

hLegend = legend('Total','Solid','Liquid+Gas','Location','NorthWest','FontName','Georgia','FontSize',14);
hXLabel = xlabel('Time (seconds)');
hYLabel = ylabel('Volume (m^{3})');
title({'Volume Leaving'},'FontName','Georgia','FontSize',14)
set( gca,'FontName' , 'garamond' ); % Legend Font
set([hLegend, gca],'FontSize',12); % Legend font-size
set(hLegend,'box','off') % Hiding legend box
set([hXLabel, hYLabel],'FontName', 'Georgia','FontSize',16); % Axes Label font name & size
set(gca, 'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth' ,1)
% Leaving______________________

%% DEBUGGING
% Liquid Content, Internal Liquid Content, Porosity Calculation
liquidContent(1:timeIndex,1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
internalLiquidContent(1:timeIndex,1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
porosity(1:timeIndex,1:numberOfCompartments,1:numberOfFirstSolidBins,1:numberOfSecondSolidBins)=0;
for n=1:timeIndex-1
    for c=1:numberOfCompartments
        for s=1:numberOfFirstSolidBins
            for ss=1:numberOfSecondSolidBins
                liquidContent(n,c,s,ss)=(liquidBinsAllCompartmentsOverTime(n,c,s,ss)/(vs(s)+vss(ss)+liquidBinsAllCompartmentsOverTime(n,c,s,ss)))*100;
                internalLiquidContent(n,c,s,ss)=((internalVolumeBinsAllCompartmentsOverTime(n,c,s,ss)-(gasBinsAllCompartmentsOverTime(n,c,s,ss)+vs(s)+vss(ss)))/internalVolumeBinsAllCompartmentsOverTime(n,c,s,ss))*100;
                porosity(n,c,s,ss)=(gasBinsAllCompartmentsOverTime(n,c,s,ss)/(vs(s)+vss(ss)+gasBinsAllCompartmentsOverTime(n,c,s,ss)))*100;
            end
        end
    end
end
% Liquid Content, Internal Liquid Content, Porosity Calculation

% Stuff Coming In------------------------
StuffEntering=0;
for s=1:numberOfFirstSolidBins
    for ss=1:numberOfSecondSolidBins
        StuffEntering=StuffEntering+fIn(s,ss)*(vs(s)+vss(ss));
    end
end
% Stuff Coming In_________________

% Reshaping F ---------------------------
f1=reshape(fAllCompartments(1,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
f2=reshape(fAllCompartments(2,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
f3=reshape(fAllCompartments(3,:,:),[numberOfFirstSolidBins numberOfSecondSolidBins]);
% Reshaping F___________________
