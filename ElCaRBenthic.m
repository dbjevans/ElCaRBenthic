%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ElCaRBenthic v1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute temperature and/or Omega+/-[CO32-] from benthic foraminifera
% Mg/Ca, B/Ca, Sr/Ca, or Mg/Li data
% Monte Carlo uncertainty propagation can account for proxy value, 
% calibration coefficient, and secondary parameter (MgCaSW, T, Omega, 
% salinity, depth) uncertainties
%
% [TempOut,OmegaOut,CO3Out] = ...
%    ElCaRBenthic(sampleage,sampledata,speciesID,elementID,Hin,temperature,...
%    Omega,salinity,depth,plotFig)
%
% All input fields are required.
% - sample age: in ka
% - sample data: either one, two, four, or six columns (e.g. [Mg/Ca Sr/Ca], 
%   Mg/Ca only, [Mg/Ca Mg/CaUncertainty B/Ca B/CaUncertainty], or 
%   [Mg/Ca Mg/CaUncertainty Li/Ca Li/CaUncertainty Sr/Ca SrCaUncertainty]. 
%   The order is unimportant EXCEPT that Mg/Ca must be followed directly by
%   Li/Ca (rowise), e.g. [Li/Ca Mg/Ca Sr/Ca] will cause an error. Mg/Li
%   data can be input directly, but changes in seawater chemistry will not
%   be accounted for in this case.
% - species ID: (1 = L. wuellerstorfi, 2 = Uvigerina spp., 3 = O. umbonatus, 
%   4 = C. mundulus, 5 = C. pachyderma, 6 = N. umbonifera)
% - elementID: (1 = B/Ca, 2 = Sr/Ca, 3 = Mg/Ca, 4 = Mg/Li. Note: separating 
%   Mg/Ca and Li/Ca data should still be denoted using 4, e.g. passing 
%   Mg/Ca Li/Ca and Sr/Ca data would be denoted [4 1])
% - Hin: The value of H, where H is the power coefficient of the 
%   relationship between seawater and shell Mg/Ca:
%       Mg/Ca(corrected) = Mg/Ca(measured)*5.2^H/MgCaSW^H
%   And may be a single value, or a value and uncertainty. Enter NaN for
%   records with the oldest sample <800 ka.
% - temperature: required if calculating Omega from Sr/Ca using a species
%   that has a Sr/Ca-T sensitivity (and no Mg/Ca data is provided) OR if 
%   calculation of [CO32-] from Omega is desired. Can be a single number 
%   (with or without uncertainty), or a vector the same size as sample
%   data. Set as NaN when calculating temperature using Mg/Ca.
% - Omega: required if calculating temperature from Mg/Ca using a species
%   with a Mg/Ca-Omega sensitivity if no B/Ca or Sr/Ca data are provided,
%   or to calculate temperature-corrected [CO32-] assuming constant Omega.
% - salinity: required to calculate [CO32-] from Omega + T. Can be a single
%   value with or without uncertainty, or a vector the same size as sample
%   data.
% - depth: required to calculate [CO32-] from Omega + T. Can be a single
%   value with or without uncertainty, or a vector the same size as sample
%   data.
% - plotFig: plot the results (1/0)
%
%
% EXAMPLES:
%
% [TempOut,OmegaOut,CO3Out] = ...
%    ElCaRBenthic(sampleage,sampledata,1,[1 3],NaN,NaN,NaN,[35 1],[2500 100],1)
%
% Calculate temperature, Omega, and [CO32-] using L. wuellerstorfi B/Ca and
% Mg/Ca. Derive [CO32-] from Omega and temperature assuming a salinity of
% 35+/-1, and a pressure of 2500+/-100 dbar, and plot the results.
%
%
% [TempOut,OmegaOut,CO3Out] = ...
%    ElCaRBenthic(sampleage,sampledata,5,2,NaN,tempIn,NaN,[35 1],[2500 100],0)
%
% Calculate Omega, and [CO32-] using C. pachyderma Sr/Ca using an
% existing temperature record to correct for the temperature contol on 
% Sr/Ca for this species. Derive [CO32-] from Omega and temperature 
% assuming a salinity of 35+/-1, and a pressure of 2500+/-100 dbar, don't 
% plot the results.
%
%
% [TempOut,OmegaOut,CO3Out] = ...
%    ElCaRBenthic(sampleage,sampledata,4,[1 3],[0.54 0.1],tempIn,NaN,[35 1],...
%    [2500 100],1)
%
% Calculate temperature, Omega, and [CO32-] using O. umbonatus B/Ca and 
% Mg/Ca. Use an O. umbonatus H = 0.54+/-0.1 (Evans et al., 2016). Derive 
% [CO32-] from Omega and temperature assuming a salinity of 35+/-1, and a 
% pressure of 2500+/-100 dbar, and plot the results.


function [TempOut,OmegaOut,CO3Out] = ...
    ElCaRBenthic(sampleage,sampledata,speciesID,elementID,Hin,temperature,...
    Omega,salinity,depth,plotFig)

% Coming in a future release?
%   Lookup modern T/S
%   relative T corrected to coretop
%   Correct [CO32-] for [Casw]


% Useful for troubleshooting
% sampleage = ;
% ampledata = ;
% speciesID = 5;
% elementID = [4 2];
% Hin = [0.5 0.2];
% temperature = NaN;
% Omega = [1 1];
% salinity = [35 1];
% depth = [3500 500];
% plotFig = 1;

% input error/warning messages
if size(sampleage,1)~=size(sampledata,1) || size(temperature,1)~=1 && ...
        size(sampleage,1)~=size(temperature,1)
    error('input datasets are not the same length')
end
if size(unique(sum(isnan(sampledata))),2)~=1
    warning(['unequal nunber of NaN values in the sample data ',...
        '- only complete rows will be computed'])
end
if speciesID>6 || any(elementID>4)
    error('error in speciesID or elementID input')
end
if size(elementID,2)>1 && (size(sampledata,2)<2 || size(sampledata,2)==5)
    error('incorrect number of element data files')
end

% remove NaNs before plotting
sampleageCalc = sampleage(sum(isnan(sampledata),2)==0,:);
sampledata = sampledata(sum(isnan(sampledata),2)==0,:);

% load required lookup data
load("ElCaRBenthic_lookup_data.mat","co2sysOut","foramData","foramDataWOAT",...
    "swMgCaComposite","LiCompOut","CaCompOut")

% seawater Mg/Ca data only required for samples at least 800 ka
if max(sampleage)>800
    mgcaSW = interp1(swMgCaComposite(:,1),swMgCaComposite(:,2:4),sampleageCalc);
    caSW = interp1((0:0.1:66)'.*1e3,CaCompOut,sampleageCalc);
    licaSW = interp1((0:0.1:66)'.*1e3,LiCompOut,sampleageCalc)./...
        [caSW(:,3) caSW(:,2) caSW(:,1)].*1e6;
else
    mgcaSW = ones(size(sampleageCalc,1),3).*5.2;
    licaSW = ones(size(sampleageCalc,1),3).*2.82;
    caSW = ones(size(sampleageCalc,1),3).*10.3;
end

% Generate dataset uncertainties if none given
% This will fail if data inputs are inconsistent
%   (e.g. three column input consisting of Mg/Ca+unc + Li/Ca)
% ! Either all or no data uncertainties must be given !
if size(elementID,2)==1
    % if no data uncertainty given, set to zero
    if size(sampledata,2)==1
        sampledata(:,2) = 0;
    end
    % If Mg/Ca and Li/Ca fed in but no uncertainty given
    if elementID==4
        % If Mg/Ca and Li/Ca fed in seperately, sampledata
        % could be two or four columns (if with uncertainty), distinguished
        % from input of Mg/Li with uncertainty in the former case by
        % comparing the means of the first two columns
        if size(sampledata,2)==2 && ...
            mean(sampledata(:,2),'omitnan') > mean(sampledata(:,1),'omitnan')
    
            sampledata = [sampledata(:,1) zeros(size(sampledata,1),1) ...
                sampledata(:,2) zeros(size(sampledata,1),1)];
    
        end
    end
    % Mg/Li is now either two column (Mg/Li+unc) or four column

else
    % If two datasets (e.g. Mg/Ca and Sr/Ca)
    if size(sampledata,2)==2
        sampledata = [sampledata(:,1) zeros(size(sampledata,1),1) ...
            sampledata(:,2) zeros(size(sampledata,1),1)];
    % If three datasets (e.g. Mg/Ca, Li/Ca, and Sr/Ca)
    elseif size(sampledata,2)==3
        sampledata = [sampledata(:,1) zeros(size(sampledata,1),1) ...
            sampledata(:,2) zeros(size(sampledata,1),1) ...
            sampledata(:,3) zeros(size(sampledata,1),1)];
    end
    % If size(sampledata,2)==4 | 6 then uncertainties are already given
end
% Sample data is now either 2, 4, or 6 columns

% add in missing datasets, if missing
% if no salinity uncertainty, set to zero
if size(salinity,2)==1
    salinity(:,2) = 0;
end
% if only one salinity given, fill down
if size(salinity,1)==1
    salinity = ones(size(sampledata,1),2).*salinity;
end
% if no depth uncertainty, set to zero
if size(depth,2)==1
    depth(:,2) = 0;
end
% if only one depth given, fill down
if size(depth,1)==1
    depth = ones(size(sampledata,1),2).*depth;
end
% If no H uncertainty given, set to zero
if size(Hin,2)==1
    Hin(1,2) = 0;
end
% Rename others
tempOmega = Omega;
tempIn = temperature;


%%% Perform calculations
if size(elementID,2)==1 % if only one dataset present
    if elementID==1 || elementID==2 % If it is B/Ca or Sr/Ca data

        % when passing Sr or B data, and if no temperature uncertainty, set to zero
        if size(tempIn,2)==1
            tempIn(:,2) = 0;
        end
        % when passing Sr or B data, if only one temperature given, fill down
        if size(tempIn,1)==1
            tempIn = ones(size(sampledata,1),2).*tempIn;
        end

        % Calculate Omega and CO3 from shell data and constant T
        [OmegaOut,CO3Out] = calcOmega(speciesID,elementID,foramDataWOAT,...
            foramData,co2sysOut,tempIn,depth,salinity,caSW,sampledata);

        % Reformat input temperature data for plotting
        TempOut = [NaN(size(OmegaOut,1),1) tempIn(:,1)-tempIn(:,2) ....
            NaN(size(OmegaOut,1),1) tempIn(:,1) NaN(size(OmegaOut,1),1) ...
            tempIn(:,1)+tempIn(:,2) NaN(size(OmegaOut,1),1)];

    else
        % when passing Mg(Li) data, and if no Omega uncertainty, set to zero
        if size(tempOmega,2)==1
            tempOmega(:,2) = 0;
        end
        % when passing Mg(Li) data, if only one Omega given, fill down
        if size(tempOmega,1)==1
            tempOmega = ones(size(sampledata,1),2).*tempOmega;
        end

        if elementID==3 % If it is Mg/Ca data
    
            % Calculate Temperature from shell data
            % Note this requires 1/Omega^2 rather than Omega as an input
            TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                co2sysOut,tempOmega,sampledata,3,mgcaSW,licaSW,Hin);
    
            % Reformat input Omega data for plotting
            OmegaOut = [NaN(size(TempOut,1),1) tempOmega(:,1)-tempOmega(:,2) ....
                NaN(size(TempOut,1),1) tempOmega NaN(size(TempOut,1),1) ...
                tempOmega(:,1)+tempOmega(:,2) NaN(size(TempOut,1),1)];
            CO3Out = NaN(size(OmegaOut));
    
        elseif elementID==4 % If it is Mg/Li data
    
            % Calculate Temperature from shell data
            % Note this requires 1/Omega^2 rather than Omega as an input
            TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                co2sysOut,tempOmega,sampledata,4,mgcaSW,licaSW,Hin);
    
            % Reformat input Omega data for plotting
            OmegaOut = [NaN(size(TempOut,1),1) tempOmega(:,1)-tempOmega(:,2) ....
                NaN(size(TempOut,1),1) tempOmega NaN(size(TempOut,1),1) ...
                tempOmega(:,1)+tempOmega(:,2) NaN(size(TempOut,1),1)];
            CO3Out = NaN(size(OmegaOut));
        end
    end

else % if more than one dataset present (Mg/Ca or Mg/Li + Sr/Ca or B/Ca)

    % Find location of B or Sr data
    elLoc = [0 0];
    elLoc(elementID<3) = 1;
    elLoc(2,:) = [1 2];
    elLoc = sortrows(elLoc',1,"descend")';
    if elLoc(2,1)==2
        if size(sampledata,2)==4 % If two datasets were fed in
            OmegaElData = sampledata(:,3:4);
            TempElData = sampledata(:,1:2);
        else % If three (Mg/Ca, Li/Ca, and B or Sr/Ca)
            OmegaElData = sampledata(:,5:6);
            TempElData = sampledata(:,1:4);
        end
    else
        if size(sampledata,2)==4
            OmegaElData = sampledata(:,1:2);
            TempElData = sampledata(:,3:4);
        else
            OmegaElData = sampledata(:,1:2);
            TempElData = sampledata(:,3:6);
        end
    end
    % Extract ID of relevant Omega-sensitive element
    elID = elementID(1,elLoc(2,1));

    % if B/Ca + Mg/Ca or Mg/Li are input
    if sum(elementID==1)~=0
        % B/Ca not sensitive to T, so calculate Omega first
        
        % Initially calculate Omega and CO3 from shell data and constant T
        [OmegaOut,~] = calcOmega(speciesID,elID,foramDataWOAT,foramData,...
            co2sysOut,ones(size(sampledata,1),2),depth,salinity,caSW,OmegaElData);
        % Then calculate T using Omega where relevant
        % Note this requires 1/Omega^2 rather than Omega as an input
        % Note this uses 95% calibration-derived CI as the Omega unc.
        if sum(elementID==3)~=0 % Mg/Ca
            TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                co2sysOut,[OmegaOut(:,4) OmegaOut(:,6)-OmegaOut(:,4)],TempElData,...
                3,mgcaSW,licaSW,Hin);
        else % Mg/Li
            TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                co2sysOut,[OmegaOut(:,4) OmegaOut(:,6)-OmegaOut(:,4)],TempElData,...
                4,mgcaSW,licaSW,Hin);
        end
        % Recalculate CO3 using Mg/Ca-T
        [OmegaOut,CO3Out] = calcOmega(speciesID,elID,foramDataWOAT,foramData,...
            co2sysOut,[TempOut(:,4) TempOut(:,6)-TempOut(:,4)],depth,...
            salinity,caSW,OmegaElData);
    
    else % If Sr/Ca + Mg/Ca or Mg/Li are input

        % If C. pachyderma, Sr/Ca is sensitive to T, Mg/Ca is not 
        % If O. umbonatus or Uvigerina, order makes no difference
        % sensitive to Omega --> calculate T first
        if speciesID==2 || speciesID==3 || speciesID==5
            % Calculate T
            if sum(elementID==3)~=0 % Mg/Ca
                TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                    co2sysOut,[ones(size(TempElData,1),2)],TempElData,...
                    3,mgcaSW,licaSW,Hin);
            else % Mg/Li
                TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                    co2sysOut,[ones(size(TempElData,1),2)],TempElData,...
                    4,mgcaSW,licaSW,Hin);
            end
            % Calculate Omega and CO3 using Mg/Ca-T
            [OmegaOut,CO3Out] = calcOmega(speciesID,elID,foramDataWOAT,foramData,...
                co2sysOut,[TempOut(:,4) TempOut(:,6)-TempOut(:,4)],depth,...
                salinity,caSW,OmegaElData);

        % If L. wuellerstorfi or C. mundulus, Mg/Ca is sensitive to Omega,
        % Sr/Ca is not sensitive to T, so calculate Omega first
        else

            % Initially calculate Omega and CO3 from shell data and constant T
            [OmegaOut,~] = calcOmega(speciesID,elID,foramDataWOAT,foramData,...
                co2sysOut,ones(size(sampledata,1),2),depth,salinity,...
                caSW,OmegaElData);
            % Then calculate T using Omega
            % Note this requires 1/Omega^2 rather than Omega as an input
            % Note this uses 95% calibration-derived CI as the Omega unc.
            if sum(elementID==3)~=0 % Mg/Ca
                TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                    co2sysOut,[OmegaOut(:,4) OmegaOut(:,6)-OmegaOut(:,4)],...
                    TempElData,3,mgcaSW,licaSW,Hin);
            else % Mg/Li
                TempOut = calcTemp(sampleageCalc,speciesID,foramDataWOAT,foramData,...
                    co2sysOut,[OmegaOut(:,4) OmegaOut(:,6)-OmegaOut(:,4)],...
                    TempElData,4,mgcaSW,licaSW,Hin);
            end
            % Recalculate CO3 using Mg/Ca-T (necessary for C. pachyderma)
            [OmegaOut,CO3Out] = calcOmega(speciesID,elID,foramDataWOAT,foramData,...
                co2sysOut,[TempOut(:,4) TempOut(:,6)-TempOut(:,4)],depth,...
                salinity,caSW,OmegaElData);

        end
    end
end


% Plot figure
if plotFig==1
    cmap = parula(12);

    close(figure)
    F1 = figure;
    set(F1,'PaperUnits','centimeter','units','centimeter',...
        'papersize',[20 6.5],'Position',[20 16 20 6.5],'color',[1 1 1])
    t = tiledlayout(1,3);
    
    % Temperature
    t1 = nexttile;
    %plot uncertainties and best estimate
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [TempOut(:,2) ; flipud(TempOut(:,6))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    hold on
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [TempOut(:,3) ; flipud(TempOut(:,5))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    % fill([sampleageCalc ; flipud(sampleageCalc)],...
    %     [OmegaOut(:,1) ; flipud(OmegaOut(:,7))], ...
    %     cmap(4,:),'facealpha',0.1,'edgecolor','none')
    plot(sampleageCalc,TempOut(:,4),'-','color',cmap(4,:))
    plot(sampleageCalc,TempOut(:,1),'--','color',cmap(4,:))
    plot(sampleageCalc,TempOut(:,7),'--','color',cmap(4,:))

    xlabel('age (ka)')
    ylabel('temperature (\circC)')

    % Omega
    t2 = nexttile;
    %plot uncertainties and best estimate
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [OmegaOut(:,2) ; flipud(OmegaOut(:,6))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    hold on
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [OmegaOut(:,3) ; flipud(OmegaOut(:,5))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    % fill([sampleageCalc ; flipud(sampleageCalc)],...
    %     [OmegaOut(:,1) ; flipud(OmegaOut(:,7))], ...
    %     cmap(4,:),'facealpha',0.1,'edgecolor','none')
    plot(sampleageCalc,OmegaOut(:,4),'-','color',cmap(4,:))
    plot(sampleageCalc,OmegaOut(:,1),'--','color',cmap(4,:))
    plot(sampleageCalc,OmegaOut(:,7),'--','color',cmap(4,:))

    legend('','calibration CI','best estimate','IPI','',...
        'location','northeast','fontsize',6,'box','on')
    
    xlabel('age (ka)')
    ylabel('\Omega_{calcite}')
    yL = get(gca,'ylim');
    ylim([0 yL(2)])
    
    % CO3
    t3 = nexttile;
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [CO3Out(:,2) ; flipud(CO3Out(:,6))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    hold on
    fill([sampleageCalc ; flipud(sampleageCalc)],...
        [CO3Out(:,3) ; flipud(CO3Out(:,5))], ...
        cmap(4,:),'facealpha',0.25,'edgecolor','none')
    % fill([sampleageCalc ; flipud(sampleageCalc)],...
    %     [CO3Out(:,1) ; flipud(CO3Out(:,7))], ...
    %     cmap(4,:),'facealpha',0.1,'edgecolor','none')
    plot(sampleageCalc,CO3Out(:,4),'-','color',cmap(4,:))
    plot(sampleageCalc,CO3Out(:,1),'--','color',cmap(4,:))
    plot(sampleageCalc,CO3Out(:,7),'--','color',cmap(4,:))

    xlabel('age (ka)')
    ylabel('[CO_3^{2-}] (\mumol/kg)')
    yL = get(gca,'ylim');
    ylim([0 yL(2)])
    
    set(t1,'box','on','layer','top','fontsize',8)
    set(t2,'box','on','layer','top','fontsize',8)
    set(t3,'box','on','layer','top','fontsize',8)
    
    t.TileSpacing = "Tight";
    t.Padding = "tight";
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions

% Define a custom function that fits a linear model and returns coefficients
function coefficients = fitlm_coefficients(data)
    X = data(:, 1:end-1); % All columns except the last one
    y = data(:, end);     % The last column
    mdl = fitlm(X, y, 'linear');
    coefficients = mdl.Coefficients.Estimate'; % Return coefficients as a row vector
end

% Derive Ksp for calcite for all samples (from co2sys)
% pressure in dbar, temperature in oC, salinity on usual scale
function Ksp = CaSolubility(Sal,TempC,Pdbar)
    RGasConstant = 83.14462618; % ml bar-1 K-1 mol-1,
    RT = RGasConstant.*(TempC + 273.15);
    TempK = TempC + 273.15;
    
    % CalciteSolubility:
    % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKCa = -171.9065 - 0.077993.*TempK + 2839.319./TempK;
    logKCa = logKCa + 71.595.*log(TempK)./log(10);
    logKCa = logKCa + (-0.77712 + 0.0028426.*TempK + 178.34./TempK).*Sal.^0.5;
    logKCa = logKCa - 0.07711.*Sal+ 0.0041249.*Sal.^0.5.*Sal;
    % '       sd fit = .01 (for Sal part, not part independent of Sal)
    KCa = 10.^(logKCa);% ' this is in (mol/kg-SW)^2
    
    % PressureCorrectionForCalcite:
    % '       Ingle, Marine Chemistry 3:301-319, 1975
    deltaVKCa = -48.76 + 0.5304.*TempC;
    KappaKCa  = (-11.76 + 0.3692.*TempC)./1000;
    lnKCafac  = (-deltaVKCa + 0.5.*KappaKCa.*(Pdbar./10)).*(Pdbar./10)./RT;
    Ksp       = KCa.*exp(lnKCafac);
end

%%% Calculate Omega & [CO32-] from B/Ca or Sr/Ca
function [OmegaOut,CO3Out] = calcOmega(speciesID,elementID,foramDataWOAT,...
    foramData,co2sysOut,temperature,depth,salinity,CaSW,sampledata)

    tempC = temperature;
    sal = salinity;
    % Only include unflagged analyses in calibration
    isFlag = isnan(foramData.flagged_mgca);

    %%% find relevant portion of the database
    if speciesID==1 % L. wuellerstorfi
        DataSubset = any((strncmp(foramData.species,"l.wue",5) & isFlag),2);
    elseif speciesID==2 % Uvigerina spp.
        DataSubset = any((strncmp(foramData.species,"u.per",5) | ...
            strncmp(foramData.species,"uvi",3) & isFlag),2);
    elseif speciesID==3 % O. umbonatus
        DataSubset = any((strncmp(foramData.species,"o.umb",5) & isFlag),2);
    elseif speciesID==4 % C. mundulus
        DataSubset = any((strncmp(foramData.species,"c.mun",5) & isFlag),2);
    elseif speciesID==5 % C. pachyderma
        DataSubset = any((strncmp(foramData.species,"c.pac",5) & isFlag),2);
    elseif speciesID==6 % N. umbonifera
        DataSubset = any((strncmp(foramData.species,"n.umb",5) & isFlag),2);
    end

    %%% Calculate coefficient covariance for relevant species
    % If C. pachyderma Sr/Ca, include T in regression
    if speciesID==5 && elementID==2
        tempx = [foramDataWOAT(DataSubset) co2sysOut(DataSubset,17)];
    else
        tempx = co2sysOut(DataSubset,17);
    end
    if elementID==1
        tempy = foramData.bca(DataSubset);
    else
        tempy = foramData.srca(DataSubset);
    end
    % Combine X and y into a single matrix for bootstrapping
    data = [tempx, tempy];
    % best model estimate
    refFit = fitlm(tempx,tempy,'linear');
    % Perform bootstrapping
    bootCoefficients = bootstrp(1000, @fitlm_coefficients, data);
    % fit to bootstrap coefficients
    % If C. pachyderma Sr/Ca, Omega is the third coefficient
    if speciesID==5 && elementID==2
        % covariation is strongest between T and Omega coefficients rather 
        % than Omega and intercept in this case
        fitboostrap = ...
            fitlm(bootCoefficients(:,2),bootCoefficients(:,3)); % m = x1 + x2*int

        % sample coefficients 10000 times across MgCa range and output results
        ElCaRand = repmat(sampledata(:,1),1,10000);
        % Add random measurement uncertainties
        ElCaRand = ElCaRand + normrnd(0,0.5,size(ElCaRand)).*...
            repmat(sampledata(:,2),1,size(ElCaRand,2));
        
        % sample model coefficients, use SD=1 as fitlm returns 1s
        % uncertainties
        coefRand = normrnd(0,1,[4 size(ElCaRand,2)]);
        % calculate all possible Omega
        tempTsens = refFit.Coefficients{2,1} + ... % Tsens = x1 + X2*Osens
            refFit.Coefficients{2,2}.*repmat(coefRand(1,:),size(ElCaRand,1),1);
        tempOsens = tempTsens.*(fitboostrap.Coefficients{2,1} + ...
            fitboostrap.Coefficients{2,2}.*repmat(coefRand(2,:),size(ElCaRand,1),1)) + ...
            (fitboostrap.Coefficients{1,1} + ...
            fitboostrap.Coefficients{1,2}.*repmat(coefRand(3,:),size(ElCaRand,1),1));
        
        %%% Omega calculations
        % Randomly sample T within uncertainty
        tempCresamp = tempC(:,1) + tempC(:,2).*normrnd(0,0.5,size(ElCaRand));
        % Omega, with uncertainty
        Omegaoutlinear = (ElCaRand - ...
            ((refFit.Coefficients{1,1} + refFit.Coefficients{1,2}.*...
            repmat(coefRand(4,:),size(ElCaRand,1),1)) + ...
            tempCresamp.*tempTsens))./tempOsens;
        % Omega, no uncertainty
        Refoutlinear = (sampledata(:,1) - ...
            (refFit.Coefficients{1,1} + refFit.Coefficients{2,1}.*tempC(:,1)))./...
            refFit.Coefficients{3,1};
        % Calculate inverse prediction intervals
        tempY = (tempy - ...
            (refFit.Coefficients{1,1} + refFit.Coefficients{2,1}.*...
            foramDataWOAT(DataSubset)))./...
            refFit.Coefficients{3,1};
        tempX = tempx(:,2);
        IPIOmega = calculate_IPI(0,[tempX(~isnan(tempY)) tempY(~isnan(tempY))]);
    else
        fitboostrap = ...
            fitlm(bootCoefficients(:,1),bootCoefficients(:,2)); % m = x1 + x2*int
        
        % sample coefficients 10000 times across MgCa range and output results
        ElCaRand = repmat(sampledata(:,1),1,10000);
        % Add random measurement uncertainties
        ElCaRand = ElCaRand + normrnd(0,0.5,size(ElCaRand)).*...
            repmat(sampledata(:,2),1,size(ElCaRand,2));
        
        % sample model coefficients
        coefRand = normrnd(0,1,[3 size(ElCaRand,2)]);
        % calculate all possible Omega
        tempInt = refFit.Coefficients{1,1} + ... % Omega coef = x1 + X2*B
            refFit.Coefficients{1,2}.*repmat(coefRand(1,:),size(ElCaRand,1),1);
        tempSens = tempInt.*(fitboostrap.Coefficients{2,1} + ...
            fitboostrap.Coefficients{2,2}.*repmat(coefRand(2,:),size(ElCaRand,1),1)) + ...
            (fitboostrap.Coefficients{1,1} + ...
            fitboostrap.Coefficients{1,2}.*repmat(coefRand(3,:),size(ElCaRand,1),1));
        
        %%% Omega calculations
        % Omega, with uncertainty
        Omegaoutlinear = (ElCaRand - tempInt)./tempSens;
        % Omega, no uncertainty
        Refoutlinear = (sampledata(:,1) - refFit.Coefficients{1,1})./...
            refFit.Coefficients{2,1};
        % Calculate inverse prediction intervals
        tempY = (tempy - ...
            refFit.Coefficients{1,1})./refFit.Coefficients{2,1};
        tempX = tempx;
        IPIOmega = calculate_IPI(0,[tempX(~isnan(tempY)) tempY(~isnan(tempY))]);
    end

    % Extract output data
    OmegaOut = [Refoutlinear-IPIOmega ...
        prctile(Omegaoutlinear,2.5,2) prctile(Omegaoutlinear,16,2) ...
        Refoutlinear prctile(Omegaoutlinear,84,2) ...
        prctile(Omegaoutlinear,97.5,2) Refoutlinear+IPIOmega];
    
    %%% Convert Omega to [CO32-]
    % Add random uncertainties
    KspSal = repmat(sal(:,1),1,size(ElCaRand,2)) + ...
        repmat(sal(:,2),1,size(ElCaRand,2)).*...
        repmat(normrnd(0,0.5,1,size(ElCaRand,2)),size(sal,1),1);
    % if C. pachyderma Sr/Ca, T is required in Sr/Ca-Omega regression and
    % has therefore already been resampled above
    if speciesID==5 && elementID==2
        KspTemp = tempCresamp;
    else
        KspTemp = repmat(tempC(:,1),1,size(ElCaRand,2)) + ...
            repmat(tempC(:,2),1,size(ElCaRand,2)).*...
            repmat(normrnd(0,0.5,1,size(ElCaRand,2)),size(tempC,1),1);
    end
    KspPres = repmat(depth(:,1),1,size(ElCaRand,2)) + ...
        repmat(depth(:,2),1,size(ElCaRand,2)).*...
        repmat(normrnd(0,0.5,1,size(ElCaRand,2)),size(depth,1),1);
    % Calculate Ksp
    Ksp = CaSolubility(KspSal,KspTemp,KspPres);
    % Convert Omega to [CO32-] using Ksp and modern [Ca2+]
    CO3 = Omegaoutlinear.*Ksp./0.0103*1e6;
    % Ksp with no uncertainty
    KspBest = CaSolubility(sal(:,1),tempC(:,1),depth(:,1));
    % Best estimate and best estimate +/-IPI converted to [CO32-]
    CO3Best = [Refoutlinear.*KspBest./0.0103*1e6 ...
        (Refoutlinear-IPIOmega).*KspBest./0.0103*1e6 ...
        (Refoutlinear+IPIOmega).*KspBest./0.0103*1e6];
    % Avoids negative root issues
    CO3 = real(CO3);
    
    % Extract output data
    CO3Out = [CO3Best(:,2) ...
        prctile(CO3,2.5,2) prctile(CO3,16,2) CO3Best(:,1) ...
        prctile(CO3,84,2) ...
        prctile(CO3,97.5,2) CO3Best(:,3)];

    % Adjust for seawater [Ca] changes
    CO3Out(:,1:3) = CO3Out(:,1:3).*10.3./CaSW(:,3);
    CO3Out(:,4) = CO3Out(:,4).*10.3./CaSW(:,2);
    CO3Out(:,5:7) = CO3Out(:,5:7).*10.3./CaSW(:,1);

end

%%% Calculate T from Mg/Ca
function TempOut = calcTemp(sampleage,speciesID,foramDataWOAT,foramData,...
    co2sysOut,Omega,sampledata,TelID,MgCaSW,LiCaSW,Hin)

    % Convert uncertainty back to absolute value before transposing to
    % calculate uncertainty 
    tempOmega = [Omega(1,1) Omega(1,1)+Omega(1,2)].^-2;
    tempOmega = [tempOmega(:,1) tempOmega(:,2)-tempOmega(:,1)];

    % Only include unflagged analyses in calibration
    isFlag = isnan(foramData.flagged_mgca);

    %%% find relevant portion of the database
    if speciesID==1 % L. wuellerstorfi
        DataSubset = any((strncmp(foramData.species,"l.wue",5) & isFlag),2);
    elseif speciesID==2 % Uvigerina spp.
        DataSubset = any((strncmp(foramData.species,"u.per",5) | ...
            strncmp(foramData.species,"uvi",3) & isFlag),2);
    elseif speciesID==3 % O. umbonatus
        DataSubset = any((strncmp(foramData.species,"o.umb",5) & isFlag),2);
    elseif speciesID==4 % C. mundulus
        DataSubset = any((strncmp(foramData.species,"c.mun",5) & isFlag),2);
    elseif speciesID==5 % C. pachyderma
        DataSubset = any((strncmp(foramData.species,"c.pac",5) & isFlag),2);
    elseif speciesID==6 % N. umbonifera
        DataSubset = any((strncmp(foramData.species,"c.pac",5) & isFlag),2);
    end

    %%% Calculate coefficient covariance for relevant species
    % If C. pachyderma Mg/Ca, include Omega in regression
    if TelID==3 && (speciesID==1 || speciesID==4)
        tempx = [foramDataWOAT(DataSubset) co2sysOut(DataSubset,17).^-2];
    else
        tempx = foramDataWOAT(DataSubset);
    end
    if TelID==3 % Mg/Ca
        tempy = log(foramData.mgca(DataSubset));
    else % Mg/Li
        tempy = foramData.mgca(DataSubset)./foramData.lica(DataSubset);
    end

    % Combine X and y into a single matrix for bootstrapping
    data = [tempx, tempy];
    % best model estimate
    refFit = fitlm(tempx,tempy,'linear');
    % Perform bootstrapping
    bootCoefficients = bootstrp(1000, @fitlm_coefficients, data);

    if TelID==3 % If Mg/Ca
        % sample data 10000 times including uncertainty
        ElCaRand = repmat(sampledata(:,1),1,10000);
        % Add random measurement uncertainties
        ElCaRand = ElCaRand + normrnd(0,0.5,size(ElCaRand)).*...
            repmat(sampledata(:,2),1,size(ElCaRand,2));

    else % If Mg/Li

        % Treat Mg/Ca data (first two columns in the same way)
        % sample data 10000 times including uncertainty
        ElCaRand = repmat(sampledata(:,1),1,10000);
        % Add random measurement uncertainties
        ElCaRand = ElCaRand + normrnd(0,0.5,size(ElCaRand)).*...
            repmat(sampledata(:,2),1,size(ElCaRand,2));
        
        if size(sampledata,2)==4 % If Mg/Ca and Li/Ca fed in seperately
            % Add uncertainties to Li/Ca data too
            % sample data 10000 times including uncertainty
            ElCaRand2 = repmat(sampledata(:,3),1,10000);
            % Add random measurement uncertainties
            ElCaRand2 = ElCaRand2 + normrnd(0,0.5,size(ElCaRand2)).*...
                repmat(sampledata(:,4),1,size(ElCaRand2,2));
           
        end
    end
    
    % Add seawater Mg/Ca + uncertainties in MgCaSW and H, if processing
    % samples older than 0.8 Ma
    if max(sampleage)>800
        % If Mg/Ca data or seperate Mg/Ca and Li/Ca datasets
        if TelID==3 || (TelID==4 && size(sampledata,2)==4)
            % adjust sample data for MgCaSW change
            sampledata(:,1) = sampledata(:,1).*5.2^Hin(1)./MgCaSW(:,1).^Hin(1);
    
            % Randomly resample the value of H
            Hin = repmat(normrnd(Hin(1),Hin(2)/2,1,size(ElCaRand,2)),...
                size(ElCaRand,1),1);
            % Randomly resample Mg/CaSW
            MgCaSW = repmat(MgCaSW(:,3),1,size(ElCaRand,2)) + ...
                repmat(normrnd(0.5,0.25,1,size(ElCaRand,2)),size(ElCaRand,1),1).*...
                repmat(MgCaSW(:,2) - MgCaSW(:,3),1,size(ElCaRand,2));
            % Adjust resampled Mg/Ca data
            ElCaRand = ElCaRand.*5.2.^Hin./MgCaSW.^Hin;
        end

        % If Mg/Ca and Li/Ca fed in seperately, Li/Ca also needs correcting
        if TelID==4 && size(sampledata,2)==4

            % Reformat seawater Li/Ca data
            %LiCaSW = ones(size(sampledata,1),3).*[4.43/10.3 29.23/10.3 6.63/10.3];
            LiCaSW = [LiCaSW(:,2)-LiCaSW(:,1) LiCaSW(:,2) LiCaSW(:,3)-LiCaSW(:,2)];

            % adjust sample data for LiCaSW change
            sampledata(:,3) = sampledata(:,3).*(29.23/10.3)./LiCaSW(:,2);

            % Randomly resample Li
            LiRand = normrnd(0,0.5,1,size(ElCaRand2,2));
            LiCaSWhigh = repmat(LiCaSW(:,2),1,size(ElCaRand2,2)) + ...
                repmat(LiRand,size(ElCaRand2,1),1).*...
                repmat(LiCaSW(:,3),1,size(ElCaRand2,2));
            LiCaSWlow = repmat(LiCaSW(:,2),1,size(ElCaRand2,2)) - ...
                repmat(LiRand,size(ElCaRand2,1),1).*...
                repmat(LiCaSW(:,1),1,size(ElCaRand2,2));
            LiCaSW(:,LiRand>=0) = LiCaSWhigh(:,LiRand>=0);
            LiCaSW(:,LiRand<0) = LiCaSWlow(:,LiRand<0);
            % Adjust resampled Li/Ca data assuming linear seawater-shell Li/Ca
            ElCaRand2 = ElCaRand2.*(29.23/10.3)./LiCaSW;

        end
    end

    if size(sampledata,2)==4 % If Mg/Ca and Li/Ca fed in seperately
        % Combine Mg/Ca and Li/Ca datasets
        ElCaRand = ElCaRand./ElCaRand2;
    end
    
    % If processing Mg/Ca data
    % If L. wuellerstorfi or C. mundulus, Mg/Ca is sensitive to T & Omega
    if TelID==3 && (speciesID==1 || speciesID==4)
        % covariation is strongest between T and Omega coefficients rather 
        % than Omega and intercept in this case
        fitboostrap = ...
            fitlm(bootCoefficients(:,2),bootCoefficients(:,3)); % m = x1 + x2*int

        % sample model coefficients, use SD=1 as fitlm returns 1s
        % uncertainties
        coefRand = normrnd(0,1,[4 size(ElCaRand,2)]);
        % calculate all possible Omega
        tempTsens = refFit.Coefficients{2,1} + ... % Tsens = x1 + X2*Osens
            refFit.Coefficients{2,2}.*repmat(coefRand(1,:),size(ElCaRand,1),1);
        tempOsens = tempTsens.*(fitboostrap.Coefficients{2,1} + ...
            fitboostrap.Coefficients{2,2}.*repmat(coefRand(2,:),size(ElCaRand,1),1)) + ...
            (fitboostrap.Coefficients{1,1} + ...
            fitboostrap.Coefficients{1,2}.*repmat(coefRand(3,:),size(ElCaRand,1),1));

        %%% Temperature calculations
        % Randomly sample T within uncertainty
        tempOmegaresamp = tempOmega(:,1) + ...
            tempOmega(:,2).*normrnd(0,0.5,size(ElCaRand));
        % Temperature, with uncertainty
        Tout = (log(ElCaRand) - ...
            ((refFit.Coefficients{1,1} + refFit.Coefficients{1,2}.*...
            repmat(coefRand(4,:),size(ElCaRand,1),1)) + ...
            tempOmegaresamp.*tempOsens))./tempTsens;
        % Temperature, no uncertainty
        Refout = (log(sampledata(:,1)) - ...
            (refFit.Coefficients{1,1} + refFit.Coefficients{3,1}.*tempOmega(:,1)))./...
            refFit.Coefficients{2,1};
        % Calculate inverse prediction intervals
        tempY = (tempy - ...
            (refFit.Coefficients{1,1} + refFit.Coefficients{3,1}.*...
            co2sysOut(DataSubset,17).^-2))./...
            refFit.Coefficients{2,1};
        tempX = tempx(:,2);
        IPItemp = calculate_IPI(0,[tempX(~isnan(tempY)) tempY(~isnan(tempY))]);
    else % If Mg/Li or Mg/Ca in the case of species without an Omega sensitivity
        fitboostrap = ...
            fitlm(bootCoefficients(:,1),bootCoefficients(:,2)); % m = x1 + x2*int

        % sample model coefficients
        coefRand = normrnd(0,1,[3 size(ElCaRand,2)]);
        % calculate all possible Temperature
        tempInt = refFit.Coefficients{1,1} + ... % Temp. coef = x1 + X2*B
            refFit.Coefficients{1,2}.*repmat(coefRand(1,:),size(ElCaRand,1),1);
        tempSens = tempInt.*(fitboostrap.Coefficients{2,1} + ...
            fitboostrap.Coefficients{2,2}.*repmat(coefRand(2,:),size(ElCaRand,1),1)) + ...
            (fitboostrap.Coefficients{1,1} + ...
            fitboostrap.Coefficients{1,2}.*repmat(coefRand(3,:),size(ElCaRand,1),1));

        %%% Temperature calculations
        if TelID==3
            % Mg/Ca temperature, with uncertainty
            Tout = (log(ElCaRand) - tempInt)./tempSens;
            % Temperature, no uncertainty
            Refout = (log(sampledata(:,1)) - refFit.Coefficients{1,1})./...
                refFit.Coefficients{2,1};
        else
            % Mg/Li temperature, with uncertainty
            Tout = (ElCaRand - tempInt)./tempSens;
            % Temperature, no uncertainty
            if size(sampledata,2)==4
                Refout = (sampledata(:,1)./sampledata(:,3) - refFit.Coefficients{1,1})./...
                    refFit.Coefficients{2,1};
            else
                Refout = (sampledata(:,1) - refFit.Coefficients{1,1})./...
                    refFit.Coefficients{2,1};
            end
        end
        % Calculate inverse prediction intervals
        tempY = (tempy - ...
            refFit.Coefficients{1,1})./refFit.Coefficients{2,1};
        tempX = tempx;
        IPItemp = calculate_IPI(0,[tempX(~isnan(tempY)) tempY(~isnan(tempY))]);
    end

    % Avoids negative root issues
    Tout = real(Tout);
    
    % Extract output data
    TempOut = [Refout-IPItemp ...
        prctile(Tout,2.5,2) prctile(Tout,16,2) ...
        Refout prctile(Tout,84,2) ...
        prctile(Tout,97.5,2) Refout+IPItemp];

end