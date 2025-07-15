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
