function [IPIout,regressionCoeff] = calculate_IPI(plotQ,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate simple IPI (Method 1 in McClelland et al. 2021)
%
% plotQ (1/0) will produce a graph of the data with IPI plotted
% The second input can either be:
%   not entered (function brings up a user interface to select data file)
%   a filename in the current directory ('my_data.csv')
%   an n*2 matrix (if the data are already loaded into the workspace)
%
% If the data are loaded through a user interface or via a filename
% this data must be in the form of a csv file with the predictor values
% in the first column and the measured values in the second column,
% with one header row
%
% Examples: 
% Calculate IPI & OLS regression coefficients, and produce a graph
% after manually selecting the data file:
% [IPI,regressionCoeff] = calculate_IPI(1)
% Read data from an exisiting variable in the workspace:
% [IPI,regressionCoeff] = calculate_IPI(1,calData)
% Load data using the filename:
% [IPI,regressionCoeff] = calculate_IPI(1,'my_calibration_data.csv')


% functions to calculate critical t values manually to avoid use of statistics toolbox
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));  % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;   % 1-tailed t-distribution
t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ‘lpha & Degrees-Of-Freedom v
Clev = 0.95;    % define confidence interval
Cqts = (1-Clev)/2;

if nargin==1    % if no input file is specified
    getfid = uigetfile('*.*');
    xydata = csvread(getfid,1,0);   % read data assuming top row is header
else
    if isa(varargin{:},'double')    % if a matrix is passed to the function
        xydata = varargin{:};
    else                            % load the filename pased to the function
    getfid = varargin{:};
    xydata = csvread(getfid,1,0);
    end
end
degF = size(xydata,1) - 2; %Degrees of freedom


B10hat = [ones(size(xydata,1),1) xydata(:,1)]\xydata(:,2);  % OLS regression coefficients

critT = t_inv(Cqts, degF) * [-1 1];    % critical t value

Yhat = B10hat(1) + xydata(:,1).*B10hat(2);
ysd = std(xydata(:,2) - Yhat);
Xel = ysd*critT(2)/B10hat(2);

IPIout = Xel;
regressionCoeff = B10hat;

X_new_hat = (xydata(:,2) - B10hat(1))./B10hat(2);

%display(['the inverse prediction interval of this dataset is ',num2str(round(Xel,2))])

if plotQ==1
    xPlotVals = round(min(xydata(:,1)),1):0.1:round(max(xydata(:,1)),1);
    close(figure(1))
    figure(1)
    hold on
    p0 = scatter(xydata(:,1),X_new_hat,'.');
    p1 = plot(xPlotVals,xPlotVals,'-k');
    p2 = plot(xPlotVals,xPlotVals+IPIout,'--k');
    plot(xPlotVals,xPlotVals-IPIout,'--k')
    xlabel('measured parameter')
    ylabel('predicted parameter')
    set(gcf,'color','w')
    set(gca,'box','on')
    legend([p0 p1 p2],'data','1:1 line','IPI','location','northwest')
    axislimits = [xlim ylim];
    text((axislimits(2)-axislimits(1))*0.8+axislimits(1),...
        (axislimits(4)-axislimits(3))*0.05+axislimits(3),['IPI = ',num2str(round(Xel,2))],'fontsize',8)
end
