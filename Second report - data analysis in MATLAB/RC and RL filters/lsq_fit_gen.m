function [fit_out,dfit_out,C,chi2,N_DOF]=lsq_fit_gen(y,F,varargin)
% [fit_out,dfit_out,C,chi2,N_DOF]=lsq_fit_gen(y,F,varargin)
%
% general LSQ fit routine that fits data x(t) as a sum of functions given
% in the matrix F:
%       model: y = SUM_j(a_j * f_j)
%       found by minimizing SUM_i((y_i - SUM_j(a_j * f_j(i)).^2 / sigma_i^2)
%       size(F) = (N_pts x N_func), where N_pts = length(y) and N_func is the 
%       number of functions that we are fitting to.  
%       y is the column vector with the dependent variable data. 
%
% output:   fitout (fit parameters): [a_1; a_2; ... ...a_N_func]
%
%           dfit_out (fit uncertainties):   [da_ord1; da_ord2; ... ]
%           C   full covariance matrix
%           chi2    fit chi^2 (meaningless if fit errors not given)
%           N_DOF number of degrees of freedom
%
% options:
% 'err'     specify data uncertainties (next arg either a vector of errors
%           or a single number to be applied to all points ... otherwise
%           unity errors and a good fit are assumed)
% 'x'       allows option of plotting against some independent variable
%           (such as time).   next argument is column vector x (length(x) = length(y))
%
% 'logplot' asks for logaritmic plots
% 'nopl', 'nobs' turn off the plotting and text summary outputs
%
% bw 27/07/2010

plt=0;
talk=1;
err_given=0; % assume errors are NOT specified, and
% thus we assign unity error to each point
first_time=1; % internal variable to turn off certain things the second time round
% when evaluating the x errors
last_time=1;
sigtest = 0;
logplot = 0;

if ~isempty(varargin)
    j = 1;
    while j<=length(varargin)
        switch lower(repr(varargin{j}))
            case 'err'
                err = varargin{j+1};
                err_given = 1;
                j = j+1;
            case 'x'
                plt=1;
                x = varargin{j+1};
                if length(x) == length(y)
                    plt = 1;
                    j = j+1;
                end
            case 'nobs'
                talk=0;
            case 'one_more_time'
                first_time=0;
            case 'sigtest'
                sigtest = 1;
            case 'logplot'
                logplot = 1;
            otherwise
                warning('Unknown option "%s"',repr(varargin{j}));
        end
        j=j+1;
    end
end


if ~exist('err')
    err=ones(length(y),1);
    if talk
        disp(' Assuming all points have unity error');
    end
elseif length(err)==1
    if talk
        disp([' Assuming all points have error dy = ' num2str(err)]);
    end
    err=err*ones(length(y),1);
else
    if talk
        disp(' Errors entered point by point');
    end
end


M=size(F,2);
N_func = M;
N_pts = size(F,1);

G=[];
V=[]; % data vector
for i=1:M
    V(i) = sum (F(:,i) .* y ./ err.^2);
    for j=1:M
        G(i,j)=sum(F(:,i) .* F(:,j) ./ err.^2);
    end
end
V=V'; % it automatically makes V a row vector, we want a column vector

% "function" matrix
% "data" vector:

C=inv(G);
fit_out= C * V;

y_fit= F*fit_out;
dy_res=y-y_fit;
dy_mean=sum(dy_res)/(length(y)-length(fit_out));

% assign fit uncertainties, if necessary assuming a good fit (chi^2 = 1)
N_DOF=length(y) - length(fit_out);
chi2=sum(dy_res.^2 ./ err.^2) / N_DOF;

if (err_given==0)
    % here we assume a good fit, and scale the covariance matrix such that chi^2 would be = 1
    % this allows us to get an idea of the fit errors
    C=C*chi2;
end

dfit_out=[];
for i=1:length(fit_out)
    dfit_out(i)=sqrt(C(i,i));
end
dfit_out=dfit_out';  % make it a column vector


sigma=sqrt(sum(dy_res.^2)/N_DOF);

if talk && first_time
    disp([' Fitting ' int2str(N_pts) ' data points to model with ' ... 
        int2str(N_func) ' functions']);
end

if talk
    disp('   ');
    disp(' Fit results: ');
    %
    if err_given
        disp([' Fit chi^2 = ' num2str(chi2) ' per DOF (' int2str(N_DOF) ' DOF)']);
    else
        disp(' Assume good fit for fit parameter uncertainties');
    end
    disp([' Average fit residual = ' num2str(dy_mean)]);
    disp([' RMS residual = ' num2str(sigma)]);
    disp(' Fit parameters: ');
    if N_func >=0
        for j=1:N_func
            disp([' Coefficient ' int2str(j) ': ' num2str(fit_out(j)) ' +/- ' ...
                num2str(dfit_out(j))]);
        end
    end
    disp('   ');
end

if plt && last_time
    figure
    axtop=axes('position',[.15 .5 .75 .4],'fontsize',14);
    if err_given
        errorbar(x,y,err,'b.');
    else
        plot(x,y,'b.');
    end
    if logplot
        set(gca,'XScale','log')
        set(gca,'YScale','log')
    end

    x_mod=x;
    y_mod=F * fit_out;
    hold on

    plot (x_mod,y_mod,'r');
    grid on;
    axbot=axes('position',[.15 .1 .75 .3],'fontsize',14);
    plot (x,dy_res,'r');
    if err_given
        errorbar(x,dy_res,err,'b.');
    else
        plot(x,dy_res,'b.');
    end
    if logplot
        set(gca,'XScale','log')
        set(gca,'YScale','log')
    end
    grid on;
end

if talk
    disp('   ');
end


