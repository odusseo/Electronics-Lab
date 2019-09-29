% A few useful command-line examples to print a figure into a file

% M Hueller 15/04/2014


%% Prepare a dummy figure

% Make sure we capture the figure handle for later use
hfig = figure();

% Plot some simple data
plot([1:10], [2:2:20], 'rx', 'MarkerSize', 12);

% Add labels
xlabel('Frequency [Hz]');
ylabel('Output [V]');

% Add legend
legend('My nice data');

% Grid
grid on;

%% Paper settings
% Set paper orientation
set(hfig, 'PaperOrientation', 'landscape');

% Set paper type
set(hfig, 'PaperType', 'A4');

% Set paper positioning units
set(hfig, 'PaperUnits', 'centimeters');

% Set paper positioning to be manual
set(hfig, 'PaperPositionMode', 'manual');

% Set paper positioning to be filling better the page
set(hfig, 'PaperPosition', [0 0 29 21]);


%% Printing to file

% exports a figure to an EPS color format file, myfigure01.eps, and includes a color TIFF preview.
print(hfig, '-depsc', '-tiff', 'myfigure01.eps')

% exports a figure to an EPS color format file, myfigure02.eps
print(hfig, '-depsc', 'myfigure02.eps')

% exports a figure to an EPS black-and-white format file, myfigure03.eps, and includes a color TIFF preview.
print(hfig, '-depsc', '-tiff', 'myfigure03.eps')

% exports a figure to an EPS black-and-white format file, myfigure04.eps
print(hfig, '-depsc', 'myfigure04.eps')

% exports a figure to a PDF color format file, myfigure05.pdf
print(hfig, '-dpdf', 'myfigure05.pdf')

% exports a figure to a PNG format file, myfigure06.png
print(hfig, '-dpng', 'myfigure06.png')

% Version info:
% $Id: printing.m 4736 2014-04-15 16:53:39Z mauro.hueller $
