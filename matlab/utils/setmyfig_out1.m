function setmyfig_out1(fig)

%Default parameters
width  = 13.0;    % Width in inches
height = 12.0;    % Height in inches
alw    = 0.75;   % AxesLineWidth
fsz    = 20;     % Fontsize

% Set figures
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); 
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(fig, 'color', 'w');

% Set Tick Marks
% set(gca,'XTick',0:0.03:0.12);

% Prepare output for printing
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width+2, height];
set(gcf,'PaperPosition', myfiguresize);
set(gca, 'FontSize', fsz, 'LineWidth', alw);

