function save_fig(F,PATHOUT, NAME)
% save_fig: Function to save figure with specific parameters
%     This function helps to save a certain figure.
%     
% INPUT: 
%     F: any figure handle e.g.: gcf)
%     PATHOUT: string to where figure should be saved
%     NAME: Name as string under which figure should be saved

% OUTPUT:
%     No output, as figure is saved in specified location
% 
% Figure format is .png as default, can be changed to any type
% 
% 
% Author: (Julius Welzel & Mareike Daeglau, University of Oldenburg, 2018)

% check inputs
if ~ischar(PATHOUT) | ~ischar(PATHOUT)
    error('Specifiy PATHOUT and NAME as string')
end

% set figure options
%set(F, 'Units','normalized','outerposition',[0 0 1 1]); % full screen
%set(F,'color','white'); % set figure background color in rgb mode
%print('-bestfit','-dpdf')
saveas(F,[PATHOUT NAME '.png']); % save figure 
saveas(F,[PATHOUT NAME '.svg']); 

set(F,'Units','Inches');
pos = get(F,'Position');
set(F,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(F,[PATHOUT NAME],'-dpdf','-r0')

close;


end

