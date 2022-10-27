function ch = smallcolorbar(axs,location,axiscl)

% ch = smallcolorbar(axs,location,axiscl)

if nargin < 1; axs      = gca;           end
if nargin < 2; location = 'eastoutside'; end
if nargin < 3; axiscl   = [0 0 0];       end

axspos = get(axs, 'position');
ratio  = [.05 .33];

switch location
  case 'eastoutside'
      x_size = axspos(3)*ratio(1);
      y_size = axspos(4)*ratio(2);
      pos    = [sum(axspos([1 3])) + x_size*.15, axspos(2), x_size, y_size];
    case 'northeastoutside'
        x_size = axspos(3)*ratio(1);
        y_size = (axspos(4)*ratio(2));
        pos    = [sum(axspos([1 3])) + (x_size*.15), .8-axspos(2)-y_size, x_size, y_size];
        location = 'eastoutside';
    case 'southoutside'
        x_size = axspos(3)*ratio(2);
        y_size = axspos(4)*ratio(1);
        pos    = [axspos(1), axspos(2) - y_size*1.15, x_size, y_size];
    
  otherwise
    error('position not recognized')
end

ch = colorbar('location',location,'position',pos,'Color',axiscl);