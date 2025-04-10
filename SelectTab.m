function [h] = SelectTab()
%SelectTab Selects tab in current figure with multiple panes
%   Input arguments:
%       []
%   Output arguments:
%       (1) h: handle to new figure

h = gcf;
Pos = get(h,'Position');
ax = h.Children.SelectedTab.Children;
copyobj(ax,figure())
h = gcf;
set(h,'Position',Pos);
end