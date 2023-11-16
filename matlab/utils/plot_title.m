function plot_title(title_str,varargin)
%
% Plot a title
%
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:50, h{i}.first_flag = true; end; end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'tfs',15);               % title font size
addParameter(ps,'tfn','consolas');       % title font name 
addParameter(ps,'tfc','k');              % title font color
addParameter(ps,'interpreter','latex');  % 'none', 'latex'
addParameter(ps,'ALLOW_UNDERBAR',0);     % allow $a_1$
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
tfs             = ps.Results.tfs;
tfn             = ps.Results.tfn;
tfc             = ps.Results.tfc;
interpreter     = ps.Results.interpreter;
ALLOW_UNDERBAR  = ps.Results.ALLOW_UNDERBAR;

% Plot title 
if isequal(interpreter,'latex') && (ALLOW_UNDERBAR==0)
    % change '_' to '-' for latex formatting
    title_str =  strrep(title_str,'_','-'); 
end
if h{fig_idx}.first_flag || ~ishandle(h{fig_idx}.fig)
    h{fig_idx}.first_flag = false;
    h{fig_idx}.fig = figure(fig_idx);
    h{fig_idx}.title = title(title_str,...
        'fontsize',tfs,'fontname',tfn,'interpreter',interpreter,'color',tfc);
else
    h{fig_idx}.title.String = title_str;
    h{fig_idx}.title.Color  = tfc;
end
