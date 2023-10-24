% Create figure
f = figure(1);
clf;

% Create plot axes
a = axes(f);

% scale plot axes
xlim([-0.25 1.25]);
ylim([-0.25 1.25]);

% Create vectors
createVector([0 1],[0 0],'Headwidth',5,'Headlength',5);
createVector([0 0.5*cos(pi/3)],[0 0.5*sin(pi/3)],'Headwidth',5,'Headlength',5);
createVector([0.5*cos(pi/3) 1],[0.5*sin(pi/3) 0],'Headwidth',5,'Headlength',5);

% Add text labels to plot
text(0.05,0.05,'\theta');
text(0.5,-0.05,'u')
text(0.25*cos(pi/3) - sqrt(0.05^2/2), 0.25*sin(pi/3) + sqrt(0.05^2/2),'v')
text(0.5*cos(pi/3) + (1 - 0.5*cos(pi/3))/2, 0.25*sin(pi/3) + sqrt(0.05^2/2), 'u - v')

% Function creates a vector
% (note vector is not scaled when axes changes sizes)
function createVector(x,y,varargin)

    % get position of current axes
    a = gca;
    p = a.Position;

    % determine normalized position of x coordinates
    xl = xlim;
    xs = p(3)/(xl(2) - xl(1));
    x = p(1) + xs*(x - xl(1));

    % determine normalized position of y coordinates
    yl = ylim;
    ys = p(4)/(yl(2) - yl(1));
    y = p(2) + ys*(y - yl(1));

    % create arrow
    annotation('arrow',x,y,varargin{:});
end