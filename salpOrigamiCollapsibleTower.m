%Generates figures used for lasercutting origami collapsing cylinder
%Save output figures as .svg then pop over to a cutter
%Nathan Justus - Jan 30, 2024

clear all

%% Define design parameters and cutfile storage location

%Polygon side length
L = 7/8;
%Polygon number of sides
N = 6;
%Number of collapse sections
M = 3;
%Nozzle Angle (larger angle looks like a cylinder, smaller like a plane)
nozzleTheta = pi/4;
%Exit area ratio (ratio of nozzle area to cylinder area)
A_exit = 0.1;

%Rhombus angle
theta = pi/N;
theta = pi/2-(2*theta);

%Line width for plotting
lw = 1;
%Check if figures folder exists, make it if not
if ~exist('origamiFigures','dir')
    mkdir('origamiFigures');
end

Lstr = '0p875';
%Make folder to save figures
mkdir('origamiFigures',['N',num2str(N),'M',num2str(M),'L',Lstr]);
%File preamble
pbl = ['origamiFigures\',['N',num2str(N),'M',num2str(M),'L',Lstr],'\'];
%File ending
fend = ['_N',num2str(N),'M',num2str(M),'.svg'];
%Boolean to close figures after generating and saving
closeTheFigures = false;

%% Generate etch-pattern for main cylinder

%Get x coordinates for rows in the etch
midBits = L:L*cos(theta):2*M*L*cos(theta)+L;
xs = [0,midBits,2*M*L*cos(theta)+2*L];
%y values for a vertical stack of nodes
ys_v1 = [0:L:N*L];
%y values shifted up for the alternating bits
ys_v2 = ys_v1+L*sin(theta);

%Generate y values for a horizontal row of nodes
stepUp = [0,L*sin(theta)];
steps = repmat(stepUp,[1,M]);
ys_h = [0,steps,0,0];

%Set up figure for plotting
fig1 = figure(1);
clf;
hold on;
set(fig1,'Units','inches');

%Plot all the rows
for i = 0:L:N*L
    plot(xs,ys_h+i,'k','LineWidth',lw);
end

%Plot all the columns
for i = 1:numel(xs)
    if ys_h(i)
        ys = ys_v2;
    else
        ys = ys_v1;
    end
    these_xs = ones(N+1,1)*xs(i);
    if i == numel(xs)
        dottedLine = plot(these_xs,ys,'k','LineWidth',lw);
    else
        plot(these_xs,ys,'k','LineWidth',lw);
    end
end

%Plot the diagonal bits
for i = 2:2:2*M
    for j = 1:(numel(ys_v1)-1)
        plot([xs(i),xs(i+1)],[ys_v1(j),ys_v2(j+1)],'k','LineWidth',lw);
        if i == 2*M
            dashedLine = plot([xs(i+1),xs(i+2)],[ys_v2(j+1),ys_v1(j)],'k','LineWidth',lw);
        else
            plot([xs(i+1),xs(i+2)],[ys_v2(j+1),ys_v1(j)],'k','LineWidth',lw);
        end
    end
end

%Plot the bits used to make the structural rings at each end
xs_stripe = [L/2];
xs_stripe = [xs_stripe,2*M*L*cos(theta)+2*L-fliplr(xs_stripe)];
for i = 1:numel(xs_stripe)
    these_xs = ones(size(ys_v1))*xs_stripe(i);
    plot(these_xs,ys_v1,'k','LineWidth',lw);
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([min(xs),max(xs),min(ys_v1),max(ys_v2)]);
fig1.Position = [1,1,max(xs)-min(xs),max(ys_v2)-min(ys_v1)];
saveas(fig1,[pbl,'CylinderEtch',fend]);

%% Generate the cut pattern outline for the main cylinder

%Can do a separate figure or just plot a different color on the same plot
% fig2 = figure(2);
% clf;
% hold on;
% set(fig2,'Units','inches');

%Bottom edge
plot(xs,ys_h,'r','LineWidth',lw);
%Top edge
plot(xs,ys_h+N*L,'r','LineWidth',lw);
%Left edge
plot(xs(1)*ones(size(ys_v1)),ys_v1,'r','LineWidth',lw);
%Right edge
redLine = plot(xs(end)*ones(size(ys_v1)),ys_v1,'r','LineWidth',lw);

holeR = .175/2;
thetas = linspace(0,2*pi,50);
holeX = holeR*cos(thetas);
holeY = holeR*sin(thetas);

xHoles = [L/4,3*L/4];
xHoles = [xHoles,2*M*L*cos(theta)+2*L-fliplr(xHoles)];
for i = 1:numel(xHoles)
    for j = 1:N
        thisX = xHoles(i);
        thisY = (j-1)*L + L/2;
        plot(holeX+thisX,holeY+thisY,'r','LineWidth',lw);
    end
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([min(xs),max(xs),min(ys_v1),max(ys_v2)]);
fig2.Position = [1,1,max(xs)-min(xs),max(ys_v2)-min(ys_v1)];
saveas(fig2,[pbl,'CylinderCut_Combined',fend]);
