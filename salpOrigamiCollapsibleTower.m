%Generates figures used for lasercutting origami collapsing cylinder
%Save output figures as .svg then pop over to a cutter
%Nathan Justus - Jan 30, 2024

clear all

%% Define design parameters and cutfile storage location

%Polygon side length
L = 1;
%Polygon number of sides
N = 4;
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
%Make folder to save figures
mkdir('origamiFigures',['N',num2str(N),'M',num2str(M)]);
%File preamble
pbl = ['origamiFigures\',['N',num2str(N),'M',num2str(M)],'\'];
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
    plot(these_xs,ys,'k','LineWidth',lw);
end

%Plot the diagonal bits
for i = 2:2:2*M
    for j = 1:(numel(ys_v1)-1)
        plot([xs(i),xs(i+1)],[ys_v1(j),ys_v2(j+1)],'k','LineWidth',lw);
        plot([xs(i+1),xs(i+2)],[ys_v2(j+1),ys_v1(j)],'k','LineWidth',lw);
    end
end

%Plot the bits used to make the structural rings at each end
xs_stripe = [L/3-L/20,2*L/3-L/15];
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
fig2 = figure(2);
clf;
hold on;
set(fig2,'Units','inches');

%Bottom edge
plot(xs,ys_h,'r','LineWidth',lw);
%Top edge
plot(xs,ys_h+N*L,'r','LineWidth',lw);
%Left edge
plot(xs(1)*ones(size(ys_v1)),ys_v1,'r','LineWidth',lw);
%Right edge
plot(xs(end)*ones(size(ys_v1)),ys_v1,'r','LineWidth',lw);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([min(xs),max(xs),min(ys_v1),max(ys_v2)]);
fig2.Position = [1,1,max(xs)-min(xs),max(ys_v2)-min(ys_v1)];
saveas(fig2,[pbl,'CylinderCut',fend]);

%% Generate the etch pattern for the nozzle

%Figure prep
fig3 = figure(3);
clf;
hold on;
set(fig3,'Units','inches');

%Radius of polygon inscribed circle
R = L/(2*tan(pi/N));
%Length from center of a polygon edge to the tip of nozzle cylinder
%(pretending that it extends all the way to the thruster central axis)
nozzleL = R/cos(nozzleTheta);

%Angle to rotate around for rolled-out nozzle
dtheta = 2*atan2(L/2,nozzleL);
thetas = [0:dtheta:N*dtheta];

%Like nozzleL but for a polygon vertex instead of edge center
nozzleH = norm([nozzleL,L/2]);
%Add a little bit to make the tags to attach the nozzle to the cylinder
tagL = nozzleH + L/3;

%Exit area ratio is squared linear ratio of polygon sides
linearRatio = sqrt(A_exit);
nozzleR = linearRatio*nozzleH;

%Save points for relevant rings for the rolled-out nozzle
inner = zeros(2,N+1);
outer = zeros(2,N+1);
flap = zeros(2,N+1);
for i = 1:numel(thetas)
    inner(1,i) = nozzleR*cos(thetas(i));
    inner(2,i) = nozzleR*sin(thetas(i));
    outer(1,i) = nozzleH*cos(thetas(i));
    outer(2,i) = nozzleH*sin(thetas(i));
    flap(1,i) = tagL*cos(thetas(i));
    flap(2,i) = tagL*sin(thetas(i));
end

%Plot the inner and outer rings
for i = 1:N
    plot(inner(1,i:(i+1)),inner(2,i:(i+1)),'k','LineWidth',lw);
    plot(outer(1,i:(i+1)),outer(2,i:(i+1)),'k','LineWidth',lw);
end
%Plot the radial lines
for i = 1:N
    plot([inner(1,i),flap(1,i)],[inner(2,i),flap(2,i)],'k','LineWidth',lw);
end
if mod(N,2) == 0
    plot([inner(1,N+1),outer(1,N+1)],[inner(2,N+1),outer(2,N+1)],'k','LineWidth',lw);
else
    plot([inner(1,N+1),flap(1,N+1)],[inner(2,N+1),flap(2,N+1)],'k','LineWidth',lw);
end
%Plot the flap outer edges
for i = 1:2:N
    plot(flap(1,i:(i+1)),flap(2,i:(i+1)),'k','LineWidth',lw);
end

lows = min([inner,outer,flap]')';
highs = max([inner,outer,flap]')';

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([lows(1),highs(1),lows(2),highs(2)]);
fig3.Position = [1,1,highs(1)-lows(1),highs(2)-lows(2)];
saveas(fig3,[pbl,'NozzleEtch',fend]);

%% Generate nozzle cut pattern outline

%Can do a separate figure or just plot a different color on the same plot
fig4 = figure(4);
clf;
hold on;
set(fig4,'Units','inches');

nFlapUnits = floor(N/2);

edgePoints = [];
edgePoints = [edgePoints,fliplr(inner)];

for i = 1:nFlapUnits
    index = 2*i-1;
    edgePoints = [edgePoints,flap(:,index),flap(:,index+1),outer(:,index+1),outer(:,index+2)];
end

if mod(N,2)
    edgePoints = [edgePoints,flap(:,end-1),flap(:,end)];
end

edgePoints = [edgePoints,inner(:,end)];

plot(edgePoints(1,:),edgePoints(2,:),'r','LineWidth',lw);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([lows(1),highs(1),lows(2),highs(2)]);
fig4.Position = [1,1,highs(1)-lows(1),highs(2)-lows(2)];
saveas(fig4,[pbl,'NozzleCut',fend]);

%% Plot Checkvalve Cap Etch

fig5 = figure(5);
clf;
hold on;
set(fig5,'Units','inches');

H1 = L/(2*sin(pi/N));
H2 = H1 + L/3;

R = 7*H1/12;

thetas = linspace(0,2*pi,N+1);

edgePoints = [];
outEdgePoints = [];

for i = 1:numel(thetas)
    edgePoints = [edgePoints,[H1*cos(thetas(i));H1*sin(thetas(i))]];
    outEdgePoints = [outEdgePoints,[H2*cos(thetas(i));H2*sin(thetas(i))]];
end

for i = 1:floor(N/2)
    index = 2*i-1;
    edgePoints = [edgePoints,outEdgePoints(:,index),outEdgePoints(:,index+1),edgePoints(:,index+1),edgePoints(:,index+2)];
end

if mod(N,2)
    edgePoints = [edgePoints,outEdgePoints(:,end-1),outEdgePoints(:,end)];
end

thetas_circle = linspace(0,2*pi,101);

plot(edgePoints(1,:),edgePoints(2,:),'k','LineWidth',lw);
plot(R*cos(thetas_circle),R*sin(thetas_circle),'k','LineWidth',lw);

lows = min(edgePoints')';
highs = max(edgePoints')';

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([lows(1),highs(1),lows(2),highs(2)]);
fig5.Position = [1,1,highs(1)-lows(1),highs(2)-lows(2)];
saveas(fig5,[pbl,'CheckValveEtch',fend]);

%% Plot Checkvalve Cap Cut

fig6 = figure(6);
clf;
hold on;
set(fig6,'Units','inches');
xticks();
yticks();
xticklabels();
yticklabels();
axis equal;
axis off;

H1 = L/(2*sin(pi/N));
H2 = H1 + L/3;

R = 7*H1/12;

thetas = linspace(0,2*pi,N+1);

edgePoints = [];
outEdgePoints = [];

for i = 1:numel(thetas)
    edgePoints = [edgePoints,[H1*cos(thetas(i));H1*sin(thetas(i))]];
    outEdgePoints = [outEdgePoints,[H2*cos(thetas(i));H2*sin(thetas(i))]];
end

for i = 1:floor(N/2)
    index = 2*i-1;
    edgePoints = [edgePoints,outEdgePoints(:,index),outEdgePoints(:,index+1),edgePoints(:,index+1),edgePoints(:,index+2)];
end

if mod(N,2)
    edgePoints = [edgePoints,outEdgePoints(:,end-1),outEdgePoints(:,end)];
end

thetas_circle = linspace(0,2*pi,101);

plot(edgePoints(1,N+1:end),edgePoints(2,N+1:end),'r','LineWidth',lw);
plot(R*cos(thetas_circle),R*sin(thetas_circle),'r','LineWidth',lw);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([lows(1),highs(1),lows(2),highs(2)]);
fig6.Position = [1,1,highs(1)-lows(1),highs(2)-lows(2)];
saveas(fig6,[pbl,'CheckValveCut',fend]);

%% Plot Valve

fig7 = figure(7);
clf;
hold on;
set(fig7,'Units','inches');
xticks();
yticks();
xticklabels();
yticklabels();
axis equal;
axis off;

H1 = L/(2*sin(pi/N));
H2 = H1 + L/3;

R = 7*H1/12;

ValveR = (H1+R)/2;

thetas = linspace(0,2*pi,N+1);

edgePoints = [];

for i = 1:numel(thetas)
    edgePoints = [edgePoints,[ValveR*cos(thetas(i));ValveR*sin(thetas(i))]];
end

plot(edgePoints(1,:),edgePoints(2,:),'r','LineWidth',lw);

lows = min(edgePoints')';
highs = max(edgePoints')';

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([lows(1),highs(1),lows(2),highs(2)]);
fig7.Position = [1,1,highs(1)-lows(1),highs(2)-lows(2)];
saveas(fig7,[pbl,'ValveCut',fend]);

if closeTheFigures
    close all
end

