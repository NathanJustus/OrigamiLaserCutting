%Generates figures used for lasercutting origami collapsing sphere
%Save output figures as .svg then pop over to a cutter
%Nathan Justus - Feb 1, 2024

clear all

%% Define design parameters and cutfile storage location

%Unit length
L = 1;
%Number of rows
R = 6;
%Number of columns
C = 9;


%Line width for plotting
lw = 1;
%Check if figures folder exists, make it if not
if ~exist('origamiFigures','dir')
    mkdir('origamiFigures');
end
%Make folder to save figures
mkdir('origamiFigures',['Sphere_','R',num2str(R),'C',num2str(C)]);
%File preamble
pbl = ['origamiFigures\','Sphere_','R',num2str(R),'C',num2str(C),'\'];
%File ending
fend = ['_R',num2str(R),'C',num2str(C),'.svg'];
%Boolean to close figures after generating and saving
closeTheFigures = false;

%% Generate etch-pattern for main cylinder

%Set up figure for plotting
fig1 = figure(1);
clf;
hold on;
set(fig1,'Units','inches');

for row = 1:R
    for col = 1:C

        xL = (col-1)*L;
        xm = xL + L/2;
        xR = xL + L;
        yB = (row-1)*L;
        ym = yB + L/2;
        yT = yB + L;

        if mod(row,2)
            xs = [xL,xm,xR,xm,xL];
            ys = [ym,yT,ym,yB,ym];
            plot(xs,ys,'k','LineWidth',lw);
        else
            x1s = [xL,xm,xR];
            y1s = [yT,ym,yT];
            x2s = [xL,xm,xR];
            y2s = [yB,ym,yB];
            plot(x1s,y1s,'k','LineWidth',lw);
            plot(x2s,y2s,'k','LineWidth',lw);
        end

    end
end

for i = 1:(R+1)
    y = L*(i-1);
    plot([0,C*L],[y,y],'k','LineWidth',lw);
end
for j = 1:.5:(C+1)
    x = L*(j-1);
    plot([x,x],[0,R*L],'k','LineWidth',lw);
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([0,C*L,0,R*L]);
fig1.Position = [1,1,C*L,R*L];
saveas(fig1,[pbl,'SphereEtch',fend]);

%% Generate the cut pattern outline for the sphere
%Can do a separate figure or just plot a different color on the same plot
fig2 = figure(2);
clf;
hold on;
set(fig2,'Units','inches');

%Draw the perimeter box
xs = [0,C*L,C*L,0,0];
ys = [0,0,R*L,R*L,0];
plot(xs,ys,'r','LineWidth',lw);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis equal;
axis off;

axis([0,C*L,0,R*L]);
fig2.Position = [1,1,C*L,R*L];
saveas(fig2,[pbl,'SphereCut',fend]);

if closeTheFigures
    close all
end

