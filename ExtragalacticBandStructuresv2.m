%Extragalactic Band Structure Research Code
%Version 2
%Adapted by Alexander Criswell from ExponentialGalacticDisks.m by Dr.
%Curtis Struck and Dr. Bruce Elmegreen

%To set up particles in a disk on circular orbits in a galactic disk/halo
%potential, and perturb them with orbiting disturbers (Dwarf Galaxies)
%in a band structure around the disk. Additionally, output disk snapshots
%at later times to watch galactic evolution.

clear,clc
numrun = 4;

%First, set up disk using a variation of the SPH setup code.

%General Parameters
%UNITS:
%1 Mass Unit = 2.762810^8 Mo
%1 Radius Unit = .5 kpc
%1 Time Unit = 9.8 Myr
 
mcl = .007; %Particle Mass
RMAX = 9.0; %Galactic Radius
%Galactic Mass
GMtot = 0.0;
    Gcon = 4*pi*pi;   % G x M
    GM(1) = 0.99*Gcon;
    GMtot = GMtot + GM(1);
    
%Satellite Dwarf Galaxies
Ndist = 20; %Number of disturbers
Mdist = 100; %Disturber mass (multiple of base mass)
DZeros = zeros(1,Ndist+1);
GM = [GM DZeros];
itheta = zeros(1,Ndist+1);
iphi = zeros(1,Ndist+1);
ivtan = zeros(1,Ndist+1);
irad = zeros(1,Ndist+1);

%Give each SDG mass, initial radius, initial theta, initial velocity.
%Define radial, azimuthal band width.
BandRad = 30;
BandRadWidth = 5;
BandAngle = pi/4;
BandAngleWidth = pi/8;
for DST = 2:Ndist+1
    GM(DST) = Mdist*mcl*Gcon;
    GMtot = GMtot + abs(GM(DST));
    itheta(DST) = rand*2*pi; %Azimuth
    iphi(DST) = BandAngle + rand*BandAngleWidth; %Inclination
    irad(DST) = BandRad + rand*BandRadWidth;
    ivtan(DST) = sqrt(GM(1)/(irad(DST))^(2*-0.3)); %v1 - circular orbital velocity at radius irad
end

%Stellar Disk Parameters
H = 0.1;
N1 = 5;
del = -0.3;
NREM = 1; %Inner Edge
NRNG = 70;
NMID =  NRNG; %Radius of const particle rings

%Distribute Particles in Radius and Longitude
II = 0;
    for J = NREM+1:NRNG
        NPER = N1*(J/NRNG)*(RMAX/H);        % CONSTANT DEN.
        for I = 1:NPER
            II = II+1;
            R(II)=RMAX*J/NRNG;
            PHI(II)=2.0*pi*I/NPER;
            VTHETA(II)= sqrt(GM(1)/R(II)^(2*del));
        end
	end
ntot = II;

%Distribute Randomly in z-velocity
vz0 = 1.5;
z20 = 0.5;
    for k=1:ntot
        rsign = round(10*rand);
        vz2(k) = vz0*rand*(-1)^rsign;
        rsign2 = round(10*rand);
        z2(k) = z20*rand*(-1)^rsign2;
    end
    
%Polar to Cartesian
%Particles
    q = R;
    x2 = q.*cos(PHI);
    y2 = q.*sin(PHI);
    x0 = x2;
    y0 = y2;
    vx2 = -VTHETA.*sin(PHI);
    vy2 = VTHETA.*cos(PHI);
%SDGs
    phioftheta = zeros(1,Ndist+1);
    %loop to find starting inclination of SDG with initial theta
    for nphi = 2:Ndist+1
        if itheta(nphi) <= pi
            phioftheta(nphi) = iphi(nphi) + (1-(2/pi)*iphi(nphi))*itheta(nphi);
        else
            phioftheta(nphi) = pi - iphi(nphi) - (1-(2/pi)*iphi(nphi))*(itheta(nphi)-pi);
        end
    end
    xdist2 = irad.*cos(itheta).*sin(iphi);
    ydist2 = irad.*sin(itheta).*sin(iphi);
    zdist2 = irad.*cos(phioftheta);
    phioftheta = zeros(1,Ndist+1);
%     vyd2 = ivtan;
%     vxd2 = zeros(1, Ndist+1);
%     vzd2 = zeros(1,Ndist+1);
    vxd2 = -1.*ivtan.*sin(itheta).*sin(phioftheta); %v terms assume circular orbits
    vyd2 = ivtan.*cos(itheta); 
    vzd2 = ivtan.*sin(itheta).*cos(phioftheta); 
    

%Graphics Positioning Stuff
hgap = 0.35;
vgap = 0.32;
side = 0.27;
%
lf1 = 0.2;
lf2 = 0.2 + hgap;
left = ([lf1, lf2, lf1, lf2, lf1, lf2]);
%
bt1 = 0.7;
bt2 = 0.7 - vgap;
bt3 = bt2 - vgap;
bottom = ([bt1, bt1, bt2, bt2, bt3, bt3]);

%Compute Particle Positions at Later Times
tstep = 5; %Timestep
tend = 100; %End Time
tstart = 0.0001; %Start Time
tout = [tstart, tstep:tstep:tend];

nsub = size(tout);

%Numerical Orbits
saveiter = 0; %Counter variable for save code
t1 = 0;
for isub=1:nsub(2)
    t2 = tout(isub);
    %SDG Orbits
    for d = 2:(Ndist+1)
        dd = d;
    %Initial Values
        xdist1 = xdist2(d);
        ydist1 = ydist2(d);
        zdist1 = zdist2(d);
        vxd1 = vxd2(d);
        vyd1 = vyd2(d);
        vzd1 = vzd2(d);
    %Integrator
        [t,sdgn] = integrosdg(GM,xdist1,ydist1,zdist1,vxd1,vyd1,vzd1,Ndist,t1,t2,del);
        lsdg = sdgn(:,1);
        lastd = length(lsdg); %Find height of sdgn
        xdist2(dd) = sdgn(lastd,1);
        ydist2(dd) = sdgn(lastd,2);
        zdist2(dd) = sdgn(lastd,3);
        vxd2(dd) = sdgn(lastd,4);
        vyd2(dd) = sdgn(lastd,5);
        vzd2(dd) = sdgn(lastd,6);
        clear sdgn
    end
    %Advance particles over timestep.
    jj = 0;
    for j=1:1:ntot
        jj=jj+1;
    %Initial Values
        x1 = x2(j);
        y1 = y2(j);
        vx1 = vx2(j);
        vy1 = vy2(j);
        z1 = z2(j);
        vz1 = vz2(j);
    %Integrator
        [t,uvn] = integro3dstars(del, Ndist, GM, t1, t2, x1, y1, vx1, vy1, z1, vz1, xdist2, ydist2, zdist2);
        un = uvn(:,1);
        lastd=length(un);
        x2(jj) = uvn(lastd,1);
        y2(jj) = uvn(lastd,2);
        vx2(jj) = uvn(lastd,3);
        vy2(jj) = uvn(lastd,4);
        z2(jj) = uvn(lastd,5);
        vz2(jj) = uvn(lastd,6);
    %Star Formation/Death
          if j==2000
            t3 = t;
            x3 = un;
            y3 = uvn(:,2);
            elseif j==6000
            t4=t;
            x4 = un;
            y4 = uvn(:,2);
            elseif j==12000
            t5=t;
            x5 = un;
            y5 = uvn(:,2);
            elseif j==15000
            t6=t;
            x6 = un;
            y6 = uvn(:,2);
          end %If
    end %j for
clear uvn

%Dwarf Galaxy Position Vector Components (for orbital path tracing)
% xvec(isub) = xdist2(Ndist+1);
% yvec(isub) = ydist2(Ndist+1);
% zvec(isub) = zdist2(Ndist+1);

figure(isub)
subplot(2,2,1) %XY (Top-down)
    plot(x2,y2,'k+','MarkerSize', 1)
    hold on
    for idx = 2:(Ndist+1)
        plot(xdist2(idx), ydist2(idx), 'go', 'MarkerSize', 10)
        hold on
%         plot(vxd2(idx),vyd2(idx),'go','MarkerSize',10)
%         hold on
    end %XY Disturber plot loop
    hold on
    grid on
    axis ([-12 12 -12 12])
    axis ('square')
    xlabel('X')
    ylabel('Y')
    set(gca,'FontSize',[16])
    set(gca, 'FontName', 'Helvetica', 'FontWeight', 'Bold')
    set(gca, 'TickLength', [0.05 0.025])
    ttt = num2str(t(lastd));
    text(-9, 9, ttt,'Fontsize',[18])
    masslabel = num2str(Mdist);
    text(8, 9, masslabel,'Fontsize',[18])
    hold off
subplot(2,2,2) %XZ (Side 1)
    plot(x2,z2,'k+','MarkerSize', 1)
    hold on
    for idx = 2:(Ndist+1)
        plot(xdist2(idx), zdist2(idx), 'go', 'MarkerSize', 10)
        hold on
%         plot(vxd2(idx),vyd2(idx),'go','MarkerSize',10)
%         hold on
    end %XZ Disturber plot loop
    grid on
    axis ([-12 12 -12 12])
    axis ('square')
    xlabel('X')
    ylabel('Z')
    set(gca,'FontSize',[16])
    set(gca, 'FontName', 'Helvetica', 'FontWeight', 'Bold')
    set(gca, 'TickLength', [0.05 0.025])
    ttt = num2str(t(lastd));
    text(-9, 9, ttt,'Fontsize',[18])
    hold off
subplot(2,2,3) %YZ (Side 2)
    plot(y2,z2,'k+','MarkerSize', 1)
    hold on
    for idx = 2:(Ndist+1)
        plot(ydist2(idx), zdist2(idx), 'go', 'MarkerSize', 10)
        hold on
        %         plot(vxd2(idx),vyd2(idx),'go','MarkerSize',10)
        %         hold on
    end %YZ Disturber plot loop
    grid on
    axis ([-12 12 -12 12])
    axis ('square')
    xlabel('Y')
    ylabel('Z')
    set(gca,'FontSize',[16])
    set(gca, 'FontName', 'Helvetica', 'FontWeight', 'Bold')
    set(gca, 'TickLength', [0.05 0.025])
    ttt = num2str(t(lastd));
    text(-9, 9, ttt,'Fontsize',[18])
    hold off
subplot(2,2,4)
    %stairs(Rr,nR,'r-','MarkerSize',3)
    hold off
    axis ([0 15 0 3.0])
    axis ('square')
    xlabel('R')
    ylabel('log10(Surface Density)')
    set(gca,'FontSize',[16])
    set(gca, 'FontName', 'Helvetica', 'FontWeight', 'Bold')
%Maximize and save figures
ih = figure(isub);
set(ih,'WindowState','maximized')
savelocale = 'C:\Users\Alexander\Documents\ISU Senior Year\Research\Recent Research Runs';
foldername = strcat('SDGGalaxyTest',num2str(numrun));
subsavelocale = fullfile(savelocale,foldername);
subfoldername = strcat(foldername,'figures');
if saveiter==0
    locname = fullfile(savelocale,foldername);
    sublocname = fullfile(subsavelocale,subfoldername);
    mkdir(locname)
    mkdir(sublocname)
    saveiter = 1;
end
timestr = num2str(t(lastd));
TF = contains(timestr,'.');
if TF == 1
    timestr = strrep(timestr,'.','p');
end
jpegsavename = strcat(timestr,'.jpg');
jpegfullsave = fullfile(savelocale,foldername,jpegsavename);
figsavename = strcat(timestr,'.fig');
figfullsave = fullfile(savelocale,foldername,subfoldername,figsavename);
saveas(ih, jpegfullsave)
saveas(ih, figfullsave)

%3D Plot
figure(isub+100)
    scatter3(x2,y2,z2,1,'k+')
    hold on
    for idx = 2:(Ndist+1) %3D Plot loop
        scatter3(xdist2(idx), ydist2(idx),zdist2(idx),'b')
        hold on
    end
    axis([-40 40 -40 40 -40 40])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    tlab = strcat('Disturber Mass =',{' '},num2str(Mdist),{'   '},'Time =',{' '},num2str(t(lastd)));
    title(tlab)
    set(gca,'FontSize',[16])
    set(gca, 'FontName', 'Helvetica', 'FontWeight', 'Bold')
    set(gca, 'TickLength', [0.05 0.025])
    grid on
    hold off
isoplot = figure(isub+100);
set(isoplot,'WindowState','maximized')
opfigsavename = strcat('3dplot',timestr,'.fig');
opjpegsavename = strcat('3dplot',timestr,'.jpg');
opfigfullsave = fullfile(savelocale,foldername,subfoldername,opfigsavename);
opjpegfullsave = fullfile(savelocale,foldername,opjpegsavename);
saveas(isoplot,opfigfullsave)
saveas(isoplot,opjpegfullsave)

t1 = t2;
end %Loop over particles.

% figure(100)
%     plot3(xvec,yvec,zvec)
%     axis([-15 15 -15 15 -15 15])
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     grid on
% orbitplot = figure(100);
% set(orbitplot,'WindowState','maximized')
% opfigsavename = strcat('orbitplot',num2str(numrun),'.fig');
% opjpegsavename = strcat('orbitplot',num2str(numrun),'.jpg');
% opfigfullsave = fullfile(savelocale,foldername,subfoldername,opfigsavename);
% opjpegfullsave = fullfile(savelocale,foldername,opjpegsavename);
% saveas(orbitplot,opfigfullsave)
% saveas(orbitplot,opjpegfullsave)
    






