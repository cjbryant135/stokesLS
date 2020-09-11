clear all; close all;
set(0, 'DefaultTextInterpreter', 'latex')
%script for plotting the level set data

%copy data (REMOVE THIS BEFORE HANDING OFF)
system('scp colton@98.213.50.185:/home/colton/Documents/Northwestern/Research/stokesLS/phi\*.out ..');

outDir = '../'; %location of data


for k=0:100:1000
    %open file
    phiOut = fopen([outDir 'phi_' num2str(k) '.out']); 

    %read file
    tmaxi   = fread(phiOut, 1, 'int');
    tmaxj   = fread(phiOut, 1, 'int');
    klen    = fread(phiOut, 1, 'int');
    phi     = fread(phiOut, [tmaxj, tmaxi], 'double');
    coords  = fread(phiOut, tmaxi+tmaxj, 'double');

    %plot it
    xPhys = coords(1:tmaxi);
    yPhys = coords(tmaxi+1:tmaxi+tmaxj);
    [xGrid, yGrid] = meshgrid(xPhys, yPhys);
    figure(1) 
    cla()
    contourf(xGrid, yGrid, phi, -1:0.1:1); %contour plot of phi
    hold on
    contour(xGrid(3:end-2,3:end-2), yGrid(3:end-2,3:end-2), phi(3:end-2,3:end-2), [0,0], '-r', 'LineWidth', 2) %overlay the interface
    set(gca, 'FontSize', 22)
    %title(['$t = $' num2str(k*0.001)]) %current time
    colorbar
    caxis([-1 1])
    
    
end
fclose('all');
