clear all; close all;

set(0, 'DefaultTextInterpreter', 'latex')
%script for plotting the level set data

%copy data (REMOVE THIS BEFORE HANDING OFF)
system('scp colton@98.213.50.185:/home/colton/Documents/Northwestern/Research/stokesLS/u\*.out ..');
system('scp colton@98.213.50.185:/home/colton/Documents/Northwestern/Research/stokesLS/v\*.out ..');
system('scp colton@98.213.50.185:/home/colton/Documents/Northwestern/Research/stokesLS/p\*.out ..');

outDir = '../'; %location of data

for k=0:100:1000

    uOut = fopen([outDir 'u_' num2str(k) '.out']);
    vOut = fopen([outDir 'v_' num2str(k) '.out']);
    pOut = fopen([outDir 'p_' num2str(k) '.out']);

    umaxi   = fread(uOut, 1, 'int');
    umaxj   = fread(uOut, 1, 'int');
    u       = fread(uOut, [umaxj, umaxi], 'double');
    ucoords = fread(uOut, umaxi+umaxj, 'double');

    uxp = ucoords(1:umaxi);
    uyp = ucoords(umaxi+1:end);
    [uX, uY] = meshgrid(uxp, uyp);

    vmaxi   = fread(vOut, 1, 'int');
    vmaxj   = fread(vOut, 1, 'int');
    v       = fread(vOut, [vmaxj, vmaxi], 'double');
    vcoords = fread(vOut, vmaxi+vmaxj, 'double');

    vxp = ucoords(1:vmaxi);
    vyp = ucoords(vmaxi+1:end);
    [vX, vY] = meshgrid(vxp, vyp);

    pmaxi   = fread(pOut, 1, 'int');
    pmaxj   = fread(pOut, 1, 'int');
    p       = fread(pOut, [pmaxj, pmaxi], 'double');
    pcoords = fread(pOut, pmaxi+pmaxj, 'double');

    pxp = pcoords(1:pmaxi);
    pyp = pcoords(pmaxi+1:end);
    [pX, pY] = meshgrid(pxp, pyp);

    figure(1)
    subplot(1,3,1)
    contourf(uX, uY, u, 40)

    subplot(1,3,2)
    contourf(vX, vY, v, 40)

    subplot(1,3,3)
    contourf(pX, pY, p, 40)



end

fclose('all');