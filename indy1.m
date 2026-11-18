addpath('/Users/lyla/Downloads/JP1/')
addpath('/Users/lyla/Downloads/JP1/slepian_alpha-master/')

clear


% Get data
    % Earth Data
    load topo

    % CMB Data
    filename='wmap_mollweide_band_imap_r9_9yr_K_v5.fits';
    cmb=fitsread(filename, 'image');


% Do Calculations
    % Whole Earth Calculations
    lmcosi=xyz2plm(topo,120);
    [Sl,l]=plm2spec(lmcosi);
    % Partial Earth Calculations
    topo2=topo;
    oshan=topo<0;
    oshanf=sum(oshan(:))/prod(size(oshan));
    topo2(oshan)=0;
    lmcosi2=xyz2plm(topo2,120);
    [Sl2,l2]=plm2spec(lmcosi2);

    % Whole CMB Calculations
    lmcosi3=xyz2plm(cmb,120);
    [Sl3,l3]=plm2spec(lmcosi3);
    % Partial CMB Calculations
    cmb2=cmb;
    belt=cmb>0.1;
    beltf=sum(belt(:))/prod(size(belt));
    cmb2(belt)=0
    lmcosi4=xyz2plm(topo2,120);
    [Sl4,l4]=plm2spec(lmcosi4);


% Graphics
    % Whole Earth
    figure(1)
    clf
    subplot(221)
    imagesc(topo)
    colorbar;
    xticks([50 150 250 350])
    xticklabels({5725, 171450, 28571, 40000})
    xlabel('Distance (km)')
    yticks([45 90 135 180])
    yticklabels({3000, 6000, 9000, 12000})
    ylabel('Distance (km)')
    title('Whole Earth Topography');
    axis xy image
    subplot(222)
    loglog(l,Sl)
    title('Whole Earth Topography Spectral Density');
    % Partial Earth, adjusted for area
    subplot(223)
    imagesc(topo2)
    colorbar;
    clim([0,500])
    xticks([50 150 250 350])
    xticklabels({5725, 171450, 28571, 40000})
    xlabel('Distance (km)')
    yticks([45 90 135 180])
    yticklabels({3000, 6000, 9000, 12000})
    ylabel('Distance (km)')
    title('Partial Earth Topography');
    axis xy image
    subplot(224)
    loglog(l,Sl)
    hold on
    loglog(l2,Sl2/(1-oshanf)^2)
    hold off
    title('Partial Earth Topography Spectral Density');

    % Save the figure using
    print('-dpdf','indy1_1')

    % Whole CMB
    figure(2)
    clf
    subplot(221)
    imagesc(cmb);
    set(gca, 'YDir', 'normal'); 
    clim([-0.2, 0.2]);
    colormap("jet");
    colorbar;
    xticks([1000 2000 3000 4000])
    xticklabels({23, 46, 69, 92})
    xlabel('Distance (Gly)')
    yticks([500 1000 1500 2000])
    yticklabels({11, 22, 33, 44})
    ylabel('Distance (Gly)')
    title('Whole CMB Temp Data');
    subplot(222)
    loglog(l3,Sl3)
    title('Whole CMB Spectral Density')
    % Partial CMB, adjusted for area
    subplot(223)
    imagesc(cmb2)
    colorbar;
    clim([-0.1,0.1])
    xticks([1000 2000 3000 4000])
    xticklabels({23, 46, 69, 92})
    xlabel('Distance (Gly)')
    yticks([500 1000 1500 2000])
    yticklabels({11, 22, 33, 44})
    ylabel('Distance (Gly)')
    title('Partial CMB Temperatures');
    axis xy image
    subplot(224)
    loglog(l3,Sl3)
    hold on
    loglog(l4,Sl4/(1-beltf)^1/100000)
    hold off
    title('Partial CMB Temperature Spectral Density');


    % Save the figure using
    print('-dpdf','indy1_2')