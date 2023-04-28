function display_brain_img(img,background_img,z_axis,x_axis,title_name,pv)
% Displays a freground image (img) on top of a background image
% You may adjust the code if you find it necessary

figure; ax1 = axes;
imagesc(ax1,x_axis,z_axis,background_img); colormap(ax1,'gray'); hold on;
ax2 = axes; hb = imagesc(ax2,x_axis,z_axis,img);

switch pv

    case 1

    hb.AlphaData = 5*img; cmap = 'hot';
    ax2.Visible = 'off'; 
    c = colorbar;  c.FontSize = 12;
    c.Ticks = linspace(min(nonzeros(img(:))),max(img(:)),4);
    c.TickLabels = round(linspace(min(nonzeros(img(:))),max(img(:)),4),2);
    set(get(c,'Title'),'String','PCC');
    caxis([min(nonzeros(img(:))) max(img(:))]);

    case 2

    cmap = lines(1); hb.AlphaData = .4*img;

end

colormap(ax2,cmap);
ax2.Visible = 'off';
ax2.YDir = 'reverse';
linkprop([ax1,ax2],'Position');
xlabel(ax1,'Width [mm]'); ylabel(ax1,'Depth [mm]');
title(ax1,title_name);
ax1.YTick = 2:7; 
set(ax1,'FontSize',14);

end