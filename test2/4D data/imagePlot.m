function imagePlot(z, varargin)

opts = process_options(varargin, ...
                       'dims', size(z), ...
                       'cmap','default', ...
                       'coloraxis', [], ...
                       'titleStr', '', ...
                       'refVec', [], ...
                       'cbar', true, ...
                       'centercaxis',true,...
                       'newFigure', true, ...
                       'specPlot', false, ...
                       'zoomFactor', 1, ...
                       'xLabel', [], ...
                       'yLabel', [], ...
                       'xAxis', [], ...
                       'yAxis', []);

%Unpack option variables into workspace
v2struct(opts);

if(length(dims) > 1 && dims(1)==1 && dims(2) > 1)
    dims(1) = dims(2); dims(2) = 1;
    z = z(:);
end

if(newFigure)
    figure;
end

if(~isempty(refVec) && ~isreal(refVec))
    disp('imagePlot: reference vector is complex, only taking real part');
    refVec = real(refVec);
end

if(~isreal(z))
    disp('imagePlot: data is complex, only plotting real part');
    z = real(z);    
end


if(~isempty(refVec) && ~specPlot)
   relerr = norm(z(:) - refVec(:))/norm(refVec(:));
   titleStr = [titleStr, ' relerr = ' num2str(relerr)];
end

z = reshape(z,dims);
if (specPlot)
    if(size(dims,2) > 1)
        z = fftshift(abs(fft2(z)))/numel(z);
    else
        z = fftshift(abs(fft(z)))/numel(z);
    end
end
if(~isempty(xAxis) && ~isempty(yAxis))
    imagesc(xAxis, yAxis, z);
else
    imagesc(z);
end
zoom reset;
zoom(zoomFactor);
colormap(cmap);
if(~isempty(coloraxis))
    caxis(coloraxis);
else
    coloraxis = caxis;
end
if(cbar)
    colorbar;
end
if(centercaxis)
    m = min(abs(coloraxis));
    caxis([-m m]);
end
title(titleStr); xlabel(xLabel); ylabel(yLabel);

end