function interpolation4dFigures(data,interp_data, known_srcs, ...
                                varargin)
% Displays common source gathers for the interpolated data and
% compares to the known data if it exists
%
% Arguments:
%    data          - known 4d data (with zeros) in either 
%                       [src x src y rec x rec y] format  OR
%                       [src x rec x src y rec y] format
%    interp_data   - interpolated data
%    known_srcs    - (2 x num_srcs) array of 2d coordinates where the known sources are (x coordinate in the first row, y coordinate in the second row)
%    saveDir       - directory to save figures in (default: [], do not save figures)
%    isPermuted    - true if the data has dimensions [src x rec x src y rec y] (default)
%                  - false if the data has dimensions [src x src y rec x rec y]
%    dispTitles    - true if figures should have titles (default: false)
%    axisLabels    - true if axes should have labels (default: true)
%    samplePoints  - (2 x srcs) array of 2d source coordinates where the figures should be drawn (x coordinate in the first row, y coordinate in the second row)
%    upsampleRecs  - receiver grid size to display the figures in (using fourier interpolation) (default: 0, no interpolation)

   opts = process_options(varargin, ...
                          'dispTitles',false,...
                          'axisLabels',true,...
                          'saveDir',[], ...
                          'samplePoints',[],...
                          'isPermuted',true, ...
                          'upsampleRecs',0);
   v2struct(opts);
   dims = size(data);
   nsrcs = dims(1);
   if(isempty(samplePoints))
       error('Need sample points');
   end

   xsol_full = reshape(interp_data,dims);
   if(isPermuted)
       p = [1 3 2 4];
       xsol_full = permute(xsol_full,p);
       data = permute(data,p);
   end
      
   saveFigs = ~isempty(saveDir);
   if(saveFigs && saveDir(end) ~= '/')
       saveDir = [saveDir '/'];
   end
   
   if axisLabels
       xlabel = 'Receiver y'; ylabel = 'Receiver x';
   else
       xlabel = ''; ylabel = '';
   end
   
   % Plotting functions
   image_src = @(A) imagePlot(A,'cbar',true,'xLabel',xlabel, 'yLabel',ylabel,'centercaxis',true);
   image_interp = @(A,cxis) imagePlot(A,'cbar',true,'xLabel',xlabel, 'yLabel',ylabel,'centercaxis',true,'coloraxis',cxis);

   
   for i=1:length(samplePoints)
       srcx = samplePoints(1,i); 
       srcy = samplePoints(2,i); 
       
       saveFigure = @(filename,postfilename) print(gcf,'-depsc',[saveDir filename '-[' num2str(srcx) '-' num2str(srcy) ']-' postfilename '.eps']);
       
       isknown = find(srcx == known_srcs(1,:) & srcy == known_srcs(2,:));
       slice = squeeze(xsol_full(srcx,srcy,:,:));       
       if(length(isknown) > 0)
           %Known shot record
           knownSlice = squeeze(data(srcx,srcy,:,:));
           snr_slice = SNR(vec(knownSlice),vec(slice));
           if upsampleRecs > 0
               knownSlice = interpft(interpft(knownSlice,upsampleRecs,1),upsampleRecs,2);
           end
           image_src(knownSlice);
           trueaxis = caxis;
           if dispTitles
               title(['True data - (srcx,srcy) = (' num2str(srcx) ',' num2str(srcy) ')']);
           end
           if saveFigs
               saveFigure('truedata','');
           end
           
           %Interpolation result
           if upsampleRecs > 0
               slice = interpft(interpft(slice,upsampleRecs,1),upsampleRecs,2);
           end
           image_interp(slice,trueaxis);
           
           if dispTitles
               title(['Interpolated data - SNR' num2str(snr_slice,3) 'dB']);
           end
           if saveFigs
               saveFigure('interpdata',['snr-' num2str(snr_slice,3) ] ) ;
           end
           
           %Difference
           d = slice - knownSlice;
           if(upsampleRecs > 0)
               d = interpft(interpft(d,upsampleRecs,1),upsampleRecs,2);
           end
           
           image_interp(d,trueaxis);
           if dispTitles
               title('Difference');
           end
           if saveFigs
               saveFigure('difference','');
           end
           
       else
           if upsampleRecs > 0
               slice = interpft(interpft(slice,upsampleRecs,1),upsampleRecs,2);
           end
           image_src(slice);
           if dispTitles
               title(['Interpolated data - (srcx,srcy) = (' num2str(srcx) ',' num2str(srcy) ')']);
           end
           if saveFigs
               saveFigure('interpdata','');
           end
       end
   end
   
end