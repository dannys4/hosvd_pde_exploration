function printpdf(h,outfilename)
% From StackOverflow post on Cropping PDF figures
% See original post here: 
% https://stackoverflow.com/questions/3801730/get-rid-of-the-white-space-around-matlab-figures-pdf-output
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);
end