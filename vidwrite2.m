function vidwrite2(Xfilt,handles,aeR,h)
f = figure(14);
file = handles.savefigname.String;
path = handles.savefolder.String;
v = VideoWriter([path '\' file]);
% v = VideoWriter(file);
v.FrameRate = str2double(handles.framerate.String);
open(v);

n = size(Xfilt,3);
p = 'n';
imagesc(Xfilt(:,:,1))
colormap(h);
if ~isnan(aeR)
    caxis(aeR)
else
    caxis([min(Xfilt(:)) max(Xfilt(:))])
end
for i = 1:n

    title([p ' = ' num2str(i)])
    set(f.CurrentAxes.Children,'CData',Xfilt(:,:,i))
%     if ~isempty(handles.movietext.String)
%         text(1,15,handles.movietext.String,'Color','white','FontSize',14)
%     end
    
    frame = getframe;
    writeVideo(v,frame);
end

close(v)
