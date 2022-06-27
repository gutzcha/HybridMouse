function vid=mat2wav(filename,audio_file,fr)
if ~exist('fr','var')
    fr=250000;
end
hSpectro=spectro(audio_file);
frame_rate=30;
bin_len=length(audio_file);
sec_len=bin_len/fr;
num_of_frames=frame_rate*sec_len;

new_fr=fr/10;
new_audio_file=audio_file(1:10:end);
    audiowrite([filename,'.wav'],new_audio_file,new_fr);



hSpectro_axes=findobj(hSpectro,'type','axes');
max_z=max(max(hSpectro_axes.Children.ZData));
xlim=hSpectro_axes.XLim;
ylim=hSpectro_axes.YLim;
x_locations=linspace(xlim(1),xlim(2),num_of_frames);
new_axis=axes('Parent',hSpectro,...
              'Position',hSpectro_axes.Position,...
              'YAxisLocation','right',...
              'Color','none',...
              'YColor','none',...
              'XColor','none');
          new_axis.XLim=xlim;
          new_axis.YLim=ylim;
          new_axis.ZLim=hSpectro_axes.ZLim;
    k=line([x_locations(1),x_locations(1)],...
           [0,             ylim(2)],...
           [max_z,         max_z],...
           'Color',        'r',...
           'LineWidth',    3);
    for i = 1:num_of_frames
        k.XData=[x_locations(i),x_locations(i)];
        hSpectro_axes.YLim=[1 100];
        F(i)=getframe(hSpectro_axes);
    end
    new_vid = VideoWriter([filename,'.avi']);
    new_vid.FrameRate = frame_rate;
    
    open(new_vid)
        
        for u=1:num_of_frames
            % add each frame to the video
            frame = F(u).cdata;
            writeVideo(new_vid, frame);
        end
%          close(new_vid)

%     aud_file=audioread([filename,'.wav']);
    vid_aud=vision.VideoFileWriter([filename,'.avi']);
    vid_aud.FrameRate =frame_rate;
    vid_aud.AudioInputPort=true;
%     vid_aud.AudioInputPort=
close hSpectro;
end