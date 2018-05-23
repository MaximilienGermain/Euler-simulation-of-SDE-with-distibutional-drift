function video(filename,frames)

myVideo = VideoWriter([filename '.avi']);
myVideo.FrameRate = 1;  % Default 30
myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);

end