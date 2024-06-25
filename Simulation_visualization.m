% visualization
load('Output/simulation_1015.mat');

tic;
outputVideo = VideoWriter('Boids_example.mp4', 'MPEG-4');
outputVideo.FrameRate = 30;  % 设置帧率
outputVideo.Quality = 100;  % 设置视频质量（范围：0-100）
open(outputVideo);
for t=1:300
    p_angle=atan2(p_result{t,2}(:,2),p_result{t,2}(:,1));
    particles=[p_result{t,1},p_result{t,1}(:,1)-f_length*cos(p_angle),p_result{t,1}(:,2)-f_length*sin(p_angle)];
    p_x = particles(:, [1 3])';
    p_y = particles(:, [2 4])';
    hp=line(p_x, p_y, 'Color', 'k'); axis equal; axis off;
    axis([0 10 0 15]);
    drawnow;
    M=getframe;
    frame = M.cdata;
    writeVideo(outputVideo, frame);
    delete(hp);
end
close(outputVideo);
Time=toc;