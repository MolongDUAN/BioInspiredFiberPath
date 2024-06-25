% Input displacement data Ux, Uy, and interpolant function UxF(x,y), UyF(x,y)
% or NURBS surface data p_xyz,{px_knot,py_knot}

% Initialization
load('Data/ini_nurbC_1015.mat');
addpath(genpath('nurbs-1.3.7'));
vol_ratio=0.6;  % fiber volume ratio
f_length=0.35;
f_width=0.35;
part_area=150;
% boundary function: fd / xy range, for rectangle case 0<=x<=10 0<=y<=15
num_agents = ceil(part_area*vol_ratio/(f_length*f_width));   % particle #
max_speed = 0.35;
separation_radius = 0.6;   % searching radius 0.6
cohesion_radius = 1.5; % 1.5
separation_weight = 0.4; % 0.4
alignment_weight = 0.8; % 0.8
cohesion_weight = 0.3; % 0.3
fieldGrad_weight = 0.7; % 0.8, field 0.7
p_field_weight = 1;

min_speed=0.1;

% Initial position, velocity, and energy field
positions = [rand(num_agents,1)*10,rand(num_agents,1)*15];  % random position
velocities = [zeros(num_agents,1),0.1*ones(num_agents,1)];   % initial speed y=0.1
surface = nrbmak(p_xyz,{px_knot,py_knot});
dsurface = nrbderiv (surface);

% plot boundary
plot([0 10 10 0 0],[0 0 15 15 0],'-b','LineWidth',1); axis equal; axis off; hold on

% main
p_result=cell(300,2);
tic;
for t = 1:300
    % store positions and velocities
    p_result{t,1}=positions;
    p_result{t,2}=velocities;
    
    % particle directions and particle ending points
    p_angle=atan2(velocities(:,2),velocities(:,1));
    p_end=[positions(:,1)-f_length*cos(p_angle),positions(:,2)-f_length*sin(p_angle)];
    
    % neighbors of each particle
    neighbors_s = cell(num_agents, 1);
    neighbors_c = cell(num_agents, 1);
    neighbors_e = cell(num_agents, 1);
    % calculate distance
    n = size(positions, 1);
    vec_x = zeros(n, n);
    vec_y = zeros(n, n);
    end_x = zeros(n, n);
    end_y = zeros(n, n);
    for i = 1:n-1
        for j = i+1:n
            vec_x(i,j) = positions(i,1)-positions(j,1);
            vec_x(j,i)=vec_x(i,j);
            vec_y(i,j) = positions(i,2)-positions(j,2);
            vec_y(j,i)=vec_y(i,j);
            % for ending points
            end_x(i,j) = positions(i,1)-p_end(j,1);
            end_x(j,i)=end_x(i,j);
            end_y(i,j) = positions(i,2)-p_end(j,2);
            end_y(j,i)=end_y(i,j);
        end
    end
    V_x=gpuArray(vec_x);
    V_y=gpuArray(vec_y);
    dis=@(a,b) sqrt(a.^2 + b.^2);
    D=arrayfun(dis, V_x,V_y);
    d=gather(D);
    Ve_x=gpuArray(end_x);
    Ve_y=gpuArray(end_y);
    De=arrayfun(dis, Ve_x,Ve_y);
    de=gather(De);
    
    for i = 1:num_agents
        for j = 1:num_agents
            if i ~= j && d(i,j) <= separation_radius
                neighbors_s{i} = [neighbors_s{i}; j];
            end
            if i ~= j && d(i,j) <= cohesion_radius
                neighbors_c{i} = [neighbors_c{i}; j];
            end
            if i ~= j && de(i,j) <= cohesion_radius
                neighbors_e{i} = [neighbors_e{i}; j];
            end
        end
    end
    
    % velocity from displacement
%     fieldGrad = [UxF(positions),UyF(positions)];
    parfor i=1:num_agents
        pos_x=positions(i,1)/10;
        pos_y=positions(i,2)/15;
        [~,tmp_grad]=nrbdeval (surface, dsurface, {pos_x,pos_y});
        tmp_grad{1}=tmp_grad{1}/norm(tmp_grad{1});
        tmp_grad{2}=tmp_grad{2}/norm(tmp_grad{2});
        fieldGrad(i,:)=[tmp_grad{1}(3),tmp_grad{2}(3)];
    end
    fieldGrad_value=sqrt(sum(fieldGrad.^2, 2));
    fieldGrad_n = fieldGrad./fieldGrad_value;   % normalize
    
    % update velocity and position
    velocities_tmp=[];
    for i = 1:num_agents
        separation = computeSeparation(positions, positions(i,:), neighbors_s{i});
%         alignment = computeAlignment(velocities, neighbors_c{i});
%         cohesion = computeCohesion(positions, positions(i,:), neighbors_c{i});
%         p_field = computeField(f_length, positions, p_end, positions(i,:), neighbors_c{i});
        [p_field,p_grad] = computeField(f_length, positions, velocities, p_end, positions(i,:), neighbors_e{i});
        
        % compare fieldGrad and p_grad
        if abs(fieldGrad_value(i))>=p_grad
            % check fieldGrad direction
            if velocities(i,:)*fieldGrad_n(i,:)'<0 % change direction
                fieldGrad_p=-fieldGrad_n(i,:);
            else
                fieldGrad_p=fieldGrad_n(i,:);
            end
        else % follow the forward particle
            fieldGrad_p=p_field;
        end

%         velocities_tmp(i, :) = -separation_weight*separation + alignment_weight*alignment + cohesion_weight*cohesion + fieldGrad_weight*fieldGrad_n(i,:);
        velocities_tmp(i, :) = -separation_weight*separation + fieldGrad_weight*fieldGrad_p;
    end
    
    velocities = velocities + velocities_tmp;
    
    % limit speed
    norm_vel = sqrt(sum(velocities.^2,2));
    overMax = norm_vel>max_speed;
    velocities(overMax,:) = velocities(overMax,:) ./ norm_vel(overMax,:) * max_speed;
    
    positions = positions + velocities;
    
    % regenerate outside particles
    positions_new=positions;
    parfor i = 1:num_agents
        if positions(i,1)<0 || positions(i,1)>10
            % outside
            if positions(i,2)>15
                if positions(i,1)>=3.5 && positions(i,1)<=6.5
                    positions_new(i,:) = [10*(positions(i,1)-3.5)/3,0.1]; % y=0.1
                else
                    positions_new(i,:) = [rand(1,1)*10,0.1];
                end
            else
                positions_new(i,:) = [rand(1,1)*10,0.1];
            end
            velocities(i,:)=[0 0.5];
        else
            if positions(i,2)<0 || positions(i,2)>15
                positions_new(i,:) = [rand(1,1)*10,0.1];
                velocities(i,:)=[0 0.5];
            end
        end
    end
    positions=positions_new;
    
end
Time_plot=toc;

% visualization
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

% Separation fun
function separation = computeSeparation(positions, cal_p, neighbors)
    separation = zeros(1, 2);
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        separation = separation + (positions(neighbor_index, :) - cal_p) / norm(positions(neighbor_index, :) - cal_p);
    end
    if ~isempty(neighbors)
        separation = separation / length(neighbors);
    end
end

% Alignment fun
function alignment = computeAlignment(velocities, neighbors)
    alignment = zeros(1, 2);
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        alignment = alignment + velocities(neighbor_index, :);
    end
    if ~isempty(neighbors)
        alignment = alignment / length(neighbors);
    end
end

% Cohesion fun
function cohesion = computeCohesion(positions, cal_p, neighbors)
    cohesion = zeros(1, 2);
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        cohesion = cohesion + positions(neighbor_index, :);
    end
    if ~isempty(neighbors)
        cohesion = cohesion / length(neighbors);
        cohesion = (cohesion - cal_p) / norm(cohesion - cal_p);
    end
end

% Field fun
% function field = computeField(f_length, positions, p_end, cal_p, neighbors)
%     % distance to center of ellipse
%     r_ellipse=zeros(length(neighbors),1);
%     field_dir=zeros(length(neighbors),2);
%     for i = 1:length(neighbors)
%         neighbor_index = neighbors(i);
%         dis_a=norm(positions(neighbor_index, :)-cal_p);
%         dis_b=norm(p_end(neighbor_index, :)-cal_p);
%         r_ellipse(i)=(dis_a+dis_b-f_length)/2;
%         field_dir(i,:)=(positions(neighbor_index, :)-cal_p)+(p_end(neighbor_index, :)-cal_p);
%         field_dir(i,:)=field_dir(i,:)/norm(field_dir(i,:));
%     end
%     % consider the fiber width
%     r_ellipse(r_ellipse<0.4)=0;
%     r_ellipse(r_ellipse>=0.4)=r_ellipse(r_ellipse>=0.4)-0.4;
%     
%     field_value=1./(1+exp(-2.5*r_ellipse))-0.5;
% %     field_value=2./(1+exp(-2.5*(r_ellipse-0.175)))-1;
%     field=field_dir.*field_value;
%     field=sum(field,1)/length(neighbors);
% end

function [field,grad] = computeField(f_length, positions, velocities, p_end, cal_p, neighbors)
    % distance to center of ellipse
    r_ellipse=zeros(length(neighbors),1);
    field_dir=zeros(length(neighbors),2);
    
    % only look forward
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        if velocities(neighbor_index,:)*(p_end(neighbor_index, :)-cal_p)'>0 % direction check
            dis_a=norm(positions(neighbor_index, :)-cal_p);
            dis_b=norm(p_end(neighbor_index, :)-cal_p);
            r_ellipse(i)=(dis_a+dis_b-f_length)/2;
            field_dir(i,:)=(positions(neighbor_index, :)-cal_p)+(p_end(neighbor_index, :)-cal_p);
            field_dir(i,:)=field_dir(i,:)/norm(field_dir(i,:));
        end
    end
    % consider the fiber width
    r_ellipse(r_ellipse<0.4)=0;
    r_ellipse(r_ellipse>=0.4)=r_ellipse(r_ellipse>=0.4)-0.4;
    
    % find neighbor particle with max grad
    if ~isempty(neighbors)
        field_dir(r_ellipse==0,:)=[];
        r_ellipse(r_ellipse==0)=[]; % remove close particles
        if ~isempty(r_ellipse)
            field_grad=5*exp(-2.5*r_ellipse) ./ (1+exp(-2.5*r_ellipse)).^2; % gradient
            field_grad=abs(field_grad);
            [grad,grad_index]=max(field_grad);
            field=field_dir(grad_index,:);
        else
            grad=0;
            field=[0,0];
        end
    else
        grad=0;
        field=[0,0];
    end
end
