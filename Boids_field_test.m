% Flocking simulation for fiber particles based on Boids algorithm
%--------------------------------------------------------------------

%% 0. Initialization
% 0.1 Fiber
vol_ratio=0.6;  % fiber volume ratio [%]
f_length=0.35;  % fiber particle length [mm]
f_width=0.35;   % fiber particle width [mm]

% 0.2 Design domain
load('Data/ini_nurbC_1015.mat');    % input NURBS surface data
addpath(genpath('nurbs-1.3.7'));
boundary_box_x=[0 10];
boundary_box_y=[0 15];
part_area=diff(boundary_box_x)*diff(boundary_box_y);
% boundary function: fd / xy range, for rectangle case 0<=x<=10 0<=y<=15
surface = nrbmak(p_xyz,{px_knot,py_knot});
dsurface = nrbderiv (surface);

% 0.3 Particle
num_agents = ceil(part_area*vol_ratio/(f_length*f_width));   % particle #
max_speed = 0.35;           % max particle movement speed [mm/iteration]
separation_radius = 0.6;    % searching radius for separation [mm], (0.6)
cohesion_radius = 1.5;      % searching radius for cohesion [mm], (1.5)
separation_weight = 0.4;    % 0.4
alignment_weight = 0.8;     % 0.8
cohesion_weight = 0.3;      % 0.3
fieldGrad_weight = 0.7;     % 0.8, field 0.7
% p_field_weight = 1;

% 0.4 Initial position, velocity, iteration
positions = [rand(num_agents,1)*(boundary_box_x(2)-boundary_box_x(1))+boundary_box_x(1),rand(num_agents,1)*(boundary_box_y(2)-boundary_box_y(1))+boundary_box_y(1)];  % random position
velocities = [zeros(num_agents,1),0.1*ones(num_agents,1)];   % initial speed y=0.1
iteration_num = 300;

% 0.5 Options
GPU_flag=1;                 % 0 - off, 1 - on
rule_flag=2;                % 0 - Boids, 1 - separation + field gradient, 2 - field value
plot_flag=0;                % 0 - no plot, 1 - plot field
plot_iter=0;              % 0 - only plot initial field, other - plot field of that iteration
plot_dis=0.02;              % the sampling point distance [mm]

% 0.6 Rotation matrix
R45 = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];
Rn45 = [cos(-pi/4), -sin(-pi/4); sin(-pi/4), cos(-pi/4)];
R90 = [cos(pi/2), -sin(pi/2); sin(pi/2), cos(pi/2)];

%% 1. Main
% 1.0 Plot initial field
if plot_flag
    if rule_flag==1
        FieldPlot(surface,plot_dis,boundary_box_x,boundary_box_y,[],[],[],[],[],[]);
    else
        print('Field is not used in the motion rule.');
    end
end

p_result=cell(iteration_num,2); % Allocate the result store space
tic;
for t = 1:iteration_num
    % 1.1 Store positions and velocities
    p_result{t,1}=positions;
    p_result{t,2}=velocities;
    
    % 1.2 Particle directions and particle ending points
    p_angle=atan2(velocities(:,2),velocities(:,1));
    p_end=[positions(:,1)-f_length*cos(p_angle),positions(:,2)-f_length*sin(p_angle)];
    
    % 1.3 Neighbors of each particle
    % calculate coordinates of vectors between particles 
    vec_x = triu(positions(:,1) - positions(:,1)');
    vec_y = triu(positions(:,2) - positions(:,2)');
    end_x = triu(positions(:,1) - p_end(:,1)'); % from start to end point
    end_y = triu(positions(:,2) - p_end(:,2)');
    % calculate distance between particles
    dis=@(a,b) sqrt(a.^2 + b.^2);               % distance function
    if GPU_flag                                 % use GPU for calculation
        V_x=gpuArray(vec_x);
        V_y=gpuArray(vec_y);
        D=arrayfun(dis, V_x,V_y);
        d=gather(D);                            % distance between particles
        Ve_x=gpuArray(end_x);
        Ve_y=gpuArray(end_y);
        De=arrayfun(dis, Ve_x,Ve_y);
        de=gather(De);                          % distance start to end
    else
        d = zeros(num_agents,num_agents);
        de = zeros(num_agents,num_agents);
        parfor i = 1:num_agents
            for j = i:num_agents
                d(i,j) = dis(vec_x(i,j), vec_y(i,j)); 
                de(i,j) = dis(end_x(i,j), end_y(i,j));
            end
        end
    end
    % get symmetric matrix
    d=triu(d) + triu(d, 1)';
    de=triu(de) + triu(de, 1)';
    % find the neighboring particles
    neighbors_s = cell(num_agents, 1);
    neighbors_c = cell(num_agents, 1);
    neighbors_e = cell(num_agents, 1);
    for i = 1:num_agents
        for j = i+1:num_agents
            if d(i,j) <= separation_radius
                neighbors_s{i} = [neighbors_s{i}; j];
                neighbors_s{j} = [neighbors_s{j}; i]; % symmetric position
            end
            if d(i,j) <= cohesion_radius
                neighbors_c{i} = [neighbors_c{i}; j];
                neighbors_c{j} = [neighbors_s{j}; i];
            end
            if de(i,j) <= cohesion_radius
                neighbors_e{i} = [neighbors_e{i}; j];
                neighbors_e{j} = [neighbors_s{j}; i];
            end
        end
    end
    
    % 1.4 Particles motion
    % update velocity and position
    velocities_tmp=[];
    parfor i = 1:num_agents
        separation = computeSeparation(positions, positions(i,:), neighbors_s{i});
        if rule_flag==0 % use typical Boids algorithm
            alignment = computeAlignment(velocities, neighbors_c{i});
            cohesion = computeCohesion(positions, positions(i,:), neighbors_c{i});
            % obtain the change of the velocity
            velocities_tmp(i, :) = -separation_weight*separation + alignment_weight*alignment + cohesion_weight*cohesion;
            
        elseif rule_flag==1
            % field value and gradient from NURBS surface
            pos_x=positions(i,1)/10;
            pos_y=positions(i,2)/15;
            [nrb_pnt,nrb_grad]=nrbdeval (surface, dsurface, {pos_x,pos_y});
            nrb_grad{1}=nrb_grad{1}/norm(nrb_grad{1});  % normalize
            nrb_grad{2}=nrb_grad{2}/norm(nrb_grad{2});
            % particle's field value and gradient
            [p_field,p_grad] = computeField(f_length, f_width, positions, velocities, p_end, positions(i,:), neighbors_e{i});
            % calcualte field gradient
            fieldGrad=[nrb_grad{1}(3)*p_field+nrb_pnt(3)*p_grad(1),nrb_grad{2}(3)*p_field+nrb_pnt(3)*p_grad(2)];
            fieldGrad_n=fieldGrad/norm(fieldGrad);
            % check fieldGrad direction
            if velocities(i,:)*fieldGrad_n'<0 % change direction
                fieldGrad_p=-fieldGrad_n;
            else
                fieldGrad_p=fieldGrad_n;
            end
            % obtain the change of the velocity
            velocities_tmp(i, :) = -separation_weight*separation + fieldGrad_weight*fieldGrad_p;
            
        elseif rule_flag==2
            fwp = calculate_points(R45, Rn45, R90, velocities(i, :), positions(i, :), max_speed);
            pos_x=fwp(:,1)'/10;
            pos_y=fwp(:,2)'/15;
            pos_x(pos_x>1)=1;
            pos_x(pos_x<0)=0;
            pos_y(pos_y>1)=1;
            pos_y(pos_y<0)=0;
            field_value=zeros(1,5);
            for j=1:5
                [nrb_pnt,~]=nrbeval (surface, {pos_x(j),pos_y(j)});
                field = computeFieldV(f_length, f_width, positions, velocities, p_end, positions(i, :), fwp(j,:), neighbors_e{i});
                field_value(j)=nrb_pnt(3)*field;
            end
            if all(field_value==0)
                fieldGrad_p=[0 0];
            else
                [~,fwp_index]=max(field_value);
                fieldGrad_p=(fwp(fwp_index,:)-positions(i, :))/norm((fwp(fwp_index,:)-positions(i, :)));
            end
            % obtain the change of the velocity
            velocities_tmp(i, :) = -separation_weight*separation + fieldGrad_weight*fieldGrad_p;
        end
    end
    % update velocities
    velocities = velocities + velocities_tmp;
    % limit speed
    norm_vel = sqrt(sum(velocities.^2,2)); % calculate velocity value
    overMax = norm_vel>max_speed;
    velocities(overMax,:) = velocities(overMax,:) ./ norm_vel(overMax,:) * max_speed;
    % update positions
    positions = positions + velocities;
    
    % regenerate outside particles (non-universal method, only for rectangle case)
    positions_new=positions;
    parfor i = 1:num_agents
        if positions(i,1)<boundary_box_x(1) || positions(i,1)>boundary_box_x(2) 
            % outside x boundary
            if positions(i,2)>boundary_box_y(2)
                if positions(i,1)>=3.5 && positions(i,1)<=6.5
                    positions_new(i,:) = [10*(positions(i,1)-3.5)/3,0.1]; % y=0.1
                else
                    positions_new(i,:) = [rand(1,1)*10,0.1];
                end
            else
                positions_new(i,:) = [rand(1,1)*10,0.1];
            end
            velocities(i,:)=[0 max_speed];
        else
            if positions(i,2)<boundary_box_y(1) || positions(i,2)>boundary_box_y(2)
                positions_new(i,:) = [rand(1,1)*10,0.1];
                velocities(i,:)=[0 max_speed];
            end
        end
    end
    positions=positions_new;
    
    % 1.5 Plot field of an iteration
    if plot_flag
        if plot_iter==t
            FieldPlot(surface,plot_dis,boundary_box_x,boundary_box_y,f_length,f_width,cohesion_radius,positions,velocities,p_end);
            input('Check the plot, press Enter to continue');
        end
    end
end
Time=toc;

%% 2. Save simulation result
currentDir = pwd;
outputDir = fullfile(currentDir, 'Output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(fullfile(outputDir, 'simulation_1015.mat')); % set name of data

%% functions
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
        % calculate average neighboring position
        cohesion = cohesion / length(neighbors);
        cohesion = (cohesion - cal_p) / norm(cohesion - cal_p);
    end
end

% Field fun
function [field,grad] = computeField(f_length, f_width, positions, velocities, p_end, cal_p, neighbors)
    % distance to center of ellipse
    r_ellipse=zeros(length(neighbors),1);   % distance
    field_dir=zeros(length(neighbors),2);   % direction
    
    % only look forward
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        if velocities(neighbor_index,:)*(p_end(neighbor_index, :)-cal_p)'>0 % direction check
            dis_a=norm(positions(neighbor_index, :)-cal_p); % distance to starting point
            dis_b=norm(p_end(neighbor_index, :)-cal_p);     % distance to ending point
            r_ellipse(i)=(dis_a+dis_b-f_length)/2;
            field_dir(i,:)=(positions(neighbor_index, :)-cal_p)+(p_end(neighbor_index, :)-cal_p);
            field_dir(i,:)=field_dir(i,:)/norm(field_dir(i,:));
        end
    end
    % consider the fiber width (remove a circle area --> not ideal)
    r_ellipse(r_ellipse<f_width)=0;
    r_ellipse(r_ellipse>=f_width)=r_ellipse(r_ellipse>=f_width)-f_width;
    
    % find neighbor particle with max grad
    if ~isempty(neighbors)
        field_dir(r_ellipse==0,:)=[];
        r_ellipse(r_ellipse==0)=[]; % remove close particles
        if ~isempty(r_ellipse)
            field_grad=5*exp(-2.5*r_ellipse) ./ (1+exp(-2.5*r_ellipse)).^2; % gradient
            field_grad=abs(field_grad);
            [grad,grad_index]=max(field_grad);
            grad=grad*field_dir(grad_index,:);  % components in x and y directions
            field=2*(1./(1+exp(-2.5*r_ellipse(grad_index)))-0.5);
        else
            grad=[0 0];
            field=0.01; % assume the field value will reduce to 1% on the particle
        end
    else
        grad=[0 0];
        field=1;
    end
end

% Forward search points fun
function fwp = calculate_points(R45, Rn45, R90, velocities, positions, max_speed)

v_norm=velocities/norm(velocities);
fwp(1,:) = positions + max_speed .* v_norm;
fwp(2,:) = positions + max_speed .* (R45 * v_norm')';
fwp(3,:) = positions + max_speed .* (Rn45 * v_norm')';
fwp(4,:) = positions + max_speed .* (R90 * v_norm')';
fwp(5,:) = positions - max_speed .* (R90 * v_norm')';

end

% Field value fun
function field = computeFieldV(f_length, f_width, positions, velocities, p_end, start_p, cal_p, neighbors)
    % distance to center of ellipse
    r_ellipse=zeros(length(neighbors),1);   % distance
    
    % only look forward
    for i = 1:length(neighbors)
        neighbor_index = neighbors(i);
        if velocities(neighbor_index,:)*(p_end(neighbor_index, :)-start_p)'>0 % direction check
            dis_a=norm(positions(neighbor_index, :)-cal_p); % distance to starting point
            dis_b=norm(p_end(neighbor_index, :)-cal_p);     % distance to ending point
            r_ellipse(i)=(dis_a+dis_b-f_length)/2;
        end
    end
    % consider the fiber width (remove a circle area --> not ideal)
    r_ellipse(r_ellipse<=f_width)=0;
    r_ellipse(r_ellipse>f_width)=r_ellipse(r_ellipse>f_width)-f_width;
    
    % find neighbor particle with max grad
    if ~isempty(neighbors)
        r_ellipse(r_ellipse==0)=[]; % remove close particles
        if ~isempty(r_ellipse)
            field=prod(2*(1./(1+exp(-2.5*r_ellipse))-0.5));
        else
            field=0.01; % assume the field value will reduce to 1% on the particle
        end
    else
        field=1;
    end
end