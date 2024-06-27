% function to plot field as surface
function FieldPlot(surface,plot_dis,boundary_box_x,boundary_box_y,f_length,f_width,cohesion_radius,positions,velocities,p_end)

figure;

% number of sampling points
E_x=ceil(diff(boundary_box_x)/plot_dis);
E_y=ceil(diff(boundary_box_y)/plot_dis);

% initial nurbs surface
surface_eval = nrbeval(surface,{linspace(0,1,E_x), linspace(0,1,E_y)});
surface_eval = permute(surface_eval,[2 3 1]);

% consider the particles' field
if isempty(positions)
    % plot initial field
    hold all
    surf(surface_eval(:,:,1),surface_eval(:,:,2),surface_eval(:,:,3),'EdgeColor','none'); axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar
    set(gca,'FontSize',26,'FontName','Times New Roman');
else
    if isempty(velocities)
        error('No velocities input.');
    else
        % calculate particles field
        % distance from sampling points to particles
        dis=@(a,b) sqrt(a.^2 + b.^2);
        surf_x=surface_eval(:,:,1);
        surf_y=surface_eval(:,:,2);
        surf_z=surface_eval(:,:,3);
        for i=1:size(positions,1)
            dis_a=gather(arrayfun(dis, gpuArray(surf_x-positions(i,1)),gpuArray(surf_y-positions(i,2))));
            dis_b=gather(arrayfun(dis, gpuArray(surf_x-p_end(i,1)),gpuArray(surf_y-p_end(i,2))));
            r_ellipse=(dis_a+dis_b-f_length)./2;
            p_field=r_ellipse;
            p_field(r_ellipse<=cohesion_radius)=2.*(1./(1+exp(-2.5*(r_ellipse(r_ellipse<=cohesion_radius)-f_width)))-0.5);
            p_field(r_ellipse>cohesion_radius)=1;
            p_field(r_ellipse<f_width)=0;
            surf_z=surf_z.*p_field;
        end
        hold all
        surf(surface_eval(:,:,1),surface_eval(:,:,2),surf_z,'EdgeColor','none'); axis equal
        xlabel('x');
        ylabel('y');
        zlabel('z');
        colorbar
        set(gca,'FontSize',26,'FontName','Times New Roman');
    end
end

end