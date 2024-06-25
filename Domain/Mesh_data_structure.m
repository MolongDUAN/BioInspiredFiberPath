function [e_p,e_t,p_e,p_t,t_e]=Mesh_data_structure(p,t)

% edge to node
e_p=zeros(3*size(t,1),2);
for i=1:size(t,1)
    e_p(3*i-2,:)=[t(i,1) t(i,2)];
    e_p(3*i-1,:)=[t(i,1) t(i,3)];
    e_p(3*i,:)=[t(i,2) t(i,3)];
end
e_p=sort(e_p,2);
e_p=unique(e_p,'stable','rows');
% node to edge
p_e=zeros(size(p,1),7);
for i=1:size(p,1)
    [findp_e,~]=find(e_p==i);
    countp_e=1;
    for j=1:size(findp_e,1)
        p_e(i,countp_e)=findp_e(j);
        countp_e=countp_e+1;
    end
end
% node to triangle
p_t=zeros(size(p,1),7);
for i=1:size(p,1)
    [findp_t,~]=find(t==i);
    countp_t=1;
    for j=1:size(findp_t,1)
        p_t(i,countp_t)=findp_t(j);
        countp_t=countp_t+1;
    end
end
% edge to triangle
e_t=zeros(size(e_p,1),2);
for i=1:size(e_p,1)
    [finde_t,~]=find(t==e_p(i,1));
    counte_t=1;
    for j=1:size(finde_t,1)
        if ismember(e_p(i,2),t(finde_t(j),:))==1
            e_t(i,counte_t)=finde_t(j);
            counte_t=counte_t+1;
        end
    end
end
% triangle to edge
t_e=zeros(size(t,1),3);
for i=1:size(t,1)
    [findt_e,~]=find(e_t==i);
    countt_e=1;
    for j=1:size(findt_e,1)
        t_e(i,countt_e)=findt_e(j);
        countt_e=countt_e+1;
    end
end

end