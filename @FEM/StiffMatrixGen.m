 function [K,U,F]=StiffMatrixGen(obj)
    % num of nodes,elements & constrains
    [nnode,~]=size(obj.p);
    [nelem,~]=size(obj.t);
%     [nconstrain,~]=size(obj.constrain);
    % define spaces of U,F,K
    F=sparse(2*nnode,1);
    U=zeros(2*nnode,1);
    K=sparse(2*nnode,2*nnode);
    % stiffness matrix assembly
    t_nof=obj.t; % store the element without fiber
    if isempty(obj.Reinforced_Element)==0 % elemental stiffness of fiber reinforced elements
        for i=1:size(obj.Reinforced_Element,1)
            [k,row_t]=KE_F(obj.thi,obj.t,obj.p,obj.Reinforced_Element,obj.E,obj.NU,obj.Df,i);
            t_nof(row_t,:)=0;
            ki=obj.t(row_t,1);
            kj=obj.t(row_t,2);
            km=obj.t(row_t,3);
            DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
            K(DOF,DOF)= K(DOF,DOF)+k;
        end
    end
    % elemental stiffness of matrix elements
    for i=1:nelem
        row_t=i;
        if t_nof(row_t,1)~=0
            [ke] = KE_M(obj.p,obj.t,obj.thi,obj.E,obj.NU,row_t);
            ki=obj.t(i,1);
            kj=obj.t(i,2);
            km=obj.t(i,3);
            DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
            K(DOF,DOF)= K(DOF,DOF)+ke;
        end
    end
    % displacement boundary condition
%     KOLD=K; % save original K
%     for i=1:nconstrain
%         m=obj.constrain(i,1);
%         n=obj.constrain(i,2);
%         U(2*(m-1)+n)=obj.constrain(i,3);
%         K(2*(m-1)+n,:)=0;
%         K(:,2*(m-1)+n)=0;
%         K(2*(m-1)+n,2*(m-1)+n)=1;
%         F=F-obj.constrain(i,3)*KOLD(:,2*(m-1)+n);
%     end
%     for i=1:nconstrain
%         m=obj.constrain(i,1);
%         n=obj.constrain(i,2);
%         F(2*(m-1)+n)=obj.constrain(i,3);
%     end

    % Load boundary condition
    for i=1:size(obj.load,1)
        j=obj.load(i,1);
        n=obj.load(i,2);
        F(2*(j-1)+n)=obj.load(i,3);
    end
    
    % solve nodal displacement
%     U=lsqminnorm(K,F);
%     K=KOLD;
    U(obj.freedofs,:)=lsqminnorm(K(obj.freedofs,obj.freedofs),F(obj.freedofs,:));
    U(obj.fixeddofs,:)= 0;
    F=K*U;
end

%% Matrix elemental stiffness matrix
function [ke]=KE_M(p,t,thi,E,NU,row_count)
% get triangle vertices
row_p=t(row_count,:);
p1=p(row_p(1),:);
p2=p(row_p(2),:);
p3=p(row_p(3),:);
% area calculation
AB=p2-p1;
AC=p3-p1;
At=abs(det([AB;AC]))/2;
% element stiffness calculation
% matrix
betai=p2(2)-p3(2);
betaj=p3(2)-p1(2);
betam=p1(2)-p2(2);
gammai=p3(1)-p2(1);
gammaj=p1(1)-p3(1);
gammam=p2(1)-p1(1);
B = [betai 0 betaj 0 betam 0 ;
    0 gammai 0 gammaj 0 gammam ;
    gammai betai gammaj betaj gammam betam]/(2*At); % Shape matrix
Dm = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
% element stiffness
ke=thi*B'*(At*Dm)*B;
end

%% Fiber elemental stiffness matrix
function [ke,row_t]=KE_F(thi,t,p,Reinforced_Element,E,NU,Df,row_count)
% get triangle vertices
row_t=Reinforced_Element(row_count,1);
row_p=t(row_t,:);
p1=p(row_p(1),:);
p2=p(row_p(2),:);
p3=p(row_p(3),:);

% get fiber area & direction
Reinforced_info=Reinforced_Element(row_count,:);
Reinforced_info(:,1)=[];
Reinforced_info(Reinforced_info==-10)=[];
Af=zeros(1,size(Reinforced_info,2)/2);
theta=zeros(1,size(Reinforced_info,2)/2);
T=zeros(3,3,size(Reinforced_info,2)/2);

for i=1:size(Reinforced_info,2)/2
    Af(i)=Reinforced_info(1,2*i-1);
    theta(i)=Reinforced_info(1,2*i);
    % rotation matrix
%     theta(i)=theta(i)-0.5*pi; % angle between fiber direction 1 and x axis
    theta(i)=theta(i);
    st = sin(theta(i));
    ct = cos(theta(i));
    T(:,:,i) = [ct^2 st^2 2*st*ct;
                st^2 ct^2 -2*st*ct; 
                -st*ct st*ct (ct^2-st^2)];
end

% area calculation
AB=p2-p1;
AC=p3-p1;
At=abs(det([AB;AC]))/2;

% element stiffness calculation
% matrix
betai=p2(2)-p3(2);
betaj=p3(2)-p1(2);
betam=p1(2)-p2(2);
gammai=p3(1)-p2(1);
gammaj=p1(1)-p3(1);
gammam=p2(1)-p1(1);
B = [betai 0 betaj 0 betam 0 ;
    0 gammai 0 gammaj 0 gammam ;
    gammai betai gammaj betaj gammam betam]/(2*At); % Shape matrix
Dm = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
% fiber
Dc=zeros(3,3,size(Reinforced_info,2)/2);
for i=1:size(Reinforced_info,2)/2
    Dc(:,:,i) = T(:,:,i)'*Df*T(:,:,i);
end

% element stiffness
fiber_area=Af(1);
fiber_D=Af(1)*Dc(:,:,1);
for i=2:size(Reinforced_info,2)/2
    fiber_area=fiber_area+Af(i);
    fiber_D=fiber_D+Af(i)*Dc(:,:,i);
end

ke=thi*B'*(At*Dm-fiber_area*Dm+fiber_D)*B;
end