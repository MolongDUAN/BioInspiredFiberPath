 function dK=ElemIncrementK(obj,D_fiber)
%% Initialization
    % plastic elastic matrix
    Dm = (obj.E/(1-obj.NU*obj.NU))*[1 obj.NU 0 ; obj.NU 1 0 ; 0 0 (1-obj.NU)/2];
    
%% Fiber reinforcement to the stiffness matrix
    % stiffness matrix assembly
    if isempty(obj.Reinforced_Element)==0 % elemental stiffness of fiber reinforced elements
        for i=1:size(obj.Reinforced_Element,1)
            [dk,~]=KE_F(obj.thi,obj.t,obj.p,obj.Reinforced_Element,Dm,D_fiber,i);
            dK{i}=dk;
        end
    end
end

%% Fiber elemental stiffness matrix
function [dk,row_t]=KE_F(thi,t,p,Reinforced_Element,Dm,D_fiber,row_count)
% get triangle vertices
row_t=Reinforced_Element(row_count,1);
row_p=t(row_t,:);
p1=p(row_p(1),:);
p2=p(row_p(2),:);
p3=p(row_p(3),:);

% get fiber area & direction
Af=Reinforced_Element(row_count,2);
theta=Reinforced_Element(row_count,3);

% element stiffness calculation
% area calculation
AB=p2-p1;
AC=p3-p1;
At=abs(det([AB;AC]))/2;
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

% fiber
Dfiber=D_fiber(theta);

% increment of element stiffness
dk=thi*Af*B'*(Dfiber-Dm)*B;
end