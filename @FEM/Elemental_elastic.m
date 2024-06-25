 function ElemD=Elemental_elastic(obj,D_fiber,initElemD)
%% Initialization
    % num of nodes,elements
    [nelem,~]=size(obj.t);
    % define spaces of U,F,K
    ElemD=cell(nelem,1);
    % plastic elastic matrix
    Dm = (obj.E/(1-obj.NU*obj.NU))*[1 obj.NU 0 ; obj.NU 1 0 ; 0 0 (1-obj.NU)/2];
    
    % generate initial stiffness matrix (no fiber)
    for i=1:nelem
        if isempty(initElemD{i})==1
            ElemD{i} = Dm;
        else
            ElemD{i}=initElemD{i};
        end
    end
    
%% Fiber reinforcement to the stiffness matrix
    % stiffness matrix assembly
    if isempty(obj.Reinforced_Element)==0 % elemental stiffness of fiber reinforced elements
        for i=1:size(obj.Reinforced_Element,1)
            [dD,row_t]=KE_F(obj.t,obj.p,obj.Reinforced_Element,Dm,D_fiber,i);
            ElemD{row_t}=ElemD{row_t}+dD;
        end
    end
    
end

%% Fiber elemental stiffness matrix
function [dD,row_t]=KE_F(t,p,Reinforced_Element,Dm,D_fiber,row_count)
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

% fiber
Dfiber=D_fiber(theta);

% increment of element stiffness
dD=(Af/At)*(Dfiber-Dm);
end