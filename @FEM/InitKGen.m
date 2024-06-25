function initK=InitKGen(obj)

% Initialization
    % num of nodes,elements & constrains
    [nnode,~]=size(obj.meshobj.p);
    [nelem,~]=size(obj.meshobj.t);
    initK=sparse(2*nnode,2*nnode);
    
    % plastic elastic matrix
    Dm = (obj.caseobj.E/(1-obj.caseobj.NU*obj.caseobj.NU))*[1 obj.caseobj.NU 0 ; obj.caseobj.NU 1 0 ; 0 0 (1-obj.caseobj.NU)/2];
    
    % generate initial stiffness matrix (no fiber)
    % elemental stiffness of matrix elements
    for i=1:nelem
        ke = KE_M(obj.meshobj.p,obj.meshobj.t,obj.caseobj.thi,Dm,i);
        ki=obj.meshobj.t(i,1);
        kj=obj.meshobj.t(i,2);
        km=obj.meshobj.t(i,3);
        DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
        initK(DOF,DOF)= initK(DOF,DOF)+ke;
    end
    
end

%% Matrix elemental stiffness matrix
function [ke]=KE_M(p,t,thi,Dm,row_count)
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
% element stiffness
ke=thi*B'*(At*Dm)*B;
end