function dK=PartialGen_new(obj)
%% Initialization
    % num of nodes,elements & constrains
    [nnode,~]=size(obj.meshobj.p);
    [nelem,~]=size(obj.meshobj.t);
    
    % elastic matrix of the matrix material
    Dm = (obj.caseobj.E/(1-obj.caseobj.NU*obj.caseobj.NU))*[1 obj.caseobj.NU 0 ; obj.caseobj.NU 1 0 ; 0 0 (1-obj.caseobj.NU)/2];
    
    % calculate partial and assembly
    dK=cell(nelem,2);
    for i=1:2*nelem
        dK{i}=sparse(2*nnode,2*nnode);
    end
    for k=1:size(obj.Reinforced_Element,1)
        row_t=obj.Reinforced_Element(k,1);
        % shape matrix
        row_p=obj.meshobj.t(row_t,:);
        p1=obj.meshobj.p(row_p(1),:);
        p2=obj.meshobj.p(row_p(2),:);
        p3=obj.meshobj.p(row_p(3),:);
        AB=p2-p1;
        AC=p3-p1;
        At=abs(det([AB;AC]))/2;
        betai=p2(2)-p3(2);
        betaj=p3(2)-p1(2);
        betam=p1(2)-p2(2);
        gammai=p3(1)-p2(1);
        gammaj=p1(1)-p3(1);
        gammam=p2(1)-p1(1);
        B = [betai 0 betaj 0 betam 0 ;
            0 gammai 0 gammaj 0 gammam ;
            gammai betai gammaj betaj gammam betam]/(2*At);
        rho=obj.Reinforced_Element(k,2);
        theta=obj.Reinforced_Element(k,3);
        elem_par=Par_k(At,B,Dm,obj.caseobj.thi,rho,theta,obj.D_fiber,obj.PDc);
        ki=row_p(1);
        kj=row_p(2);
        km=row_p(3);
        DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
        for n=1:2
            dK{k,n}(DOF,DOF) = dK{k,n}(DOF,DOF) + elem_par(:,:,n);
        end
    end

end

%% Element partial derivative
function elem_par=Par_k(At,B,Dm,thi,rho,theta,D_fiber,PDc)
    
    % elastic matrix of the fiber
    Dfiber=D_fiber(theta);
    
    % partial derivative of the elastic matrix of the fiber
    Par_Dc=PDc(theta);
    
    % partial delta K
    elem_par=zeros(6,6,2);
    elem_par(:,:,1)=thi*At*B'*(Dfiber-Dm)*B;
    elem_par(:,:,2)=thi*rho*At*B'*Par_Dc*B;
    
end