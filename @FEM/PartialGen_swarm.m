function dKdPc=PartialGen_swarm(obj,delta_u,bsp,bsp_paru,par_bsp,par_bsp_paru,swarm_l,Reinforced_element,Theta_der,C_prob,C_der,D_fiber,PDc)

% Initialization
    % num of nodes,elements & constrains
    [nnode,~]=size(obj.p);
    % elastic matrix of the matrix material
    Dm = (obj.E/(1-obj.NU*obj.NU))*[1 obj.NU 0 ; obj.NU 1 0 ; 0 0 (1-obj.NU)/2];
    
% partial
    dKdPc_tmp=cell(5,1);
    for i=1:5
        dKdPc_tmp{i,1}=sparse(2*nnode,2*nnode);
    end
    % reinforced elements in local convex hull
    hull_element=Reinforced_element;
    theta_der=Theta_der;
    for k=1:size(hull_element,1)
        row_t=hull_element(k,1);
        % shape matrix
        row_p=obj.t(row_t,:);
        p1=obj.p(row_p(1),:);
        p2=obj.p(row_p(2),:);
        p3=obj.p(row_p(3),:);
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
        Lf=hull_element(k,2);
        theta=hull_element(k,3);
        elem_par=Par_k(obj.p,obj.t,row_t,B,Dm,obj.thi,obj.w,delta_u,bsp,bsp_paru,par_bsp,par_bsp_paru,swarm_l,Lf,theta,theta_der(k,:),C_prob,C_der,D_fiber,PDc);
        ki=obj.t(row_t,1);
        kj=obj.t(row_t,2);
        km=obj.t(row_t,3);
        DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
        for n=1:5
            dKdPc_tmp{n,1}(DOF,DOF) = dKdPc_tmp{n,1}(DOF,DOF) + elem_par(:,:,n);
        end
    end
    dKdPc=dKdPc_tmp;

end

%% Element partial derivative
function elem_par=Par_k(p,t,row_t,B,Dm,thi,w,delta_u,bsp,bsp_paru,par_bsp,par_bsp_paru,swarm_l,Lf,theta,theta_der,C_prob,C_der,D_fiber,PDc)
    % element
    triangle=t(row_t,:);
    vertices=p(triangle,:);
    % Ensure the vertices order is clockwise
    tp1=vertices(3,:);
    tp2=vertices(2,:);
    tp3=vertices(1,:);
    
    % probability and its derivative calculation
    Prob=C_prob(tp1(1),tp2(1),tp3(1),tp1(2),tp2(2),tp3(2),bsp(:,1),bsp(:,2));
    Prob_der=C_der(tp1(1),tp2(1),tp3(1),tp1(2),tp2(2),tp3(2),bsp(:,1),bsp(:,2));
    
    % elastic matrix of the fiber
    Dfiber=D_fiber(theta);
    
    % partial derivative of the elastic matrix of the fiber
    Par_Dc=PDc(theta);
    Par_TheA=theta_der;
    % numerical integration
    parfor i=1:size(Prob,1)
        % Par_A
        par_A_1=bsp_paru(i,:)'*Prob_der(i,:)*par_bsp(:,:,i);
        par_A_2=Prob(i,:)*par_bsp_paru(:,:,i);
        par_A(:,:,i)=(par_A_1+par_A_2)*delta_u;
        % partial derivative of the fiber length
        par_Lf(:,:,i)=Prob_der(i,:)*par_bsp(:,:,i)*swarm_l*delta_u;
    end
    Par_A=sum(par_A,3);
    Par_theta=Par_TheA*Par_A;
    Par_Dfiber=reshape(Par_Dc,1,[]).*Par_theta(:);
    Par_Dfiber=Par_Dfiber';
    Par_Dfiber=reshape(Par_Dfiber,3,3,5);
    Par_Lf=sum(par_Lf,3);
    
    % partial delta K
    elem_par=zeros(6,6,5);
    for i=1:5
        elem_par(:,:,i)=thi*w*(Par_Lf(i)*B'*(Dfiber-Dm)*B+Lf*B'*Par_Dfiber(:,:,i)*B);
    end
    
end