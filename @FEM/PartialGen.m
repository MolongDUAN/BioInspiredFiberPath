function dKdPc=PartialGen(obj,m,En,cp,C_prob,C_der,D_fiber,PDc,Theta,Theta_der)
% function dKdPc=PartialGen(obj,p_t,m,En,cp,C_prob,C_der,D_fiber,PDc,Theta,Theta_der)
    
    addpath(genpath('nurbs-1.3.7'))
    % num of nodes,elements & constrains
    [nnode,~]=size(obj.p);
    % define spaces of partial K
    dKdPc=cell(size(cp,2),1);
    for j=1:size(cp,2)
        n=2*size(cp{j},1);
        for i=1:n
            dKdPc{j,i}=sparse(2*nnode,2*nnode);
        end
    end
    % numerical integration step size
    delta_u=1/(En-1);
    
    % reinforced element for every spline
    p=parallel.pool.Constant(obj.p);
    t=parallel.pool.Constant(obj.t);
    w=parallel.pool.Constant(obj.w);
    parfor j=1:size(cp,2)
        % b-spline
        n_p=size(cp{j},1);   % number of control points
        % Basis matrix N generation
        N=BasisGen(En-1,n_p-1,m);
        Nd=diff(N,1,1)/delta_u;
        Nd(En,:)=(N(En,:)-N(En-1,:))/delta_u;
        % spline points
        bsp=N*cp{j};
        % identify the reinforced elements
        [partial_element{j,1},theta_der{j}]=FiberElemIdentify(p.Value,t.Value,w.Value,cp{j},En,delta_u,Nd,bsp,C_prob,Theta,Theta_der);
%         [partial_element{j,1},theta_der{j}]=FiberElemIdentify_New(p.Value,t.Value,p_t,w.Value,cp{j},m,En,n_p,delta_u,Nd,bsp,C_prob,Theta,Theta_der);
    end
    
    % elastic matrix of the matrix material
    Dm = (obj.E/(1-obj.NU*obj.NU))*[1 obj.NU 0 ; obj.NU 1 0 ; 0 0 (1-obj.NU)/2];
    
    % partial for every spline
    for j=1:size(cp,2)
        % b-spline
        cp_forder=[cp{j}(:,1);cp{j}(:,2)];
        n_p=size(cp{j},1);   % number of control points
        % Basis matrix N generation
        N=BasisGen(En-1,n_p-1,m);
        % spline points
        Nd=diff(N,1,1)/delta_u;
        Nd(En,:)=(N(En,:)-N(En-1,:))/delta_u;

        % calculate partial and assembly
        for i=1:size(partial_element{j},1)
            row_t=partial_element{j}(i,1);
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
            Lf=partial_element{j}(i,2)/obj.w;
            theta=partial_element{j}(i,3);
            elem_par=Par_k(obj.p,obj.t,row_t,B,Dm,obj.thi,obj.w,cp{j},En,cp_forder,delta_u,n_p,N,Nd,Lf,theta,theta_der{j}(i,:),C_prob,C_der,D_fiber,PDc);
            ki=obj.t(row_t,1);
            kj=obj.t(row_t,2);
            km=obj.t(row_t,3);
            DOF = [2*ki-1;2*ki;2*kj-1;2*kj;2*km-1;2*km];
            for n=1:2*size(cp{j},1)
                dKdPc{j,n}(DOF,DOF) = dKdPc{j,n}(DOF,DOF) + elem_par(:,:,n);
            end
%             Elem_Par{i}=elem_par;
        end
    end

end

%% Element partial derivative
function elem_par=Par_k(p,t,row_t,B,Dm,thi,w,cp,En,cp_forder,delta_u,n_p,N,Nd,Lf,theta,theta_der,C_prob,C_der,D_fiber,PDc)
    % element
    triangle=t(row_t,:);
    vertices=p(triangle,:);
    % Ensure the vertices order is clockwise (increases program runtime)
    Triangle=polyshape(flipud(vertices));
    tp1=Triangle.Vertices(1,:);
    tp2=Triangle.Vertices(2,:);
    tp3=Triangle.Vertices(3,:);

    % B-spline
    bsp=N*cp;
    
    % probability and its derivative calculation
    num_p=1:En;
    Prob=C_prob(tp1(1),tp2(1),tp3(1),tp1(2),tp2(2),tp3(2),bsp(num_p,1),bsp(num_p,2));
    Prob_der=C_der(tp1(1),tp2(1),tp3(1),tp1(2),tp2(2),tp3(2),bsp(num_p,1),bsp(num_p,2));
    
    % elastic matrix of the fiber
    Dfiber=D_fiber(theta);
    
    % partial derivative of the elastic matrix of the fiber
    Par_Dc=PDc(theta);
    Par_TheA=theta_der;
    % numerical integration
    parfor i=1:En
        N_forder=[N(i,:),zeros(1,n_p);zeros(1,n_p),N(i,:)];
        Nd_forder=[Nd(i,:),zeros(1,n_p);zeros(1,n_p),Nd(i,:)];
        % Par_A
        par_A_1=Nd_forder*cp_forder*Prob_der(i,:)*N_forder;
        par_A_2=Prob(i,:)*Nd_forder;
        par_A(:,:,i)=(par_A_1+par_A_2)*delta_u;
        % partial derivative of the fiber length
        par_l=sqrt((Nd(i,:)*cp)*(Nd(i,:)*cp)');
        par_Lf_1=Prob_der(i,:)*N_forder*par_l;
        par_Lf_2=Prob(i,:)*(1/(2*par_l))*(Nd_forder'*Nd_forder+(Nd_forder'*Nd_forder)')*cp_forder;
        par_Lf(:,:,i)=(par_Lf_1+par_Lf_2')*delta_u;
    end
    Par_A=sum(par_A,3);
    Par_theta=Par_TheA*Par_A;
    Par_Dfiber=reshape(Par_Dc,1,[]).*Par_theta(:);
    Par_Dfiber=Par_Dfiber';
    Par_Dfiber=reshape(Par_Dfiber,3,3,2*n_p);
    Par_Lf=sum(par_Lf,3);
    
    % partial delta K
    elem_par=zeros(6,6,2*n_p);
    for i=1:2*n_p
        elem_par(:,:,i)=thi*w*(Par_Lf(i)*B'*(Dfiber-Dm)*B+Lf*B'*Par_Dfiber(:,:,i)*B);
    end
    
end