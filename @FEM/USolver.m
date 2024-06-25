function U=USolver(obj)

    [nnode,~]=size(obj.meshobj.p);
    F=zeros(2*nnode,1);
    U=zeros(2*nnode,1);

% Solve KU=F
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
% 
    % Load boundary condition
    for i=1:size(obj.load,1)
        j=obj.load(i,1);
        n=obj.load(i,2);
        F(2*(j-1)+n)=obj.load(i,3);
    end
    
    % solve nodal displacement
%     U=lsqminnorm(K,F);
%     K=KOLD;
    U(obj.freedofs,:)=lsqminnorm(obj.K(obj.freedofs,obj.freedofs),F(obj.freedofs,:));
    U(obj.fixeddofs,:)= 0;
%     F=obj.K*U;

end