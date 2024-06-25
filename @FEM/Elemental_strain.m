function ElemStrain=Elemental_strain(obj)

    % num of nodes,elements
    [nelem,~]=size(obj.meshobj.t);
    
    ElemStrain=zeros(nelem,3);
    
    parfor i=1:nelem
        % get nodal displacement
        l=obj.meshobj.t(i,1);m=obj.meshobj.t(i,2);n=obj.meshobj.t(i,3);
        u=[obj.U(2*l-1),obj.U(2*l),obj.U(2*m-1),obj.U(2*m),obj.U(2*n-1),obj.U(2*n)]';
        % get triangle vertices
        row_p=obj.meshobj.t(i,:);
        p1=obj.meshobj.p(row_p(1),:);
        p2=obj.meshobj.p(row_p(2),:);
        p3=obj.meshobj.p(row_p(3),:);
        % get shape matrix
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
        % elemental stress
        strain_tmp=B*u;
        ElemStrain(i,:)=strain_tmp';
    end

end