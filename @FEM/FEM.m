classdef FEM < handle
    % FEM solver
    properties
        caseobj % defined case
        meshobj % caculated mesh
        load
        fixeddofs
        freedofs
        D_fiber % symbolic elastic matrix of fiber
        PDc % symbolic partial with respect to theta
        Reinforced_Element % elements reinforced by fiber, [element number, fiber area, fiber direction]
    end
    
    properties(Dependent)
        K % stiffness matrix
        U % displacement
        F % force
        PartStiff % stiffness of the part
        dK % partial stiffness matrix
        ElemStrain % elemental strain vector
    end
    
    properties(Access=private)
        initK=[]; % stiffness of matrix
        K_cache=[];
        U_cache=[];
        F_cache=[];
    end
    
    methods
        function obj=FEM(caseobj,meshobj,load,fixeddofs,freedofs,D_fiber,PDc,Reinforced_Element)
            % initialization
            obj.caseobj=caseobj;
            obj.meshobj=meshobj;
            obj.load=load;
            obj.fixeddofs=fixeddofs;
            obj.freedofs=freedofs;
            obj.D_fiber=D_fiber;
            obj.PDc=PDc;
            obj.Reinforced_Element=Reinforced_Element;
        end
        
        % update Reinforced_Element
        function set.Reinforced_Element(obj,value)
            obj.Reinforced_Element=value;
            obj.K_cache=[];
            obj.U_cache=[];
            obj.F_cache=[];
        end
        
        % Solve FE problem
        function initK=get.initK(obj)
            if isempty(obj.initK)
                initK=InitKGen(obj);
            else
                initK=obj.initK;
            end
        end
        
        function K=get.K(obj)
            if isempty(obj.K_cache)
                K=StiffMatrixGen_new(obj);
                obj.K_cache=K;
            else
                K=obj.K_cache;
            end
        end
        
        function U=get.U(obj)
            if isempty(obj.U_cache)
                U=USolver(obj);
                obj.U_cache=U;
            else
                U=obj.U_cache;
            end
        end
        
        function F=get.F(obj)
            if isempty(obj.F_cache)
                F=obj.K*obj.U;
                obj.F_cache=F;
            else
                F=obj.F_cache;
            end
        end
        
        % Solve partial K
        function dK=get.dK(obj)
            dK=PartialGen_new(obj);
        end
        
        % Solve elemental strain
        function ElemStrain=get.ElemStrain(obj)
            ElemStrain=Elemental_strain(obj);
        end

%         % Solve elemental elastic matrix
%         ElemD=Elemental_elastic(obj,D_fiber,initElemD);
%         
%         % Solve elemental stress
%         ElemStress=Elemental_stress(obj,ElemD,U);
        
        function PartStiff=get.PartStiff(obj)
            % Specimen stiffness output
            % load boundary condition
            Force=sum(obj.load(:,3));
            [nload,~]=size(obj.load);
            displacement=0;
            for i=1:nload
                j=obj.load(i,1);
                n=obj.load(i,2);
                displacement=displacement+obj.U(2*(j-1)+n,:);
            end
            displacement=displacement/nload;
            PartStiff=Force/displacement;
        end
        
    end
end