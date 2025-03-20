%----------------------------------------------------------------------------
%   plane stress analysis of a solid using linear quadratic elements           
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   disp = system nodal displacement vector
%   eldisp = element nodal displacement vector
%   stress = matrix containing stresses
%   strain = matrix containing strains
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

[gcoord, nodes] = abaqusToMATLAB('Job-2.inp');

nel=length(nodes);                   % number of elements
nnel=4;                  % number of nodes per element
ndof=2;                  % number of dofs per node
nnode=length(gcoord);                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=2e11;        % elastic modulus
poisson=0.3;             % Poisson's ratio

%---------------------------------------------
%  Gather nodes from the left and right side 
%  of the beam
%---------------------------------------------
min_x = min(gcoord(:,1));  % Left side
max_x = max(gcoord(:,1));  % Right side

% Find the nodes on the left (x == min_x)
left_nodes = find(gcoord(:,1) == min_x);

% Find the nodes on the right (x == max_x)
right_nodes = find(gcoord(:,1) == max_x);

left_dofs = getDOFs(left_nodes);
right_dofs = getDOFs(right_nodes);

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof= left_dofs;        % first three dofs are constrained
bcval=zeros(size(bcdof));        % whose described values are 0 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nel,3);    % matrix containing stress components
strain=zeros(nel,3);    % matrix containing strain components
index=zeros(edof,1);    % index vector
Bmtx2=zeros(3,edof);   % kinematic matrix
matmtx=zeros(3,3);      % constitutive matrix

%----------------------------
%  force vector
%----------------------------
numBoundaryNodes = length(right_nodes);
elementLength = abs(gcoord(right_nodes(2), 2) - gcoord(right_nodes(1), 2));
tractionForcePerNode = -2000 / numBoundaryNodes;

for i = 1:numBoundaryNodes
    node = right_nodes(i);
    dof_y = 2*node;  % y-displacement DOF for the node
    ff(dof_y) = ff(dof_y) + tractionForcePerNode;  % Apply traction force in the y-direction
end

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

matmtx=fematiso(1,emodule,poisson);        % compute constitutive matrix

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element
nd(4)=nodes(iel,4); % 4th connected node for (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node
x4=gcoord(nd(4),1); y4=gcoord(nd(4),2);% coord values of 4th node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%-------------------------------------------------------
%  find the derivatives of shape functions
%-------------------------------------------------------

gaussPoints = [-1/sqrt(3),  1/sqrt(3);
               -1/sqrt(3),  1/sqrt(3)];
weights = [1, 1];
k = zeros(edof, edof);
for i=1:2
        for j=1:2
            xi = gaussPoints(1, i);  % xi coordinate
            eta = gaussPoints(2, j); % eta coordinate
            weight = weights(i) * weights(j);
            
            % Compute shape function derivatives w.r.t x, y at this Gauss point
            [dhdx, dhdy, detJ] = q4_shape_function_derivatives(x1, y1, x2, y2, x3, y3, x4, y4, xi, eta);
            
            % Compute the kinematic matrix (B-matrix)
            Bmtx2 = fekine2d(nnel, dhdx, dhdy);

            k = k + (Bmtx2' * matmtx * Bmtx2 * detJ * weight);
        end
end

%k=Bmtx2'*matmtx*Bmtx2*area;      % element stiffnes matrix

kk=feasmbl1(kk,k,index);  % assemble element matrices 

end

%-----------------------------
%   apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------

disp=kk\ff;   

%---------------------------------------
%  element stress computation
%---------------------------------------

for ielp=1:nel           % loop for the total number of elements

nd(1)=nodes(ielp,1); % 1st connected node for (iel)-th element
nd(2)=nodes(ielp,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(ielp,3); % 3rd connected node for (iel)-th element
nd(4)=nodes(ielp,4); % 4th connected node for (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node
x4=gcoord(nd(4),1); y4=gcoord(nd(4),2);% coord values of 4th node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------

for i=1:edof
eldisp(i)=disp(index(i));
end

for i=1:2
        for j=1:2
            xi = gaussPoints(1, i);  % xi coordinate
            eta = gaussPoints(2, j); % eta coordinate
            weight = weights(i) * weights(j);
            
            % Compute shape function derivatives w.r.t x, y at this Gauss point
            [dhdx, dhdy, detJ] = q4_shape_function_derivatives(x1, y1, x2, y2, x3, y3, x4, y4, xi, eta);
            
            % Compute the kinematic matrix (B-matrix)
            Bmtx2 = fekine2d(nnel, dhdx, dhdy);
        end
end

estrain=Bmtx2*eldisp;             % compute strains
estress=matmtx*estrain;             % compute stresses

for i=1:3
strain(ielp,i)=estrain(i);          % store for each element
stress(ielp,i)=estress(i);          % store for each element          
end

end

%------------------------------------
% print fem solutions
%------------------------------------


disp_x = disp(1:2:end);
disp_y = disp(2:2:end);
% Find the maximum displacement
[maxDisp, maxDispIdx] = max(abs(disp));  % Get the maximum displacement value and its index
maxDispNode = ceil(maxDispIdx / 2);      % Find the node corresponding to the max displacement
maxDispDOF = mod(maxDispIdx, 2);         % Determine whether it's x or y displacement
if maxDispDOF == 0
    maxDispDirection = 'y';
else
    maxDispDirection = 'x';
end

maxDispx = max(abs(disp_x));
fprintf('Maximum displacement at x: %.5e m\n', maxDispx);
maxDispy = max(abs(disp_y));
fprintf('Maximum displacement at y: %.5e m\n', -maxDispy);

% Print the maximum displacement information
fprintf('Maximum displacement: %.5e m at Node %d in the %s direction\n', maxDisp, maxDispNode, maxDispDirection);

% Find the maximum stress for each element
maxStress = max(sqrt(stress(:,1).^2 + stress(:,2).^2 + stress(:,3).^2));  % Max von Mises stress across all elements

% Print the maximum stress
fprintf('Maximum stress: %.5f Pa\n', maxStress);

% Find the x-coordinate corresponding to L/2
L_half = 1 / 2;  % Define L if not already defined
tolerance = 1;
% Find the nodes at x = L/2
nodes_at_L_half = find(abs(gcoord(:,1) - L_half) < tolerance);  % Set an appropriate tolerance for precision

% Extract the stress for the corresponding elements
stress_at_L_half = [];
for i = 1:length(nodes_at_L_half)
    node = nodes_at_L_half(i);
    
    % Identify the elements that contain this node
    for elem_idx = 1:nel
        if any(nodes(elem_idx, :) == node)
            stress_at_L_half = [stress_at_L_half; stress(elem_idx, :)];
        end
    end
end

[analytical_ux, analytical_uy, analytical_stressx] = analytical_solution();
fprintf('Analytical ux: %.5e m\n', analytical_ux);
fprintf('Analytical uy: %.5e m\n', analytical_uy);
fprintf('Analytical Sxx: %.5e m\n', analytical_stressx);


% Calculate von Mises stress or use other criteria
vonMisesStress = sqrt(stress_at_L_half(:,1).^2 + stress_at_L_half(:,2).^2 - stress_at_L_half(:,1).*stress_at_L_half(:,2) + stress_at_L_half(:,3).^2);

% Find the maximum stress
maxStress = max(vonMisesStress);
fprintf('Maximum von Mises stress at x = L/2: %.5f Pa\n', maxStress);


function dofs = getDOFs(nodes)
    dofs = [];
    for i = 1:length(nodes)
        node = nodes(i);
        dofs = [dofs, (2*node-1), (2*node)];  % x- and y-displacement DOFs
    end
end

function [dhdx, dhdy, detJ] = q4_shape_function_derivatives(x1, y1, x2, y2, x3, y3, x4, y4, xi, eta)
    % Compute shape function derivatives for Q4 element at a given Gauss point (xi, eta)

    % Shape function derivatives w.r.t. natural coordinates (xi, eta)
    dN_dxi = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
    dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];
    
    % Jacobian matrix
    J = [dN_dxi * [x1; x2; x3; x4], dN_deta * [x1; x2; x3; x4];
         dN_dxi * [y1; y2; y3; y4], dN_deta * [y1; y2; y3; y4]];

    % Determinant of the Jacobian
    detJ = det(J);
    
    % Inverse of the Jacobian matrix
    invJ = inv(J);

    % Derivatives of shape functions w.r.t. x and y
    dhdx = invJ(1,1) * dN_dxi + invJ(1,2) * dN_deta;
    dhdy = invJ(2,1) * dN_dxi + invJ(2,2) * dN_deta;
end


function [ux, uy, stress_x] = analytical_solution

    E = 200e9;
    L = 10;
    D = 1;
    v = 0.3; % Poisson Ratio
    I = D^3/12;
    P = 2000/D;
    x = L;
    y = D/2;
    ux = (P*y/(6*E*I))*(((6*L-3*x)*x)+(2+v)*(y^2-(0.25*D^2)));
    uy = -(P/(6*E*I))*(((L-x)*3*v*y^2)+((4+5*v)*0.25*x*D^2)+((3*L-x)*(x^2)));

    x = L/2;
    stress_x = P*(L-x)*y/I;
end










