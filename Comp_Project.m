% Composite Materials Project 

clc
clear
clear all
format long

%Given Parameters

a = 0.4;                      % An Edge-Length of The Square Plate [m]
%NM = [50e3;-50e3;1e3;-2;7;1]; % Load [N_x;N_y;N_xy;M_x;M_y;M_xy] [N/m;(N*m)m]
NM=[1000;1000;0;0;0;0];
    plyno = 6;                    % Ply Number of The Symmetric Laminate
t = 5e-3;                  % Ply Thickness [m]
ns = 2;                       % Factor of Safety

% The stresses in individual plies must be lower than the ply strength given in table.
% Magnitudes of global mid-plane strains must be lower than 5*10^-3 m/m.
% Magnitudes of global mid-plane curvatures must be lower than 4 m^-1.

% Chooseen Materials and Ply Orientations for A Symmetric Stacking Sequence for 6 Plies are 
% [0*GY-70Epoxy / 90*KevlarEpoxy / 0*GY-70Epoxy]_s

% Properties of Choosen Materials

E_1 (1:6) =181e9; %[262, 76, 262, 262, 76, 262] * 10^9;  % Longitudinal Elastic Modulus [Pa]
E_2 (1:6) = 10.3e9; %[8.3, 5.5, 8.3, 8.3, 5.5, 8.3] * 10^9;  % Transverse Elastic Modulus [Pa]
G_12 (1:6) = 7.17e9;%[4.1, 2.1 ,4.1, 4.1, 2.1, 4.1] *10^9;  % Shear Modulus [Pa]
v_12 (1:6) = 0.28;%[.25, .34, .25, .25, .34, .25];  % Major Poisson's Ratio

Strength_Ultimate_11 = [586, 1379, 586, 586, 1379, 586] * 10^6;  % Ultimate Longitudinal Strength [Pa]
Strength_Ultimate_22 = [41, 28, 41, 41, 28, 41] *10^6;  % Ultimate Transverse Strength [Pa]
Strength_Ultimate_12 = [97, 60, 97, 97, 60, 97] *10^6;  % Ultimate In-plane Shear Strength [Pa]
Density = [1690, 1380, 1690, 1690, 1380, 1690]; % Densities of Materials [kg/m^3]
angle = [-45, 45, 90, 0, -45, 45];  % Choosen Ply Orientations for A Symmetric Stacking Sequence for 6 Plies

% Initial Conditions

Mass = 0;
no_fail=0;
S = zeros(3,3); % Reducing Compliance Matrix
R = zeros(3,3); % Reuter Matrix
R(1,1) = 1;
R(2,2) = 1;
R(3,3) = 2;
A = zeros(3);  % A Matrix
B = zeros(3);  % B Matrix
D = zeros(3);  % D Matrix

h = plyno*t;                %Total Thickness of Laminate [m]
h_coordinate(1) = -h/2;     % [m]
Volume = plyno*t*a^2;       %Total Volume [m^3]

    for n = 1:plyno

        h_coordinate(n+1) = h_coordinate(n)+t;     % Locations of the Ply Surfaces [m]

        % Compliance Matrix
        S(1,1) = 1/E_1(n);
        S(2,2) = 1/E_2(n);
        S(3,3) = 1/G_12(n);
        S(1,2) = -v_12(n)/E_1(n);
        S(2,1) = S(1,2);
        
        v_21(n) = v_12(n)/E_1(n)*E_2(n);    % Minor Poisson's Ratio

        Q = inv(S);   % Reduced Stifness Matrix

        c(n) = cosd(angle(n));  
        s(n) = sind(angle(n));
        T=[c(n)^2 s(n)^2 2*s(n)*c(n); s(n)^2 c(n)^2 -2*s(n)*c(n); -s(n)*c(n) s(n)*c(n) c(n)^2-s(n)^2];   % Transfer Matrix

        Qbar = inv(T)*Q*R*T*inv(R);    % Transfer Reduced Stiffnes Matrix

        A = A+Qbar*(h_coordinate(n+1)-h_coordinate(n));            % A Matrix
        B = B+1/2*(Qbar*(h_coordinate(n+1)^2-h_coordinate(n)^2));  % B Matrix
        D = D+1/3*(Qbar*(h_coordinate(n+1)^3-h_coordinate(n)^3));  % D Matrix

        % ABBD Matrix

        ABBD(1:3,1:3) = A;   % Extansional Stiffness Matrix [Pa*m]
        ABBD(1:3,4:6) = B;   % Coupling Stiffness Matrix [Pa*m^2]
        ABBD(4:6,1:3) = B;
        ABBD(4:6,4:6) = D;   % Bending Stiffness Matrix [Pa*m^3]

        e0 = ABBD\NM;    %[Epsilon*0_x;Epsilon*0_x;Gamma*0_xy ; K_x;K_y;K_xy] [m/m; 1/m]
        Strain_midplane = e0(1:3,1);     % Midplane Strain [m/m]
        Curvature_midplane = e0(4:6,1);  % Midplane Curvatures [1/m]

        Mass = Mass+a^2*t*(Density(n));  % Mass of Laminate [kg]

    end

Density_laminate = Mass/Volume;  % Total Density [kg/m^3]

Astar = inv(A);  % Extansional Compliance Matrix [1/(Pa*m)]
Dstar = inv(D);  % Inverse of the Bending Stiffness Matrix [1/(Pa*m^3)]

E_x = 1/(h*Astar(1,1));   % Effective In-Plane Longitudinal Modulus [Pa]
E_y = 1/(h*Astar(2,2));   % Effective In-Plane Tranverse Modulus [Pa]
G_xy = 1/(h*Astar(3,3));  % Effective In-Plane Shear Modulus [Pa]
v_xy = -(Astar(1,2)/Astar(1,1));   % Effective In-Plane Poisson's Ratio v_xy
v_yx = -(Astar(1,2)/Astar(2,2));   % Effective In-Plane Poisson's Ratio v_yx

Constant_Inplane = [E_x,E_y,G_xy,v_xy,v_yx];    % In-Plane Engineering Constant of Laminate
Specific_Modulus_Inplane = Constant_Inplane(1:3)*(1/Density_laminate);   % Specific In-Plane Modulus [Pa*kg/m^3]

E_xf = 12/(h^3*Dstar(1,1));    % Effective Flexural Longitudinal Modulus [Pa]
E_yf = 12/(h^3*Dstar(2,2));    % Effective Flexural Transverse Modulus [Pa]
G_xyf = 12/(h^3*Dstar(3,3));   % Effective Flexural Shear Modulus [Pa]
v_xyf = -(Dstar(1,2)/Dstar(1,1));   % Effective Flexural Poisson's Ratio v_xyf
v_yxf = -(Dstar(1,2)/Dstar(2,2));   % Effective Flexural Poisson's Ratio v_yxf

Constant_Flexural = [E_xf,E_yf,G_xyf,v_xyf,v_yxf];   % Flexural Engineering Constant of Laminate

    for n = 1:plyno

        z=-plyno*t/2:t/2:plyno*t/2;

        Q=inv(S); 
        c(n)=cosd(angle(n));
        s(n)=sind(angle(n));
        T=[c(n)^2 s(n)^2 2*s(n)*c(n); s(n)^2 c(n)^2 -2*s(n)*c(n); -s(n)*c(n) s(n)*c(n) c(n)^2-s(n)^2]; %Transfer Matrix
        Qbar=inv(T)*Q*R*T*inv(R);
        
        Strain_global(1:3,3*n-2)=Strain_midplane + z (2*n-1) * Curvature_midplane;   % Global Strain at Top       [m/m]
        Strain_global(1:3,3*n-1)=Strain_midplane + z (2*n) * Curvature_midplane;     % Global Strain at Middle    [m/m]
        Strain_global(1:3,3*n)=Strain_midplane + z (2*n+1) * Curvature_midplane;     % Global Strain at Bottom    [m/m]
    
        SG=Strain_global;
        SG(3,:)=Strain_global(3,:)/2;
        SG1(:,3*n-2:3*n) = T * SG(:,3*n-2:3*n) ;
        Strain_local=SG1;
        Strain_local (3,:) = 2 * SG1(3,:);    % Local Strains   [m/m]
    
        Stress_global(:,3*n-2:3*n) = Qbar * Strain_global (:,3*n-2:3*n);      % Global Stresses  [Pa]
    
        Stress_local(:,3*n-2:3*n) = T * Stress_global(:,3*n-2:3*n);           % Local Stresses   [Pa]
      
    end
    
    Stress_global_x=Stress_global(1,:);         % Global Stresses at X Direction [N]
    Stress_global_y=Stress_global(2,:);         % Global Stresses at Y Direction [N]
    Stress_global_xy=Stress_global(3,:);        % Global Stresses at Z Direction [N]
    
    Strain_global_x=Strain_global(1,:);         % Global Strains at X Direction [m/m]
    Strain_global_y=Strain_global(2,:);         % Global Strains at Y Direction [m/m] 
    Strain_global_xy=Strain_global(3,:);        % Global Strains at Z Direction [m/m]

    Stress_local_1=Stress_local(1,:);           % Local Stresses at Direction 1 [N]
    Stress_local_2=Stress_local(2,:);           % Local Stresses at Direction 2 [N]
    Stress_local_12=Stress_local(3,:);          % Local Stresses at Direction 3 [N]

    Strain_local_1=Strain_local(1,:);           % Local Strains at Direction 1 [m/m]
    Strain_local_2=Strain_local(2,:);           % Local Strains at Direction 2 [m/m]
    Strain_local_12=Strain_local(3,:);          % Local Strains at Direction 3 [m/m]

% Failure Test with Safety Factor = 2

for n=1:plyno
    
    %Tsai-Hill Failure Theory
    S1=max(abs(Stress_local_1(1,3*n-2:3*n)));    % Max Stress at Direction 1 for Each Individual Ply [Pa]
    X=Strength_Ultimate_11(n)/ns;                % Ultimate Longitudinal Strength with Factor of Safety [Pa]

    S2=max(abs(Stress_local_2(1,3*n-2:3*n)));    % Max Stress at Direction 2 for Each Individual Ply [Pa]
    Y=Strength_Ultimate_22(n)/ns;                % Ultimate Transverse Strength with Factor of Safety [Pa]

    T12=max(abs(Stress_local_12(1,3*n-2:3*n)));  % Max Stress at Direction 12 for Each Individual Ply [Pa]
    Z=Strength_Ultimate_12(n)/ns;                % Ultimate In-plane Shear Strength with Factor of Safety [Pa]

    tsai=(S1/X)^2-((S1*S2)/(X^2))+(S2/Y)^2+(T12/Z)^2;   % Must be lower than 1 according to Tsai-Hill Failure Theory

    if tsai < 1  % Equation of Tsai-Hill Failure Theory

        no_fail=no_fail+1;
        
    % If value of nonfail is equal to 6 all plies in laminate pass the stress test
    
    end
end

    % When Top, Middle and Bottom Locations of The Surfaces of Plies are Considered Together
    Th=[z(1),z(2),z(3),z(3),z(4),z(5),z(5),z(6),z(7),z(7),z(8),z(9),z(9),z(10),z(11),z(11),z(12),z(13)]; 
    
              
    %%%----- PLOTTING -----%%%
           
    % Graphics of the Variation of Global Stresses Along the Thickness
    
    figure(1);
    plot(Stress_global_x,Th);
    title('Global Stress_x Along the Thickness')
    xlabel('Global Stress_x (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';    
    
    figure(2);
    plot(Stress_global_y,Th);
    title('Global Stress_y Along the Thickness')
    xlabel('Global Stress_y (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(3);
    plot(Stress_global_xy,Th);
    title('Global Stress_x_y Along the Thickness')
    xlabel('Global Stress_x_y (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    % Graphics of the Variation of Global Strains Along the Thickness
    
    figure(4);
    plot(Strain_global_x,Th);
    title('Global Strain_x Along the Thickness')
    xlabel('Global Strain_x (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(5);
    plot(Strain_global_y,Th);
    title('Global Strain_y Along the Thickness')
    xlabel('Global Strain_y (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir='reverse';
    
    figure(6);
    plot(Strain_global_xy,Th);
    title('Global Strain_x_y Along the Thickness')
    xlabel('Global Strain_x_y (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    %Graphics of the Variation of Local Stresses Along the Thickness
    
    figure(7);
    plot(Stress_local_1,Th);
    title('Local Stress_1 Along the Thickness')
    xlabel('Local Stress_1 (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(8);
    plot(Stress_local_2,Th);
    title('Local Stress_2 Along the Thickness')
    xlabel('Global Stress_2 (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(9);
    plot(Stress_local_12,Th);
    title('Local Stress_1_2 Along the Thickness')
    xlabel('Local Stress_1_2 (Pa)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    %Graphics of the Variation of Local Strains Along the Thickness
    
    figure(10);
    plot(Strain_local_1,Th);
    title('Local Strain_1 Along the Thickness')
    xlabel('Global Strain_1 (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(11);
    plot(Strain_local_2,Th);
    title('Local Strain_2 Along the Thickness')
    xlabel('Local Strain_2 (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir= 'reverse';
    
    figure(12);
    plot(Strain_local_12,Th);
    title('Local Strain_1_2 Along the Thickness')
    xlabel('Local Strain_1_2 (m/m)');
    ylabel('z (m)');
    grid on
    yticks (h_coordinate);
    ax=gca;
    ax.XGrid = 'off';
    ax.YDir='reverse';
    