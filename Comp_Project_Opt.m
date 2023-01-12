% Composite Materials Project
% For Choose Optimized Materials and Ply Orientations for a Symmetric Stacking Sequence for 6 Plies

clear
clc
clear all
format long

% Properties of Different Materials from Table 1: Material Selection Table
% Boron/Epoxy, AS Carbon/Epoxy, T-300/Epoxy, HMS Carbon/Epoxy, GY-70/Epoxy,
% Kevlar 49/Epoxy, E-Glass/Epoxy

% All Properties for All Materials from Table 1 

E1= [207,128,138,171,262,76,32]*10^9;   % Longitudinal Elastic Modulus [Pa]
E2= [19,9,10,13.8,8.3,5.5,4.8]*10^9;    % Transverse Elastic Modulus [Pa]
G12= [6.4,5.7,6.5,5.9,4.1,2.1,4.8]*10^9;   % Shear Modulus [Ga]
v12= [.21,.25,.21,.20,.25,.34,.30];   % Major Poisson's Ratio 
Ro= [1990,1540,1550,1630,1690,1380,1800];   % Density [kg/m^3]
Str_11= [1585,1448,1448,827,586,1379,1103]*10^6;   % Ultimate Longitudinal Strength [Pa]
Str_22= [63,62,45,86,41,28,97]*10^6;   % Ultimate Transverse Strength [Pa]
Str_12= [131,60,62,72,97,60,83]*10^6;   % Ultimate In-Plane Shear Strength [Pa]



b=0;   % Number of All Combinations
y=0;   % Number of Stacking Sequence Combinations that Provide the Desired Conditions
nth=1; % Number of Compared Combinations
Ex_spe_max= 0;   % Highest Inplane Speciﬁc Longitudinal Modulus [Pa*kg/m^3]
Mass_Min= 1; % Minumum Weight [kg]

% To Create Different Symmetric Stacking Sequences
% With Different Materials and Angles
% 
for i=1:7 % There are 7 Different Materials for 1st and 6th Plies
    for j=1:7 % There are 7 Different Materials for 2nd and 5th Plies
        for k=1:7 % There are 7 Different Materials for 3rd and 4th Plies
            for ii=-90:15:90 % There are 13 Different Angle for 1st and 6th Plies
                for jj=-90:15:90 % There are 13 Different Angle for 2nd and 5th Plies
                    for kk=-90:15:90 % There are 13 Different Angle for 3rd and 4th Plies
                        
                        angle=[ii,jj,kk,kk,jj,ii];
                        E_1=[E1(i),E1(j),E1(k),E1(k),E1(j),E1(i)];
                        E_2=[E2(i),E2(j),E2(k),E2(k),E2(j),E2(i)];
                        v_12=[v12(i),v12(j),v12(k),v12(k),v12(j),v12(i)];
                        G_12=[G12(i),G12(j),G12(k),G12(k),G12(j),G12(i)];
                        Density=[Ro(i),Ro(j),Ro(k),Ro(k),Ro(j),Ro(i)];
                        Strength_Ultimate_11=[Str_11(i),Str_11(j),Str_11(k),Str_11(k),Str_11(j),Str_11(i)];
                        Strength_Ultimate_22=[Str_22(i),Str_22(j),Str_22(k),Str_22(k),Str_22(j),Str_22(i)];
                        Strength_Ultimate_12=[Str_12(i),Str_12(j),Str_12(k),Str_12(k),Str_12(j),Str_12(i)];
                        
                        b=b+1;  % All Combinations
                        
                        
                        % Initial Conditions and Given Parameters
                        
                        a = 0.4;                      % An Edge-Length of The Square Plate [m]
                        NM = [50e3;-50e3;1e3;-2;7;1]; % Load [N_x;N_y;N_xy;M_x;M_y;M_xy] [N/m;(N*m)m]
                        plyno = 6;                    % Ply Number of The Symmetric Laminate
                        t = 0.25e-3;                  % Ply Thickness [m]
                        ns = 2;                       % Factor of Safety
                        Mass = 0;                     % Mass [kg]
                        no_fail=0;                    % Number of No-Fail Plies
                        S = zeros(3,3); % Reducing Compliance Matrix
                        R = zeros(3,3); % Reuter Matrix
                        R(1,1) = 1;
                        R(2,2) = 1;
                        R(3,3) = 2;
                        A = zeros(3);  % A Matrix
                        B = zeros(3);  % B Matrix
                        D = zeros(3);  % D Matrix
                        
                        h = plyno*t;                % Total Thickness of Laminate [m]
                        h_coordinate(1) = -h/2;     % Locations of the Ply Surfaces [m]
                        Volume = plyno*t*a^2;       % Total Volume [m^3]
                        
                        
                        for n = 1:plyno % Considering ply by ply
                            
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
                            
                            e0 = ABBD\NM;    % [Epsilon*0_x;Epsilon*0_x;Gamma*0_xy ; K_x;K_y;K_xy] [m/m; 1/m]
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
                        
                        for n = 1:plyno % Considering ply by ply
                            
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
                        
                        for n=1:plyno % Considering ply by ply
                            
                            %Tsai-Hill Failure Theory is Applied
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
                        
                        
                        % The stresses in individual plies must be lower than the ply strength given in table.
                        % Magnitudes of global mid-plane strains must be lower than 5*10^-3 m/m.
                        % Magnitudes of global mid-plane curvatures must be lower than 4 m^-1.
                        
                        if (no_fail)==6 && (max(abs( Strain_midplane))<5e-3/ns) && (max(abs(Curvature_midplane))<4/ns)
                            
                            y=y+1;    % Number of Stacking Sequence Combinations that Provide the Desired Conditions
                            
                            Inplane_Mod_spe=Specific_Modulus_Inplane;   % Specific In-Plane Modulus of Combinations that Provide the Desired Conditions
                            mass=Mass;  % Mass of Stacking Sequence Combinations
                            
                            % To Compare and Find Highest Inplane Speciﬁc Modulus with Minumum Weight
                            
                            if  Inplane_Mod_spe(1)>= Ex_spe_max && mass<= Mass_Min
                                
                                nth=nth+1;   % For Keeping the Properties of Compared Orientations
                                
                                Ex_spe_max= Inplane_Mod_spe(1);  % Maximum Specific Young's Modulus [Pa]
                                Mass_Min= mass ;  % Minumum Mass [kg]
                                
                                % For Creating a Struct Keeping Essential Properties and Values of Oriantations
                                Combination(nth).angle=[ii,jj,kk]; % Angle Orientations of Laminates
                                Combination(nth).material= [i,j,k];  % Material Orientation of Laminates
                                Combination(nth).Mass= mass;  % Mass of Laminates
                                Combination(nth).Sp_Modulus_Inplane= Inplane_Mod_spe;  % Specific Inplane Modulus of Laminates
                                Combination(nth).Strain_Mid= e0;  % Mid-plain Strain and Curvatures of Laminates
                                Combination(nth).localStress= Stress_local;  % Local Stresses of Laminates
                                Combination(nth).StiffnessMatrix= ABBD;  % Stiffness Matrix of Laminates 
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
end