%% Main code for STZ
% Kesson
% set constants
global H0 volume0 mesh voxel0 S X0 EpsAvg 
global hole
global Mu Nu
global Sigma 
global Eps Eps0 HTemp
global muo nuo lambdao alphao betao Epso
num=0;
while 1
    if num==0
%     dens = 6125;
%     hcap = 420;
%     heat_cont1 = dens*hcap;
%     heat_cont2 = 3e12;
%     del_t = 100;    % interval of integral in thermal diffusion %
%     dcht  = 1e-13;   % characteristic time for thermal diffusion %
%     divdt = dcht*1; % buffer for thermal diffusion time and kMC time %
%     dtt   = 0;
    mesh=[128 128];
    
    inputdata; 
    % for 2D %
    nu = nu*2;
    % 2D isotropic medium - not 3D %
    mu = E/(2+2*nu);
    lambda = E*nu/(1-nu^2);
    alpha = (1+nu)/2;
    beta  = nu;

    VOXEL_HEIGHT_IN_NM = 1.7;
    H0 = [1 0;0 1] * mesh(2) * VOXEL_HEIGHT_IN_NM;
    G0 = 2*pi*inv(H0)';
    
    volume0 = abs(det(H0));
    voxel0 = volume0;
    for i = 1 : length(mesh)
        s{i}   = 0 : 1/mesh(i) : 1-1/mesh(i)/2;
        voxel0 = voxel0 / mesh(i);
    end; 

%     hole = [];
%     Nu = ones(mesh(2),mesh(1)) * nu;
%     Mu = ones(mesh(2),mesh(1)) * mu;
    
    Eps0{1} = X0{1}*0;
    Eps0{2} = X0{1}*0;
    Eps0{3} = X0{1}*0;
    HTemp   = X0{1}*0;
%     EPSo = Eps0;
%     Eps  = Eps0;
%     Epso = Eps0;
    
    units;
    T = 293;
    HTemp   = Eps0{1} * 0 + T - 0.001*randn(mesh(2),mesh(1));
    % % kT = T * BOLZ * J_IN_EV;
    kT = HTemp * BOLZ * J_IN_EV;
    % trialfreq = voxel0 * 1e13;
    % ============================================
    trialfreq = 1E13;
    % ============================================
    Q0amp = 5;               % eV
    % % Q0amp = 3.9;               % eV
    Q0ampSur = 3;            % eV for surface or interface
    OmegaF = voxel0 / 10;    % Magnitude of transformation volume

    % strainrate = 10^2; %
    % strainrate = 1; %
    % strainrate = 1e-2; %
    strainrate = -1e-3;
    % strainrate = 1e-6; %
    % strainrate = 1e-8; %
    % 1/s %
    
    
    step = 0; time = 0; printperiod = 5; npic = 0; more off;
    strainincrement = -1e-4; hlast = 0; Amp_DE = 1; TempCap = 600;
    % assign the index to the voxel on the surface
    Modes = 20;
    % How much elevation, tilt, and roughness of energy landscape in strain space %
    RandomRatio = 0.2;
    Q0 = Q0amp * (1-RandomRatio + RandomRatio*rand(mesh(2),mesh(1),Modes)); 
    % pure shear transformations %
    DEps0{1} = Amp_DE * randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
    DEps0{2} = -DEps0{1};
    DEps0{3} = Amp_DE * randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
    % diffusional relaxation time: Acta Materialia 57 (2009) 881?92:Eq(16)in PYZ %
    tTemp = 1 ./ (trialfreq * exp(-0.37 ./ kT));
    % % tTemp = 1 / (trialfreq * exp(-0.37 * 2.93/ kT));
    % time of last transformation %
    tlast = Eps0{1} * 0 - 100*tTemp;
    % Deborah number distribution %
    DebDi = Eps0{1} * 0;
    % rmdir('Jpgs','s');mkdir('Jpgs','Eps');
    strainhistory = [0]; stresshistory = [0]; 
    avgtemhistory = [0]; mintemhistory = [0];
    CTtemphistory = [0]; Cstreshistory = [0];
    avgEps11history = [0]; avgEps22history = [0];
    avgEps12history = [0];
    sumEps11history = [0]; sumEps22history = [0];
    sumEps12history = [0];
    % supercell strain %
    % EpsAvg = [0 0 0];
    EpsAvg = 0.001*[1 -nu 0];
    % residual stress generator: satisfies self equilibration %
    ResidualStrainSTD = OmegaF/voxel0 * 0.1;
    Eps0{1} = randn(mesh(2),mesh(1)) * ResidualStrainSTD / 2; 
    Eps0{2} = -Eps0{1};
    Eps0{3} = randn(mesh(2),mesh(1)) * ResidualStrainSTD / 2; 
    Eps0{1} = Eps0{1} - mean(mean(Eps0{1}));
    Eps0{2} = Eps0{2} - mean(mean(Eps0{2}));
    Eps0{3} = Eps0{3} - mean(mean(Eps0{3}));
    %%% call for Abaqus to calculate stress and strain redistribution 
    Eps0{1} = X0{1}*0;
    Eps0{2} = X0{1}*0;
    Eps0{3} = X0{1}*0;
    dataoutput; %%%%output residual
    datinput;   %%%%%%get result
    %%%
    ResidualSigma = Sigma;
    % Damage trackers ps:PermSoft and TempSoft might be kp,kt in table 1%
    % PermSofteningCap might be Eq(15) %
    PermSoft = 10;
    TempSoft = 30;
    PermSoftening = X0{1}*0;
    PermSofteningCap = -log(0.8);
    TempSoftening = X0{1}*0;
    % incremental strain in reference to the present %
    Eps0{1} = X0{1}*0;
    Eps0{2} = X0{1}*0;
    Eps0{3} = X0{1}*0;
    num=num+1;
    else 
    %% KMC
    datinput;
        while 1
     % softening - mimicking free volume creation and annihilation %
     % Softening = PermSoftening + TempSoftening .* exp((tlast-time)/tTemp);
      Softening = PermSoftening + TempSoftening .* exp((tlast-time)./tTemp);
      SIGMA{1} = Sigma{1} + ResidualSigma{1};
      SIGMA{2} = Sigma{2} + ResidualSigma{2};
      SIGMA{3} = Sigma{3} + ResidualSigma{3};
      for m = 1 : Modes,
        % bias due to stress distribution %
        DWork = ( SIGMA{1} .* DEps0{1}(:,:,m) + ...
                  SIGMA{2} .* DEps0{2}(:,:,m) + ...
              2 * SIGMA{3} .* DEps0{3}(:,:,m) ) * voxel0 * 1e-18 * J_IN_EV;
        Q(:,:,m) = Q0(:,:,m) .* exp(-Softening) - DWork/2;
        Qq(:,:,m) = Q(:,:,m) ./ kT;
      end
    % %   DebDi = (time-tlast)./tTemp;
      DebDi = tTemp./(time-tlast);

      if ( min(min(min(Q))) < 0 )
        % athermal plasticity %
        [Qmin,m] = min(Q,[],3);
        [Row,Col] = find(Qmin<0);
    %     break %
        for p = 1 : length(Row)
          row = Row(p); col = Col(p); mode = m(row,col);
          VoxelShearTransform;
        end
        %  call for Abaqus %
        dataoutput;
        datinput;      
        %%%%%%%%%%%%%%%%%%%%%%%
%         HomogeneousStressSolver;
    %     InhomoStressSolver_Ju;
%         EpsAvg = EpsAvg - IsotropicInv([0 Stress(2) Stress(3)], lambda, mu);
      else
    % %     q = reshape(Q, mesh(2)*mesh(1)*Modes, 1);
    % %     rate = trialfreq * exp( - q / kT );
        q = reshape(Qq, mesh(2)*mesh(1)*Modes, 1);
        rate = trialfreq * exp( - q );    
        timeincrement = 1 / sum(rate);
%         PreHeatConduction;
        if ( timeincrement < strainincrement/strainrate )
          % thermal plasticity %
          dice = rand(1);
          roll = cumsum(rate) / sum(rate);
          idx = min(find( (roll > dice)));
          mode = ceil( idx / mesh(2) / mesh(1) );
          col = ceil( idx / mesh(2) - (mode - 1)*mesh(1) );
          row = idx - (mode-1)*mesh(2)*mesh(1) - (col-1) * mesh(2);
          VoxelShearTransform;
          time = time + timeincrement;
%           EpsAvg = EpsAvg + strainrate*timeincrement*[1 -nu 0];   
          EpsAvg = strainrate*timeincrement*[1 -nu 0]; 
          %%%%%%  call for Abaqus, need to modify .inp file
          revise
          %  call for Abaqus %
          dataoutput;
          datinput;      
        %%%%%%%%%%%%%%%%%%%%%%%
          Stress(1) = mean(mean(Sigma{1}));        
%           HomogeneousStressSolver;
    %       InhomoStressSolver_Ju;
%           EpsAvg = EpsAvg - IsotropicInv([0 Stress(2) Stress(3)], lambda, mu);


%           Stress(1) = Stress(1) + E * strainrate*timeincrement;
%           Sigma{1}  = Sigma{1}  + E * strainrate*timeincrement;
    %       break %
        else
          % pure elasticity %
          time = time + strainincrement/strainrate;
%           EpsAvg = EpsAvg + strainincrement*[1 -nu 0];
          EpsAvg = strainincrement*[1 -nu 0];
          %%%%%%  call for Abaqus, need to modify .inp file
          revise
          %%%%%%%%%%%%%%%%%%%%%%%%%
          Stress(1) = mean(mean(Sigma{1}));
%           Stress(1) = Stress(1) + E * strainincrement;
%           Sigma{1}  = Sigma{1}  + E * strainincrement;
        end
      end

        % analysis of accumulated plastic strain %
%         Accumulatedstrain;
      % heat distribution %
%         HeatConduction;

    % %     % heat-sink for kMC model %
    % %        d_HT  = minHTemp - T;
    % %        HTemp = HTemp - d_HT;
    % %     % heat-sink for large model %
% % % %            CTemp(:,1)        = T;
% % % %            CTemp(:,cmesh(1)) = T;
% % % %            CTemp(1,:)        = T;
% % % %            CTemp(cmesh(2),:) = T;      

           % cap temperature up bound %
% % % %            if ( max(max(HTemp)) > TempCap)
% % % %               [TRow,TCol] = find(HTemp > TempCap);
% % % %               for p = 1 : length(TRow)
% % % %                   trow = TRow(p); tcol = TCol(p);
% % % %                   HTemp(trow,tcol) = TempCap;
% % % %               end  
% % % %            end
           % %
%            kT    = HTemp * BOLZ * J_IN_EV;
           tTemp = 1 ./ (trialfreq * exp(-0.37 ./ kT));

      if ( EpsAvg(1)-strainhistory(length(strainhistory)) <= ...
           printperiod * strainincrement * 0.99 )
           AvgHTemp = mean(mean(HTemp));
           minHTemp = min(min(HTemp));
        strainhistory = [strainhistory EpsAvg(1)];
        stresshistory = [stresshistory Stress(1)];
        avgtemhistory = [avgtemhistory AvgHTemp]; 
    % %     maxtemhistory = [maxtemhistory ];
        mintemhistory = [mintemhistory minHTemp];  
        CTtemphistory = [CTtemphistory CTemp(crow,ccol)];   
        Cstreshistory = [Cstreshistory CStress(1)];
        avgEps11history = [avgEps11history CEps0_1]; 
        avgEps22history = [avgEps22history CEps0_2];
        avgEps12history = [avgEps12history CEps0_3];
        sumEps11history = [sumEps11history SEps0_1]; 
        sumEps22history = [sumEps22history SEps0_2];
        sumEps12history = [sumEps12history SEps0_3];    
    %     fprintf ('AvgHTemp = %g, minHTemp = %g\n',AvgHTemp,minHTemp);
% % % %         H = H0 * (eye(2) + T2VoightToFull(EpsAvg));
% % % %         X{1} = S{1} * H(1,1) + S{2} * H(2,1);
% % % %         X{2} = S{1} * H(1,2) + S{2} * H(2,2);
% % % %         for i = 1 : 2,
% % % %           U{i} = ( F{i} - alphao * G .* K{i} ) / muo / complex(0,1) ./ K{3};
% % % %     %       u{i} = real(ifft2(U{i})) / voxel0;      
% % % %           u{i} = real(ifft2(EPS{i})) / voxel0;
% % % %           X{i} = X{i} + u{i};
% % % %         end

        % position for large space model %
% % % % %         CH = CH0 * (eye(2) + T2VoightToFull(CEpsAvg));
% % % % %         CX{1} = CS{1} * CH(1,1) + CS{2} * CH(2,1);
% % % % %         CX{2} = CS{1} * CH(1,2) + CS{2} * CH(2,2);
% % % % %         for i = 1 : 2,
% % % % %           CU{i} = ( CF{i} - alphao * CG .* CK{i} ) / muo / complex(0,1) ./ CK{3};
% % % % %     %       u{i} = real(ifft2(U{i})) / voxel0;      
% % % % %           Cu{i} = real(ifft2(CEPS{i})) / Cvoxel0;
% % % % %           CX{i} = CX{i} + Cu{i};
% % % % %         end    

        npic = npic + 1;

      end;
      step = step + 1;  

      if (abs(EpsAvg(1)) > 0.1) 
    %      unix('./vidcompress');
        break 
      end
        end
		num=num+1;
    end  
end
