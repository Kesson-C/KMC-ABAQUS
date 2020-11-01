Eps0{1}(row,col) = Eps0{1}(row,col) + DEps0{1}(row,col,mode);
Eps0{2}(row,col) = Eps0{2}(row,col) + DEps0{2}(row,col,mode);
Eps0{3}(row,col) = Eps0{3}(row,col) + DEps0{3}(row,col,mode);

DeltaSoftening = DEps0{1}(row,col,mode)^2 + ...
    DEps0{2}(row,col,mode)^2 + 2*DEps0{3}(row,col,mode)^2;

PermSoftening(row,col) = PermSoftening(row,col) + DeltaSoftening * PermSoft;
if (PermSoftening(row,col) > PermSofteningCap) 
    PermSoftening(row,col) = PermSofteningCap;
end

TempSoftening(row,col) = DeltaSoftening * TempSoft;
tlast(row,col) = time;
  
Q0(row,col,:) = Q0amp * (1-RandomRatio + RandomRatio*rand(Modes,1));
% initial temperature distribution, stress GPa = 1e9 Pa, Pa = J/m^3 %
del_work = ( SIGMA{1}(row,col) * DEps0{1}(row,col,mode) + ...
             SIGMA{2}(row,col) * DEps0{2}(row,col,mode) + ...
         2 * SIGMA{3}(row,col) * DEps0{3}(row,col,mode) ) * 1e9;
del_T    = del_work/heat_cont1;
% del_T    = 0.1 * 0.9 * del_T;
del_T    = 0.0;
% % HTemp(row,col)   = HTemp(row,col) + del_T;
% % del_T            = del_T/(mesh(1)*mesh(2));
CTemp(crow,ccol) = CTemp(crow,ccol) + 0.1*del_T;

% cap temperature up bound %
  if ( CTemp(crow,ccol) > TempCap)
       CTemp(crow,ccol) = TempCap;
  end

% heat line source %
for ii = 1+3:cmesh(2)-3
    CTemp(ii,ccol-crow+ii) = CTemp(crow,ccol);
end
  
% ============================================
% Modified by PY Zhao
DEps0{1}(row,col,:) = Amp_DE * randn(Modes,1)/2*OmegaF/voxel0;
DEps0{3}(row,col,:) = Amp_DE * randn(Modes,1)/2*OmegaF/voxel0;
DEps0{2}(row,col,:) = -DEps0{1}(row,col,:);
% ============================================

