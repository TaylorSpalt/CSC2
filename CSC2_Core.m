function [Y_csc2] = ...
            CSC2_Core...
                (v_M,z_M_M,v_M_m1,M,G,ltG,z_M_M_M,z_M,DR,G_diag,es,cTes,...
                    Lcsc2,v_G_diag,wY,wG,cG,Wb,cWb)

  
% All references to equations are from "Efficiency improvements to 
%   Constrained Spectral Conditioning for spatial source level estimation", 
%   T.Spalt, 2018

%% Loop through reference sensors and form CSC2 versions of the CSM________
G_CSC2 = z_M_M_M;
Loss = z_M;
for mr = v_M
    %___Define here to speed processing within loops below___
    G_mr_mr = G(mr,mr);
    
    %___Form cross-spectral orthogonal weight vector (bracketed term in 
    %   Eq. (7))______________________________________________________
    Wb_Gx = z_M_M;
    for m = v_M_m1
        G_mr_m = G(mr,m);
        Wb_mr_m = Wb(mr,m);
        G_mr_mr_Wb_mr_m = G_mr_mr*Wb_mr_m;
        for mp = m+1:M                
            Wb_Gx(mp,m) = G_mr_m*cWb(mr,mp) + cG(mr,mp)*Wb_mr_m - ...
                            G_mr_mr_Wb_mr_m*cWb(mr,mp);
        end
    end

    %___Subtract orthogonal weight vector from cross-spectral entries
    %   and fill in CSM's upper triangular (Eq. (7))_____
    clear Gx;
    Gx = ltG - Wb_Gx;        
    Gx = Gx + ctranspose(Gx);
    
    %___Subtract orthogonal weight vector from autospectral entries if
    %   autospectra are being used_______________________________________
    Ga = z_M_M;     % Initialize CSM to 0s
    if DR == 0
        %___Compute orthogonal signal to autospectral entries (bracketed 
        %   term in Eq. (6))____________________________________________
        Wb_Ga = z_M;
        for m = v_M 
            Wb_Ga(m) = G(mr,m)*cWb(mr,m) + conj(G(mr,m)*cWb(mr,m)) - ...
                            G_mr_mr*Wb(mr,m)*cWb(mr,m); 
        end
        Wb_Ga = real(Wb_Ga);
        
        %___Fill in CSM diagonal with processed autospectra (Eq. (6))___
        Ga(v_G_diag) = G_diag - Wb_Ga;

        %___Calculate ratio of selected reference sensor magnitude to 
        %   the mean magnitude of all other sensors' autospectra 
        %   before processing__________________________________________
        diag_G1 = G_diag;
        diag_G1(mr) = NaN;   % Same sensor cannot be used as ref
        sref_mag = G_diag(mr)/nanmean(diag_G1);

        %___Use this ratio to scale the output value for the  
        %   the ref sensor since its autospectrum cannot be  
        %   processed using itself (Eq. (8))_____________
        Wb_Ga(mr) = NaN;
        Ga(mr,mr) = G_diag(mr) - sref_mag*nanmean(Wb_Ga);
    end
    
    %___Merge CSMs____________
    G_CSC2(:,:,mr) = Gx + Ga;

    %___Loss function calculation_________________________
    if strcmp(Lcsc2,'L1') == 1
        %___Calculate L1 norm________________________
        Loss(mr) = sum(sum(abs(wG.*G_CSC2(:,:,mr))));
    elseif strcmp(Lcsc2,'L2') == 1
        %___Calculate L2 norm___________________________
        Loss(mr) = sum(sum(abs(wG.*G_CSC2(:,:,mr)).^2));
    elseif strcmp(Lcsc2,'FDBF') == 1
        %___Calculate FDBF______________________
        Loss(mr) = cTes*(wG.*G_CSC2(:,:,mr))*es;
    else
        disp('No loss function specified');
        pause;
    end
end

%% Optimum CSM based on reference sensor which minimizes the chosen norm___
[~,opt_refs] = min(real(Loss));

%% Compute FDBF with optimal CSM___________________________________________
if strcmp(Lcsc2,'FDBF') == 1
    Y_csc2 = max(0,real(Loss(opt_refs)))*wY;
else
    Y_csc2 = max(0,real(cTes*(wG.*G_CSC2(:,:,opt_refs))*es))*wY; 
end


end