function [Y_csc2] = ...
            CSC2...
                (e,cTe,G,G_diag,v_M,v_M_m1,z_M_M,M,z_M,DR,z_M_M_M,...
                    G_diag0,wY,v_S,z_S,Lcsc2,v_G_diag,wG,r_M_S,...
                    RE_switch,SVE_switch,RE_SVE_switch,dm,gsxi,gsyi,...
                    gszi,K,sx,sy,sz,z_M_Si,lamda,gsx,gsy,gsz,mx,my,mz,...
                    syn_data) 


% All references to equations/figures are from "Constrained Spectral 
%   Conditioning for spatial sound level estimation", Spalt et al., 2016

%% Calculate CSC2 variables which are independent of scanning location_____
G_diag_M = repmat(G_diag,1,M);
ltG = tril(G);
ltG = ltG.*G_diag0;
angG = angle(G);
cG = conj(G);
absG = abs(G);
absG_2 = abs(absG).^2;

%% Calculate inputs needed for random error constraints____________________
if RE_switch == 1  ||  RE_SVE_switch == 1
    TK = 2*K;
    pi2 = pi/2;
    sqrt2 = sqrt(2);
    
    %___Calculate then apply normalized error in autospectral magnitude 
    %   (Eqs. (14,16)_____________
    E_Gxx = 1/sqrt(K);
    G_diag_RE = G_diag*(1+E_Gxx);
else
    TK = NaN;
    pi2 = NaN;
    sqrt2 = NaN;
    G_diag_RE = NaN;
end

%% Loop through grid points and calculate CSC2 estimate at each____________
Y_csc2 = z_S;
parfor s = v_S 
    %___Calculate inputs needed for CSC_________________________________
    [angPSF,angPSF_pi2,exp_iangPSF_pi2,SVE_sin_theta] = ...
        CSC_Inputs...
            (e(:,s),cTe(s,:),r_M_S(:,s),sx(s),sy(s),sz(s),gsxi,gsyi,...
                gszi,z_M_Si,dm,SVE_switch,RE_SVE_switch,v_M_m1,lamda,...
                gsx,gsy,gsz,mx,my,mz,M,v_G_diag,z_M_M,syn_data,v_M);
            
    %___Projected magnitude of cross-spectra onto orthogonal vector___
    cos_d_ang = cos(angPSF_pi2-angG);

    %___Based on constraint selection, constrain Wb_______________________
    if RE_switch == 1
        %___Apply random error constraints to cross-spectral magnitude and
        %   cosine projection terms of Eq. (7)___________________________
        [absG_RE,cos_d_ang_RE] =...
            RE_Constraints_CSC2...
                (absG,G_diag,K,TK,pi2,angPSF_pi2,angG,sqrt2,angPSF,...
                    cos_d_ang,v_G_diag,absG_2,z_M_M,v_M_m1,M);

        %___Form weight vector magnitude via Eq. (7) using random error 
        %   constraints (Eqs. (14-17) & Fig. 2)____
        Wb = absG_RE.*cos_d_ang_RE./G_diag_RE;    
    elseif SVE_switch == 1
        %___Apply steering vector error constraints to cosine projection 
        %   term of Eq. (7)______________________
        [cos_d_ang_SVE] = ...
            SVE_Constraints_CSC2...
                (cos_d_ang,SVE_sin_theta);
        
        %___Form weight vector magnitude via Eq. (7) using steering vector 
        %   error constraints (Eqs. (19-21) & Fig. 4)____
        Wb = absG.*cos_d_ang_SVE./G_diag;
    elseif RE_SVE_switch == 1
        %___Apply random error constraints to cross-spectral magnitude and
        %   cosine projection terms of Eq. (7)___________________________
        [absG_RE,cos_d_ang_RE] =...
            RE_Constraints_CSC2...
                (absG,G_diag,K,TK,pi2,angPSF_pi2,angG,sqrt2,angPSF,...
                    cos_d_ang,v_G_diag,absG_2,z_M_M,v_M_m1,M);

        %___Apply steering vector error constraints to cosine projection 
        %   term of Eq. (7)________________________
        [cos_d_ang_RE_SVE] = ...
            SVE_Constraints_CSC2...
                (cos_d_ang_RE,SVE_sin_theta);
        
        %___Form weight vector magnitude via Eq. (7) using random error 
        %   constraints (Eqs. (14-17) & Fig. 2) & steering vector 
        %   constraints (Eqs. (19-21) & Fig. 4)____
        Wb = absG_RE.*cos_d_ang_RE_SVE./G_diag_RE;   
    else
        %___Form unconstrained CSC weight vector______
        Wb = absG.*cos_d_ang./G_diag_M;
    end
    
    %___CSC uses NaN to remove ref sensors, but CSC2 uses 0___
    Wb(isnan(Wb)==1) = 0;    
    
    %___Remove sensors which violate conditions of Eq. (9)___
    Wb = Wb.*G_diag0;     % Autospectra cannot be refs               
    Wb(abs(Wb)>1) = 0;    % |Wb|>1 cannot be used  

    %___Join projected magnitude onto orthogonal vector to orthogonal 
    %   phase angle___________
    Wb = Wb.*exp_iangPSF_pi2;

    %___Compute complex conjugate____
    cWb = conj(Wb);

    %___Perform Constrained Spectral Conditioning on CSM_______________
    [Y_csc2(s)] = ...
        CSC2_Core...
            (v_M,z_M_M,v_M_m1,M,G,ltG,z_M_M_M,z_M,DR,G_diag,e(:,s),...
                cTe(s,:),Lcsc2,v_G_diag,wY,wG,cG,Wb,cWb);
end


end