function [absG_RE,cos_d_ang_RE] =...
            RE_Constraints_CSC2...
                (absG,G_diag,K,TK,pi2,angPSF_pi2,angG,sqrt2,angPSF,...
                    cos_d_ang,v_G_diag,absG_2,z_M_M,v_M_m1,M)
            

%------------------------------------------------------------------------  
% Refer to Section 3.1 of "Constrained Spectral Conditioning for spatial
%   sound level estimation", Spalt et al., 2016 for details. 
%------------------------------------------------------------------------

%% Calculate the error in the cross-spectral magnitude and use it to reduce
%   its component in the weight vector_____________________________________
%___Calculate ordinary coherence of the CSM______________
Cxy = z_M_M;
for m = v_M_m1
    for mp = m+1:M   
        Cxy(m,mp) = absG_2(m,mp)/(G_diag(m)*G_diag(mp));
    end
end
Cxy = Cxy + transpose(Cxy);
Cxy(v_G_diag) = NaN;           

%___Calculate normalized error in coherence (Eqs. (14,15))____
E_Cxy = sqrt2*(1-Cxy)./(sqrt(Cxy*K));
Cxy(E_Cxy>=1) = NaN;	

%___Calculate normalized error in cross-spectral mag (Eqs. (14,15))___
E_absG = 1./(sqrt(Cxy*K));
E_absG(E_absG>=1) = NaN;

%___Apply cross-spectral magnitude error (Eq. (16))___
absG_RE = absG.*(1-E_absG);    
absG_RE(absG_RE<0) = NaN;

%% Calculate the reduced cosine projection onto the orthogonal phase vector
%   due to the standard deviation in the phase angle_______________________
%___Calculate the std dev in phase angle. Use 3 std devs in order to
%   encompass full range of angle deviations. (Eq. (14))____
sigma = 3*abs(asin(sqrt(1-Cxy)./sqrt(TK*Cxy)));     
sigma(sigma>=pi2) = NaN;  

%___Determine if any of the PSF angles lie within +/-sigma of the cross-
%   spectral phase of the sensor being processed (Eq. (15) & Fig. 2)______
angG_p_sigma = angG + sigma;
angG_m_sigma = angG - sigma;
nan_refs1 = z_M_M;
nan_refs1(angG_p_sigma>angPSF) = 1;
nan_refs1(angG_m_sigma<angPSF) = nan_refs1(angG_m_sigma<angPSF)+1;
nan_refs1(nan_refs1<2) = 0;
nan_refs1(nan_refs1>0) = 1;
nan_refs2 = z_M_M;
nan_refs2(angG_m_sigma<pi) = 1;
nan_refs2(angPSF<angG_p_sigma) = nan_refs2(angPSF<angG_p_sigma)+1;
nan_refs2(angPSF>(angG_m_sigma+2*pi)) = ...
                            nan_refs2(angPSF>(angG_m_sigma+2*pi))+1;
nan_refs2(nan_refs2<2) = 0;
nan_refs2(nan_refs2>0) = 1;
nan_refs3 = z_M_M;
nan_refs3(angG_p_sigma>pi) = 1;
nan_refs3(angPSF>angG_m_sigma) = nan_refs3(angPSF>angG_m_sigma)+1;
nan_refs3(angPSF<(angG_p_sigma-2*pi)) = ...
                            nan_refs3(angPSF<(angG_p_sigma-2*pi))+1;
nan_refs3(nan_refs3<2) = 0;
nan_refs3(nan_refs3>0) = 1;
nan_refs = nan_refs1 + nan_refs2 + nan_refs3;

%___Compute the cross-spectra cosine projections onto the orthogonal phase 
%   vector, and save their signs for later use__________
cos_d_ang_p_sigma = cos(angPSF_pi2 - angG_p_sigma);   
cos_d_ang_m_sigma = cos(angPSF_pi2 - angG_m_sigma);
s_cos_d_ang = sign(cos_d_ang);          
s_cos_d_ang_p_sigma = sign(cos_d_ang_p_sigma);    
s_cos_d_ang_m_sigma = sign(cos_d_ang_m_sigma);   

%___Do not use altered projections whose sign is not equal to the nominal
%   projection's sign______________________________________
cos_d_ang_p_sigma(s_cos_d_ang_p_sigma~=s_cos_d_ang) = NaN;
cos_d_ang_m_sigma(s_cos_d_ang_m_sigma~=s_cos_d_ang) = NaN;

%___Use original cosine projection's sign with minimum projected magnitude
%   of the +/- sigma phase deviations (Eq. (17))__________________________
cos_d_ang_RE = s_cos_d_ang.*...
                    nanmin(abs(cos_d_ang_p_sigma),abs(cos_d_ang_m_sigma));

%___Set any ref sensors to NaN whose PSF phase falls within the +/-sigma
%   region around the cross-spectral phases of the sensor being
%   processed (Eq. (15) & Fig. 2)_______
cos_d_ang_RE(nan_refs==1) = NaN;
cos_d_ang_RE(isnan(sigma)==1) = NaN;
 

end