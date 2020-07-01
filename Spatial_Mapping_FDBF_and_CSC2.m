function Spatial_Mapping_FDBF_and_CSC2


%% User-defined inputs_____________________________________________________
%___Sensor array inputs (E.Sarradj, “Optimal planar microphone 
%   array arrangements,” Fortschritte der Akustik-DAGA, 220-223, 2015)__
M = 31;                 % # of sensors in array
rmax = 25;              % Max radius of any sensor, in
V = 5;                  % Sensor layout spiral parameter (default V=5,
                        %   Vogel's spiral)
H = 0;                  % Sensor layout area parameter (default H=0)

%___Processing inputs_____________________________________________________
Fs = 200000;            % Sampling frequency, Hz
nFFT = 8192;            % FFT block size, samples
K = 1000;               % # of FFT blocks. Minimum=2 due to RE constraint.
f = 4000;               % Narrowband frequency of data, Hz
SNR = 18;               % Signal-to-Noise Ratio, dB (incoherent 
                        %   between sensors)
spert = 0.1;        	% Max distance to randomly perturb user-defined
                        %	source locations, in. (Note: must be <
                        %	min[gsx,gsy,gsz]/2 in order for the perturbed
                        %   location to still pertain to the nominal grid
                        %   point).
pref = 1/((2e-5)^2);	% Reference pressure, pascals
c = 13397.2441;      	% Speed of sound in propagation medium, in/s
syn_data = 1;        	% 1 if data is computer-generated
DR = 0;             	% 1 to set diagonal of CSM to 0
Lcsc2 = 'FDBF';         % Loss function minimization used to choose
                        %   optimal reference sensor for CSC2 processing. 
                        %   'L1' = l1-norm, 'L2' = l2-norm,
                        %   'FDBF' = frequency-domain beamform
                        
%___Cartesian coordinates of scanning grid_______________
sx0 = 0;  	% Center of grid in x-coordinate, in
sy0 = 0;  	% Center of grid in x-coordinate, in
sz0 = 50;  	% Scan grid plane distance from array, in
Sx = 25;  	% # of grid points in x-coordinate        
Sy = 25;  	% # of grid points in y-coordinate    
Sz = 1;    	% # of grid points in z-coordinate
gsx = 1;   	% Spacing between grid points in x-coord, in
gsy = 1;   	% Spacing between grid points in y-coord, in
gsz = 1;   	% Spacing between grid points in z-coord, in

%___CSC constraint inputs__________________________________________
RE_switch = 0;      % 1 to use random error constraints
SVE_switch = 0;     % 1 to use steering vector error constraints
RE_SVE_switch = 1;  % 1 to use both constraints
s_intra = 10;       % Number of grid points within existing points
dm = 0.2;           % Uncertainty of sensor coordinates, in

%% Define sensor array coordinates (E.Sarradj, “Optimal planar microphone 
%   array arrangements,” Fortschritte der Akustik-DAGA, 220-223, 2015)_____
[mx,my] = ...
    Optimal_Array_Sensor_Coordinates...
        (H,V,rmax,M);
mz = zeros(M,1); 

%% Build scanning grid_____________________________________________________ 
[S,sx,sy,sz] = ...
    Build_Scan_Grid...
        (sx0,sy0,sz0,Sx,Sy,Sz,gsx,gsy,gsz);
    
%___Form intragrid point volumetric mesh____
if SVE_switch == 1  ||  RE_SVE_switch == 1
    [gsxi,gsyi,gszi,Si] = ...
        Intragrid_Inputs...
            (s_intra,gsx,gsy,gsz);
else
    gsxi = NaN;
    gsyi = NaN;
    gszi = NaN;
    Si = 0;
end
  
%% Pre-allocation variables________________________________________________
M_p1 = M + 1;
M_2 = M^2;
v_G_diag = 1:M_p1:M_2;
v_M = 1:M;
v_M_m1 = 1:M-1;
z_S = zeros(S,1); 
z_M_S = zeros(M,S);
z_M_Si = zeros(M,Si);
z_M_M = zeros(M);
o_M_M = ones(M);
z_M = zeros(M,1);
z_M_M_M = zeros(M,M,M);
v_S = 1:S; 
v_K = 1:K;
G_diag0 = o_M_M;
G_diag0(v_G_diag) = 0;
two_invK = 2/K;
df = Fs/nFFT;               

%% Calculate array CSM weighting and FDBF normalization constant___________
%___Default is unity weighting matrix. Other user-defined CSM weighting
%   matrix can be specified here.____
wG = o_M_M;     

%___Set CSM diagonal to 0 if selected____
if DR == 1
    wG = wG.*G_diag0;
end

%___Compute normalization constant for FDBF____
wY = 1/sum(wG(:));    

%% Compute frequency-related variables_____________________________________
%___Corresponding freq bin to specified narrowband freq_____
fbin = round(nFFT*f/Fs)+1;	

%___True frequency for given frequency bin___
ftrue = (fbin-1)*df;      

%___Complex angular wavenumber at true freq_____
ik = 1i*2*pi*ftrue/c;  

%___Wavelength at ftrue____
lamda = c/ftrue;

%___Calculate steering vectors (Formulation II from E. Sarradj, 
%   "Three-dimensional acoustic source mapping with different beamforming 
%   steering vector formulations", Adv. Acoust. Vib. (2012))_____________
[e,cTe,r_M_S] = ...
    Compute_Steering_Vectors...
        (sx,sy,sz,mx,my,mz,ik,z_M_S,v_M,syn_data);
    
%% Calculate CSM and true source distribution (if synthetic data being 
%   processed)_____________________________________________________________
if syn_data == 1
    [G,G_diag,Y_true] = ...  
        CSM_and_True_Source...
            (nFFT,c,S,sx,sy,sz,mx,my,mz,z_M_M,K,two_invK,M,v_G_diag,z_S,...
                SNR,spert,sz0,v_S,v_M,v_K,ftrue,e,cTe,wY,wG);
else
    %___If experimental array data is being processed, the user must 
    %   provide G and G_diag___________________________________
    Y_true = NaN;   % True source distribution assumed unknown
end
       
%% Perform spatial filtering via frequency-domain beamforming______________
[Y_fdbf] = ...
    Frequency_Domain_Beamform...
        (e,cTe,G,wY,v_S,z_S,wG);
        
%% Perform spatial filtering via CSC2______________________________________                
[Y_csc2] = ...
    CSC2...
        (e,cTe,G,G_diag,v_M,v_M_m1,z_M_M,M,z_M,DR,z_M_M_M,G_diag0,wY,...
            v_S,z_S,Lcsc2,v_G_diag,wG,r_M_S,RE_switch,SVE_switch,...
            RE_SVE_switch,dm,gsxi,gsyi,gszi,K,sx,sy,sz,z_M_Si,lamda,gsx,...
            gsy,gsz,mx,my,mz,syn_data);
    
%% Post-process and plot results___________________________________________
Plot_Results...
    (Y_true,Y_fdbf,Y_csc2,Sx,Sy,Sz,pref,sx,sy,RE_switch,SVE_switch,...
        RE_SVE_switch);
  
    
end