%{
This script takes a Sellmeier equations of different materials, extracts refractive indices at
given wavelengths, and calculates phasematching offset and group delay
over crystal length. The poling period is calculated for the central pump
wavelength, whereas the actual phasematched pump wavelengths will be offset
to either side.
%}
c = physconst('LightSpeed');

% Set crystal material and polarization axes in crystallographic coordinate
% Example: Paderborn LiNbO3 is z-cut (wafer cut axis) => x is propagation direction, y and z
% are possible
crystal = "KDP";
ax1 = 1;
ax2 = 2;

lambda_s = 810; %nm
lambda_p = 810; %nm, central wavelength meaning center between the two actual phasematched wavelengths
L = 15e-3; %m
temp = 30; % degrees centigrade

lambda_i = 1 / ((1/lambda_s)+(1/lambda_p)); %Again this is just the "average" idler wavelength, centered between the two actual ones.

% get refractive indices:
n_H = Sellmeier(lambda_s,temp,ax1,crystal);
n_V = Sellmeier(lambda_s,temp,ax2,crystal);

% get temporal walkoff for the two differently polarized signals:
temporal_sep = transit(lambda_s,n_H,L) - transit(lambda_s,n_V,L);
disp('Temporal walkoff [ps]:')
disp(temporal_sep*1e12)

% Calculate spacing of fringes:
spacing = abs(2*pi/temporal_sep); %note the temporal seperation in seconds, hence the frequency spacing is in Hertz
disp('Fringe spacing [GHz]:')
disp(spacing*1e-9)
disp('Fringe spacing @400nm [nm]:')
freq = c/400;
freq2 = freq + spacing*1e-9;
disp(c/freq-c/freq2)

% calculate poling period foir QPM. Dirty trick: Average over signal
% refractive indices
n_avg = (n_H + n_V ) / 2;
qpm = - n_avg/lambda_s - (Sellmeier(lambda_p,temp,ax1,crystal)+Sellmeier(lambda_p,temp,ax2,crystal))/(2*lambda_p) + (Sellmeier(lambda_i,temp,ax1,crystal)+Sellmeier(lambda_i,temp,ax2,crystal))/(2*lambda_i);
%qpm = - n_V/lambda_s - (Sellmeier(lambda_p,temp,1,crystal))/(lambda_p) + (Sellmeier(lambda_i,temp,1,crystal))/(lambda_i);
%qpm = - n_avg/lambda_s - Sellmeier(lambda_p,temp,1,crystal)/lambda_p + Sellmeier(lambda_i,temp,1,crystal)/lambda_i;
pp = 1/qpm; % in nm
disp('Poling period [um]:')
disp(pp*1e-3)

% Actual phasematched pump wavelengths is next.
% h-polarized:
syms l_p
eqn = phasemismatch(lambda_s,l_p,qpm,ax1,temp,crystal) == 0;
result = vpasolve(eqn, l_p);
result = vpa(result);
lambda_p_H = result(find(result>300 & result<1000));
disp('Phasematched H-pump wavelength [nm]:')
disp(lambda_p_H)
% v-polarized:
syms l_p
eqn = phasemismatch(lambda_s,l_p,qpm,ax2,temp,crystal) == 0;
result = vpasolve(eqn, l_p);
result = vpa(result);
lambda_p_V = result(find(result>300 & result<1000));
disp('Phasematched V-pump wavelength [nm]:')
disp(lambda_p_V)
disp("Wavelength shear:")
disp(abs(lambda_p_H-lambda_p_V))

% Helper functions:
function T = transit(wl,n,l)
    T = l / (physconst('LightSpeed')*n);
end

function d_k = phasemismatch(lambda_s,lambda_p,qpm,axis_s,temp,crystal)
    lambda_i = 1 / ((1/lambda_s)+(1/lambda_p));
    d_k = Sellmeier(lambda_s,temp,axis_s,crystal)/lambda_s + Sellmeier(lambda_p,temp,axis_s,crystal)/lambda_p - Sellmeier(lambda_i,temp,axis_s,crystal)/lambda_i + qpm;
    %d_k = d_k*2*pi;
end

% Here follow all the Sellmeier equations
function N = Sellmeier(wl, T, axis, crystal)
    %{
    Provides bulk Sellmeier equations for stuff.
    Congruently grown LiNbO3: The Sellmeier equations are from Edwards & Lawrence (ordinary)
    and Jundt (extraordinary)
    
    The equations are provided as functions with arguments
    being the wavelength given in [um] and the temperature given in [C]
    returning the squared refractive index.

    For Z-CUT LN (Propagation in x-direction)
    axes definition: x=0, y=2, z=2
    Crystal frame of reference! In the lab: x->z

    BBO: Bulk BBO Sellmeier, supposedly at room temperature (?), but not mentioned in the paper 
    taken from Dongxiang Zhang, Yufei Kong, Jing-yuan Zhang,
    Optics Communications; Volume 184, Issues 5–6, 15 October 2000, Pages 485-491

    BiBO: Taken from Umemura, Kentaro Miyata, Kiyoshi Kato (2007)

    LT (Lithium Tantalate): Provides bulk Sellmeier equations for congruently grown LiTaO3.
    The Sellmeier equations are from Abedin & Ito 1996
    Double check the crystal vs lab coordinate system! Not sure about this,
    but irrelevant for bulk crystals.

    KTP: Provides bulk Sellmeier equations for flux grown KTP.
    The equations are taken from Kato & Takaoka (2002)
    %}

    wl = wl*1e-6; % convert from nm to um
    
    if crystal=="LN" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oA1 = 4.9048;
        oA2 = 0.11775;
        oA3 = 0.21802;
        oA4 = 0.027153;
        oB1 = 2.2314e-8;
        oB2 = -2.9671e-8;
        oB3 = 2.1429e-8;
    
        eA1 = 5.35583; 
        eA2 = 0.100473;
        eA3 = 0.20692;
        eA4 = 100;
        eA5 = 11.34927;
        eA6 = 1.5334e-2;
        eB1 = 4.629e-7;
        eB2 = 3.862e-8;
        eB3 = -0.89e-8;
        eB4 = 2.657e-5;
    
        if axis==0 | axis==1
            oT = (T - 24.5) * (T + 570.50);
            N = oA1 + (oA2 + oB1 * oT) / (wl^2 - (oA3 + oB2 * oT)^2) + oB3 * oT - oA4 * wl^2;
        end
    
        if axis==2
            eT = (T - 24.5) * (T + 570.82);
            N = eA1 + eB1 * eT + (eA2 + eB2 * eT) / (wl^2 - (eA3 + eB3 * eT)^2) + (eA4 + eB4 * eT) / (wl^2 - eA5^2) - eA6 * wl^2;
        end
    end
    if crystal=="BBO" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if axis==0 | axis==1
            p = [2.7359, 0.01878, 0.01822, 0.01471, 0.0006081, 0.000067406];
        end
        if axis==2
            p = [2.3753, 0.01224, 0.01667, 0.01627, 0.0005716, 0.00006305];
        end
        
    	N =  p(1) + (p(2) / (wl^2 - p(3))) - p(4) * wl^2 + p(5) * wl^4 - p(6) * wl^6;
     
    end
	if crystal == "LT" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oA1 = 4.51224;
        oA2 = 0.0847522;
        oA3 = 0.19876;
        oA4 = -0.0239046;
        oB1 = -9.6649e-9;
        oB2 = 8.815e-8;
        oB3 = 4.25637e-8;
        
        eA1 = 4.52999;
        eA2 = 0.0844313;
        eA3 = 0.20344;
        eA4 = -0.0237909;
        eB1 = 1.72995e-7;
        eB2 = -4.7733e-7;
        eB3 = -8.31467e-8;
        
        if axis==0 | axis==1
            F = (T - 25.0) * (T + 25.0 + 546.0);
            N = oA1 + (oA2 + oB1 * F) / (wl^2 - (oA3 + oB2 * F)^2) + (oB3 * F) + (oA4 * wl^2);
        end
                  
        if axis==2
            F = (T - 25.0) * (T + 25.0 + 546.0);
            N = eA1 + (eA2 + eB1 * F) / (wl^2 - (eA3 + eB2 * F)^2) + (eB3 * F) + (eA4 * wl^2);
        end
    end
    if crystal == "BiBO" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if axis==0
            p = [3.07403,0.03231,0.03163,0.013376];
        end
        if axis==1
            p = [3.16940,0.03717,0.03483,0.01827];
        end
        if axis==2
            p = [3.6545,0.05112,0.03713,0.02261];
        end
        N = p(1) + p(2)/(wl^2-p(3)) - p(4)*wl^2;
    end
    if crystal == "KTP" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x1 = 3.29100;  y1 = 3.45018;  z1 = 4.59423;
        x2 = 0.04140;  y2 = 0.04341;  z2 = 0.06206;
        x3 = 0.03978;  y3 = 0.04597;  z3 = 0.04763;
        x4 = 9.35522;  y4 = 16.98825; z4 = 110.80672;
        x5 = 31.45571; y5 = 39.43799; z5 = 86.12171;
    
        if axis == 0
            N = x1 + x2 / (wl^2 - x3) + x4 / (wl^2 - x5);
        end
    
        if axis == 1
            N =  y1 + y2 / (wl^2 - y3) + y4 / (wl^2 - y5);
        end

        if axis == 2
            N =  z1 + z2 / (wl^2 - z3) + z4 / (wl^2 - z5);
        end  
    end
    if crystal == "KDP" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if axis == 0 | axis == 1 % no idea about these axes actually...
            p = [1.458524, 0.799459, 0.012741, 0.908108];
        end
        if axis == 2
            p = [1.423457, 0.708796, 0.012221, 0.225356];
        end
        N =  p(1) + p(2) / (1 - (p(3)/wl^2)) + p(4) / (1 - (5.0/wl^2)) ;
    end
end
   
    
    
    
    
    
    
    
    
    
    
    
    

