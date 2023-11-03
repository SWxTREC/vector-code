function [avg_cos_frac, avg_accom, avg_alpha_n, avg_sigma_t] = sphere_incident_angle_weighting(r,material)
%r is radius of the sphere
%break-up the surface of the sphere into incident angle regions
%returns weighted average GSI parameters based on laboratory surface scattering experiments


if material == 'SiO2'
    %range of SiO2 parameters
    incident_angles_data = [30, 45, 60]; %degrees
    alpha_n_data = [0.98, 0.99, 0.98];
    sigma_t_data = [0.49, 0.81, 0.83];
    cos_frac_data = [0.97, 0.9, 0.79];
    alpha_data = [0.71, 0.61, 0.44];

elseif material == 'aluminum'
    %range of aluminum parameters
    incident_angles_data = [30, 60]; %degrees
    alpha_n_data = [0.99, 0.98];
    sigma_t_data = [0.59, 0.85];
    cos_frac_data = [0.98, 0.88];
    alpha_data = [0.72, 0.55];

elseif material == 'Teflon'
    %range of teflon parameters
    incident_angles_data = [30, 45, 60]; %degrees
    alpha_n_data = [0.99, 0.99, 0.98];
    sigma_t_data = [0.69, 0.81, 0.84];
    cos_frac_data = [0.96, 0.87, 0.76];
    alpha_data = [0.62, 0.50, 0.30];

elseif material == 'FR4'
    %range of FR4 parameters
    incident_angles_data = [30, 45, 60]; %degrees
    alpha_n_data = [0.99, 0.99, 0.98];
    sigma_t_data = [0.59, 0.79, 0.85];
    cos_frac_data = [0.97, 0.93, 0.85];
    alpha_data = [0.68, 0.63, 0.52];

end

inc_angle_res = 2.5; %degrees
num_sections = (90/inc_angle_res)+1;

cos_frac_vals = zeros(1,num_sections);
accom_vals = zeros(1,num_sections);
alpha_n_vals = zeros(1,num_sections);
sigma_t_vals = zeros(1,num_sections);
across_vals = zeros(1,num_sections);
ascale_vals = zeros(1,num_sections);

CD_parts = zeros(1,num_sections);

for i = 1:num_sections
    inc_angle_bound1 = (i-1)*inc_angle_res;
    inc_angle_bound2 = i*inc_angle_res;
    if i == num_sections
        mid_angle = 90;
    else
        mid_angle = inc_angle_bound1 + ((inc_angle_bound2 - inc_angle_bound1)/2);
    end
    cos_frac = interp1(incident_angles_data, cos_frac_data, mid_angle, 'linear', 'extrap');
    if cos_frac > 1
        cos_frac = 1;
    elseif cos_frac < 0
        cos_frac = 0;
    end
    alpha_n = interp1(incident_angles_data, alpha_n_data, mid_angle, 'linear', 'extrap');
    if alpha_n > 1
        alpha_n = 1;
    elseif alpha_n < 0
        alpha_n = 0;
    end
    sigma_n = 1 - sqrt(1 - alpha_n);
    sigma_t = interp1(incident_angles_data, sigma_t_data, mid_angle, 'linear', 'extrap');
    if sigma_t > 1
        sigma_t = 1;
    elseif sigma_t < 0
        sigma_t = 0;
    end
    accom = interp1(incident_angles_data, alpha_data, mid_angle, 'linear', 'extrap');
    if accom > 1
        accom = 1;
    elseif accom < 0
        accom = 0;
    end
    r1 = r*sin(deg2rad(inc_angle_bound1));
    h1 = r*cos(deg2rad(inc_angle_bound1));
    r2 = r*sin(deg2rad(inc_angle_bound2));
    h2 = r*cos(deg2rad(inc_angle_bound2));
    if i == num_sections %back hemisphere
        SA = 2*pi*r^2; %surface area of half of a sphere
        a_cross = pi*r^2; %area of circle with radius r
    else %spherical segment (ring)
        SA = 2*pi*r*(h1-h2); %surface area
        a_cross = pi*(r2^2 - r1^2);
    end

    cos_frac_vals(1,i) = cos_frac;
    accom_vals(1,i) = accom;
    alpha_n_vals(1,i) = alpha_n;
    sigma_t_vals(1,i) = sigma_t;
    across_vals(1,i) = a_cross;
    ascale_vals(1,i) = SA;

end

%weighted average parameters, weighted by surface area for only the segments on the front half of the sphere
avg_cos_frac = sum(ascale_vals(1:end-1).*cos_frac_vals(1:end-1))/sum(ascale_vals(1:end-1));
avg_accom = sum(ascale_vals(1:end-1).*accom_vals(1:end-1))/sum(ascale_vals(1:end-1));
avg_alpha_n = sum(ascale_vals(1:end-1).*alpha_n_vals(1:end-1))/sum(ascale_vals(1:end-1));
avg_sigma_t = sum(ascale_vals(1:end-1).*sigma_t_vals(1:end-1))/sum(ascale_vals(1:end-1));






















