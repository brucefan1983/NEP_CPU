clear; close all;
load force_gpu.out; % GPU (GPUMD) analytical
load force_analytical.out; % CPU analytical
load force_finite_difference.out; % CPU finite difference
load descriptor.out;
load force_analytical_ref.out;
load force_finite_difference_ref.out;
load descriptor_ref.out;

figure;
plot(descriptor-descriptor_ref,'o-');hold on;
xlabel('components');
ylabel('descriptor - descriptor_ref');
set(gca,'fontsize',15);

figure;
plot(force_analytical-force_analytical_ref,'o-');hold on;
xlabel('components');
ylabel('force_analytical - force_analytical_ref');
set(gca,'fontsize',15);

figure;
plot(force_finite_difference-force_finite_difference_ref,'o-');hold on;
xlabel('components');
ylabel('force_finite_difference - force_finite_difference_ref');
set(gca,'fontsize',15);

% The difference should be of the order of 1.0e-5 (used float32 in GPU)
figure;
plot(force_gpu-force_analytical);
xlabel('force components');
ylabel('GPU (float) - CPU (double) (eV/A)');
set(gca,'fontsize',15);

figure;
plot(force_finite_difference-force_analytical);
xlabel('force components');
ylabel('finite-difference - analytical (eV/A)');
set(gca,'fontsize',15);


