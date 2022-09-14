clear; close all;
load force_gpu.out; % GPU (GPUMD) analytical
load force_analytical.out; % CPU analytical
load force_finite_difference.out; % CPU finite difference
load descriptor.out;
load virial.out; % new: correct for heat current
load force_analytical_ref.out;
load force_finite_difference_ref.out;
load descriptor_ref.out;
load virial_ref.out; % old: wrong for heat current
max_force=max(max(abs(force_gpu)));

% should be about 1e-14 or zero
figure;
subplot(2,3,1)
plot(descriptor-descriptor_ref,'o-');hold on;
xlabel('components');
ylabel('diff-descriptor');
set(gca,'fontsize',15);

% should be about 1e-14 or zero
subplot(2,3,2)
plot(force_analytical-force_analytical_ref,'o-');hold on;
xlabel('components');
ylabel('diff-analytical');
set(gca,'fontsize',15);

% should be about 1e-8 or zero
subplot(2,3,3)
plot(force_finite_difference-force_finite_difference_ref,'o-');hold on;
xlabel('components');
ylabel('diff-FD');
set(gca,'fontsize',15);

% should be about 1e-14
subplot(2,3,4)
plot(mean(virial)-mean(virial_ref),'o-');hold on;
xlabel('components');
ylabel('diff-virial');
set(gca,'fontsize',15);

% should be of the order of 1e-6 (used float32 in GPU)
subplot(2,3,5)
plot((force_gpu-force_analytical)/max_force);
xlabel('force components');
ylabel('GPU (float) - CPU (double) (eV/A)');
set(gca,'fontsize',15);

% should be of the order of 1e-8
subplot(2,3,6)
plot((force_finite_difference-force_analytical)/max_force);
xlabel('force components');
ylabel('FD - analytical (eV/A)');
set(gca,'fontsize',15);
