% Conceptual figure

clear all
close all
clc


t = 0:0.001:10;

% Sender: sum of sinusoids
s = 0.2*sin(t) + 0.2*sin(1.5*t-pi/4) + 0.3*sin(0.75*t+pi) + ...
    0.1*sin(5*t) + 0.3*sin(0.5*t + pi/2) + 0.05*sin(10*t);

% Sender activity
figure()
plot(t, s)
box off
ax1 = gca;
set(gca, 'TickDir', 'out')
pbaspect([3 1 1])
title('Sender')
ylim([-0.7 1])


% Linear receiver
r = 0.5*s;

figure()
plot(t, r)
box off
ax2 = gca;
set(gca, 'TickDir', 'out')
pbaspect([3 1 1])
title('Linear receiver')
ylim([-0.7 1])


% Non-linear receiver
r = -0.8*s.*(s-0.5).*(s-1);

figure()
plot(t, r)
box off
ax3 = gca;
set(gca, 'TickDir', 'out')
pbaspect([3 1 1])
title('Non-linear receiver')
ylim([-0.7 1])


% Noisy receiver
tau = 10;
c = 0.05;
dtovertau = 1;
dt = 0.001;
ou = NaN(1,length(s));
ou(1) = s(1);
for i = 2:length(t)
    ou(i) = ou(i-1) - (1/tau)*ou(i-1)*dt + sqrt(c)*sqrt(dt)*randn;
end
r = s + ou;

figure()
plot(t, r)
box off
ax4 = gca;
set(gca, 'TickDir', 'out')
pbaspect([3 1 1])
title('Noisy receiver')
ylim([-0.7 1])

% Link y-axis of plots
linkaxes([ax1 ax2 ax3 ax4],'xy')


