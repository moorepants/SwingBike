% File: sb_eiginfo.m
% Date: October 14, 2008
% Author: Jason Moore
% Description: Plots the phasors for the swing bike model.
function sb_eiginfo(speed,v,eigval,eigvec)
% speed = speed at which to check the eigenmodes
% v = vector of speeds
% eigval = eignenvalues at each speed
% eigvec = eigenvectors at each speed

% check to make sure the selected speed is within the calculated range
if speed > max(v)
    s=sprintf('The selected speed too high, choose a speed below %1.2f m/s.',max(v));
    disp(s)
else
    % find the speed nearest the selected speed
    a=speed*ones(size(v))-v;
    [y,index]=min(abs(a));
    mag = abs(eigval(:,index));
    pair = 0;
    for i = 1:length(mag)
        for j = (i+1):length(mag)
            test = mag(i)-mag(j);
            if test == 0
                pair = pair + 1;
            end
        end
    end
    modes=length(eigval(:,index))-pair;
    s=sprintf('There are %d modes',modes);
    disp(s)
    labels = {'q4 (roll angle)       '
              'q7 (rear steer angle) '
              'q8 (front steer angle)'
              'u4 (roll rate)        '
              'u7 (rear steer rate)  '
              'u8 (front steer rate) '};
    for i = 1:min(size(eigval))
        s1=sprintf('Speed = %1.2f m/s, Eigenvalue = %4.4f + %4.4fj ',v(index),real(eigval(i,index)),imag(eigval(i,index)));
        disp(s1)
        if real(eigval(i,index)) < 0
            s2='Stable';
            disp(s2)
        elseif real(eigval(i,index)) == 0
            s2='Marginally Stable';
            disp(s2)
        else
            s2='Unstable';
            disp(s2)
        end
        disp('Eigenvector')
        for j = 1:min(size(eigval))
            s = sprintf('%s %4.4f + %4.4fj',labels{j},real(eigvec(j,i,index)),imag(eigvec(j,i,index)));
            disp(s)
        end
        x = real(eigvec(:,i,index));
        y = imag(eigvec(:,i,index));
        figure(i)
        color = colormap(jet(length(eigval(:,1))));
        subplot(211)
        hold on
        plot([-1;1],[0,0],'k-')
        plot([0,0],[-1,1],'k-')
        for k = 1:min(size(eigval))
            plot([0;x(k)],[0;y(k)],'-o','color',color(k,:),'linewidth',3,'markeredgecolor',color(k,:),'markerfacecolor',color(k,:),'markersize',8)
        end
        %plot(0,0,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],'markersize',8)
        legend('real axis','imaginary axis',labels{1},labels{2},labels{3},labels{4},labels{5},labels{6})
        title([s1 s2])
        axis square
        box on
        hold off
        subplot(212)
        t = 0:0.01:10;
        hold on
        for k = 1:min(size(eigval))
            mag = abs(eigvec(k,i,index));
            phi = angle(eigvec(k,i,index));
            q = mag*exp(real(eigval(i,index)).*t).*cos(imag(eigval(i,index)).*t+phi);
            plot(t,q,'color',color(k,:),'linewidth',2)
        end
        axis([0 max(t) -1 1])
        box on
    end
end
