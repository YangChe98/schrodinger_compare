function printgif(printfname,psi_method,time_step,x,dt0)

fig = figure;

for i=1:(time_step+1)
    plot(x,abs(psi_method(:,i)))
    xlim([-10,10])
    ylim([0,1.2])
    drawnow
    frame = getframe(fig);
    im{i} = frame2im(frame);
end

filename = printfname; % Specify the output file name
for idx = 1:(time_step+1)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",dt0);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",dt0);
    end
end
end
