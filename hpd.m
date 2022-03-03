function [lw,up] = hpd(samps,burn,perc)
    
    try
        %matrices
        lw = zeros(size(samps,2),size(samps,3));
        up = zeros(size(samps,2),size(samps,3));
        nsamp = size(samps,1);
        for (jj = 1:size(samps,3))
    %         disp(jj)
            qt = quantile(samps((burn+1):nsamp,:,jj),[(1-perc)/2 1-(1-perc)/2]);
            lw(:,jj) = qt(1,:);
            up(:,jj) = qt(2,:);
        end
    catch 
        %vectors/scalars
        lw = zeros(size(samps,2),1);
        up = zeros(size(samps,2),1);
        nsamp = size(samps,1);
        for (jj = 1:size(samps,2))
    %         disp(jj)
            qt = quantile(samps((burn+1):nsamp,jj),[(1-perc)/2 1-(1-perc)/2]);
            lw(jj,1) = qt(1);
            up(jj,1) = qt(2);
        end
end