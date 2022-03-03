function [ret] = mcmc_smg(omega,obs,R,nsamp,app_flg,nu,burn,...
                e2i,s2i,as2,bs2,eta_flg,ae2,be2,gibbs_flg,YY_i)

%     cd sparco-1.2
%     run sparcoSetup.m
%     cd '../'
    jit = 1e-3;

    %Set dimensions
    m1 = size(omega,1);
    m2 = size(omega,2);
    sizeX = size(omega);
    n = sum(omega(:));
    [omegax,omegay] = find(omega);
    [nomegax,nomegay] = find(1-omega);
    if ((m1 > 500)||(m2 > 500))
        lrg_flg = true; %large matrices
        gibbs_flg = true;
    else 
        lrg_flg = false;
    end
    
    %Warm start via nuclear-norm completion
    idx = find(omega(:));
    idxc = find(1-omega(:));
    if (isempty(YY_i))
        KK = ones(length(omegax),1);
        msk = sparse(1:length(idx),idx,KK,length(omegax),m1*m2);
        M_op = opMatrix(msk);
        if (~lrg_flg)
            err_tol = 1e-4;
            insweep = 100;
            decfac = 0.5;
            YY_i = NonCVX_MC(obs,M_op,sizeX,1.0,1e-6,zeros(prod(sizeX),1),1,insweep,err_tol,decfac);
        else
            err_tol = 1e-4;
            insweep = 25;
            decfac = 0.25;
            YY_i = NonCVX_MC(obs,M_op,sizeX,1.0,1e-6,zeros(prod(sizeX),1),1,insweep,err_tol,decfac);
        end
    end
%     % plot nuclear norm minimisation results
%     figure,
%     %YY_ii=unpatch(YY_i,m1,m2);
%     imagesc(YY_i) 
%     caxis([-2.5 2.5])
%     colormap(hot(512))
%     colorbar
%     set(gca,'xticklabel',{[]})
%     set(gca,'yticklabel',{[]})
%     title('Matrix Completion by Nuclear-norm Minimization');
%     axis square
    
    
    %Initialize containers
    [UU,DD,VV] = svd(YY_i); %take svd for subspace estimation
    DD = diag(DD);
    UU = UU(:,1:R);
    VV = VV(:,1:R);
    DD = DD(1:R);
    
    eta2 = e2i;
    sig2 = s2i;
    YY_i(idx) = obs;
    YY = YY_i;
    
    ret = struct; %struct for MCMC samples
    ret.YY = zeros(nsamp,m1,m2);
    ret.XX = zeros(nsamp,m1,m2);
    ret.UU = zeros(nsamp,m1,R);
    ret.VV = zeros(nsamp,m2,R);
    ret.DD = zeros(nsamp,R);
    ret.eta2 = zeros(nsamp,1);
    ret.sig2 = zeros(nsamp,1);
    
    %Gibbs sampling
    for (i = 1:nsamp)
        
        disp(['MCMC sample: ' num2str(i) '/' num2str(nsamp) '...'])
        
        %Impute missing entries
        XX = UU*(diag(DD)*VV');
        ret.XX(i,:,:) = XX;
        YY(idxc) = XX(idxc) + sqrt(eta2)*randn(length(idxc),1);
        ret.YY(i,:,:) = YY;
        
        %Gibbs sampling on row subspace
%         disp('Row space sampling...')
        if (~gibbs_flg)
            UU = rmf_matrix( YY * (VV * diag(DD)) / eta2 );
        else
            UU = rmf_matrix_gibbs( YY * (VV * diag(DD)) / eta2, UU, -1);
        end
        
        while (norm(UU'*UU-eye(R)) > 1e-5)
            disp('Numerical issues in matrix Gibbs. Jittering...')
            UU = reshape(ret.UU(i-1,:,:),[m1 R]);
            MM = YY * (VV * diag(DD)) / eta2;
            MM = MM + jit*max(abs(MM(:)))*randn(size(MM,1),size(MM,2));
            UU = rmf_matrix_gibbs( MM, UU, -1);
        end
        ret.UU(i,:,:) = UU;
        
        %Gibbs sampling on column subspace
%         disp('Column space sampling...')
        if (~gibbs_flg)
            VV = rmf_matrix( YY' * (UU * diag(DD)) / eta2 );
        else
            VV = rmf_matrix_gibbs( YY' * (UU * diag(DD)) / eta2, VV, -1);
        end
        
        while (norm(VV'*VV-eye(R)) > 1e-5)
            disp('Numerical issues in matrix Gibbs. Jittering...')
            VV = reshape(ret.VV(i-1,:,:),[m1 R]);
            MM = YY' * (UU * diag(DD)) / eta2;
            MM = MM + jit*max(abs(MM(:)))*randn(size(MM,1),size(MM,2));
            VV = rmf_matrix_gibbs( MM, VV, -1);
        end
        ret.VV(i,:,:) = VV;
        
        %Singular values
        mu = sig2*diag(UU'*(YY*VV))./(eta2+sig2);
        del = sqrt(eta2*sig2/(eta2+sig2));
        if (~app_flg)
            %Independence sampler
            DD = rql(DD,mu,del,false,nu,burn);
        else
            %Gaussian approximation on RG
            DD = zeros(R,1);
            while ~(all(DD>0.0))
                DD = del*randn(R,1) + mu;
            end
        end
%         % Testing IG priors in singular values, with subspace sampling
%         DD = 1.0 ./ gamrnd( 0.01+(m1+m2)/2.0, 0.01+DD );
        ret.DD(i,:) = DD;
        
        %Gibbs sampling on sigma^2 and eta^2
        sig2 = 1.0/gamrnd( as2+R/2, 1.0/(bs2+sum(DD.^2)/2) );
        if (~eta_flg)
            eta2 = 1.0/gamrnd( ae2+m1*m2/2, 1.0/( be2+sum((YY(:)-XX(:)).^2)/2 ) );
        end
        ret.sig2(i,1) = sig2;
        ret.eta2(i,1) = eta2;
        
    end
    
%     ret    

end
